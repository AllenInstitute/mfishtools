#' Build and plot dendrogram from gene panel
#' 
#' Build and plot a dendrogram using correlation-based average linkage hierarchical clustering
#'   and only using a specified set of genes.  The output is the expected accuracy of mapping 
#'   to each node in the tree, which gives an idea of the best-case expected results for 
#'   mFISH analysis.
#'
#' @param dend dendrogram for mapping.  Ignored if medianDat is passed
#' @param refDat normalized data of the REFERENCE data set.  Ignored if medianDat is passed
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param medianDat representative value for each leaf and node.  If not entered, it is calculated
#' @param requiredGenes minimum number of genes required to be expressed in a cluster (column
#'   of medianDat) for the cluster to be included (default=2)
#' @param clusters  cluster calls for each cell
#' @param genesToMap which genes to include in the correlation mapping
#' @param plotdendro should the dendrogram be plotted (default = TRUE)
#' @param returnDendro should the dendrogram be returned (default = TRUE)
#' @param mar margins (for use with par)
#' @param use,... additional parameters for cor
#'
#' @return a list where the first entry is the resulting tree and the second entry is the 
#'   fraction of cells correctly mapping to each node using the inputted gene panel.
#'
buildTreeFromGenePanel <- function(dend = NA, refDat = NA, mapDat = refDat, medianDat = NA, 
  requiredGenes = 2, clusters = NA, genesToMap = rownames(mapDat), plotdendro = TRUE, 
  returndendro = TRUE, mar = c(12, 5, 5, 5), use = "p", ...) {
  
  library(dendextend)
  
  # Calculate the median, if needed.
  if (is.na(medianDat[1])) {
    names(clusters) = colnames(refDat)
    medianDat = do.call("cbind", tapply(names(clusters), clusters, function(x) rowMedians(refDat[, 
      x])))
    medianDat = leafToNodeMedians(dend, medianDat)
  }
  gns = intersect(genesToMap, intersect(rownames(mapDat), rownames(medianDat)))
  
  # Subset the data to relevant genes and clusters
  medianDat <- medianDat[gns, ]
  mapDat <- mapDat[gns, ]
  medianDat <- medianDat[, colSums(medianDat > 0) >= requiredGenes]
  kpDat <- (colSums(mapDat > 0) >= requiredGenes) & (is.element(clusters, colnames(medianDat)))
  mapDat <- mapDat[, kpDat]
  
  # Perform the correlation mapping
  facsCor <- corTreeMapping(medianDat = medianDat, mapDat = mapDat, use = use, ...)
  facsCl <- colnames(facsCor)[apply(facsCor, 1, which.max)]
  
  # Build a new tree based on mapping
  sCore <- function(x, use, ...) return(as.dist(1 - cor(x, use = use, ...)))
  dend <- getDend(medianDat, sCore, use = use, ...)
  
  # Which leaves have which nodes?
  has_any_labels <- function(sub_dend, the_labels) any(labels(sub_dend) %in% the_labels)
  node_labels <- NULL
  for (lab in labels(dend)) node_labels <- cbind(node_labels, noded_with_condition(dend, 
    has_any_labels, the_labels = lab))
  rownames(node_labels) <- get_nodes_attr(dend, "label")
  colnames(node_labels) <- labels(dend)
  
  # Which clusters agree at the node level?
  clTmp = as.character(clusters[kpDat])
  agreeNodes = apply(cbind(facsCl, clTmp), 1, function(lab, node_labels) {
    rowSums(node_labels[, lab]) == 2
  }, node_labels)
  colnames(agreeNodes) = clTmp
  
  # Which clusters are in each nodes?
  isInNodes = t(apply(node_labels, 1, function(node, cl, dend) {
    is.element(cl, labels(dend)[node])
  }, clTmp, dend))
  colnames(isInNodes) = clTmp
  
  # For each node, plot the fraction of cells that match if desired?
  fracAgree = rowSums(agreeNodes)/rowSums(isInNodes)
  if (plotdendro) {
    par(mar = mar)
    dend %>% set("nodes_cex", 0) %>% set("branches_col", "grey") %>% plot
    text(get_nodes_xy(dend)[, 1], get_nodes_xy(dend)[, 2], round(fracAgree * 100))
  }
  
  # Return the results (if desired)
  if (returndendro) 
    return(list(dend, fracAgree))
  
}


#' Build a dendrogram from gene panel
#' 
#' Build a dendrogram from an inputted data matrix.
#'
#' @param dat matrix of values (e.g., genes x clusters) for calculating the dendrogram
#' @param distFun function for calculating distance matrix (default is correlation-based)
#' @param ... additional variables for distFun
#'
#' @return dendrogram 
#'
getDend <- function(dat, distFun = function(x) return(as.dist(1 - cor(x))), ...) {
  distCor = distFun(dat, ...)
  distCor[is.na(distCor)] = max(distCor, na.rm = TRUE) * 1.2
  # Avoid crashing by setting NA values to hang off the side of the tree.
  avgClust = hclust(distCor, method = "average")
  dend = as.dendrogram(avgClust)
  dend = labelDend(dend)[[1]]
  return(dend)
}


#' Label dendrogram nodes
#' 
#' Add numeric node labels to a dendrogram.
#'
#' @param dend dendrogram object
#' @param distFun starting numeric node value (default=1)
#'
#' @return a list where the first item is the new dendrogram object and the second item is the final
#'   numeric node value 
#'
labelDend <- function(dend, n = 1) {
  if (is.null(attr(dend, "label"))) {
    attr(dend, "label") = paste0("n", n)
    n = n + 1
  }
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      tmp = labelDend(dend[[i]], n)
      dend[[i]] = tmp[[1]]
      n = tmp[[2]]
    }
  }
  return(list(dend, n))
}


#' Summarize matrix
#' 
#' Groups columns in a matrix by a specified group vector and summarizes using a specificed function.
#'   Optionally binarizes the matrix using a specified cutoff parameter.  This is a wrapper for tapply.
#'
#' @param mat matrix where the columns (e.g., samples) are going to be grouped
#' @param group vector of length dim(mat)[2] corresponding to the groups
#' @param scale either 'none' (default),'row', or 'column'
#' @param scaleQuantile what quantile of value should be set as 1 (default=1)
#' @param binarize should the data be binarized? (default=FALSE)
#' @param binMin minimum ON value for the binarized matrix (ignored if binarize=FALSE)
#' @param summaryFunction function (or function name) to be used for summarization
#' @param ... additional parameters for summaryFunction
#'
#' @return matrix of summarized values
#'
summarizeMatrix <- function(mat, group, scale = "none", scaleQuantile = 1, binarize = FALSE, 
  binMin = 0.5, summaryFunction = median, ...) {
  
  # Make sure the names match up
  if (is.null(colnames(mat))) 
    colnames(mat) = names(group)
  if (is.null(colnames(mat))) 
    colnames(mat) = 1:length(group)
  names(group) <- colnames(mat)
  
  # Calculate the summary
  summaryFunction <- match.fun(summaryFunction)
  runFunction <- function(x, ...) {
    if (length(x) > 1) 
      return(apply(mat[, x], 1, summaryFunction, ...))
    return(mat[, x])
  }
  summarizedMat <- do.call("cbind", tapply(names(group), group, runFunction, ...))
  if (is.factor(group)) 
    summarizedMat = summarizedMat[, levels(group)]
  rownames(summarizedMat) = rownames(mat)
  
  # Scale the data if desired
  if (substr(scale, 1, 1) == "r") 
    for (i in 1:dim(summarizedMat)[1]) summarizedMat[i, ] = summarizedMat[i, ]/max(1e-06, 
      quantile(summarizedMat[i, ], probs = scaleQuantile))
  if (substr(scale, 1, 1) == "c") 
    for (i in 1:dim(summarizedMat)[2]) summarizedMat[, i] = summarizedMat[, i]/max(1e-06, 
      quantile(summarizedMat[, i], probs = scaleQuantile))
  
  # Binarize the data if desired
  if (binarize) {
    summarizedMat = summarizedMat > binMin
    summarizedMat = summarizedMat + 1 - 1
  }
  
  return(summarizedMat)
}


#' Scale mFISH data and map to RNA-seq reference
#' 
#' This function is a wrapper for several other functions which aim to scale mFISH data to 
#'   more closely match RNA-seq data and then map the mFISH data to the closest reference
#'   classes.  There are several parameters allowing flexability in filtering and analysis.
#'
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param refSummaryDat normalized summary data of the REFERENCE data set (e.g., what to map against)
#' @param genesToMap which genes to include in the mapping (calculated in not entered)
#' @param mappingFunction which function to use for mapping (default is cellToClusterMapping_byCor)
#'   The function must include at least two parameters with the first one being mapped data and the 
#'   second data the reference.  Additional parameters are okay.  Output must be a data frame where 
#'   the first value is a mapped class.  Additional columns are okay and will be returned)
#' @param transform function for transformation of the data (default in none)
#' @param noiselevel scalar value below which all values are set to 0 (default is 0)
#' @param scaleFunction which function to use for scaling mapDat to refSummaryDat (default is setting
#'   90th quantile of mapDat to max of refSummaryDat and truncating higher mapDat values)
#' @param scaleXY should x and y coordinates be scaled from 0-1 within experiments (default = TRUE)
#' @param metadata a data frame of possible metadata (additional columns are okay and ignored):
#' \describe{
#'   \item{area}{a vector of cell areas for normalization}
#'   \item{experiment}{a vector indicating if multiple experiments should be scaled separately}
#'   \item{x,y}{x (e.g., parallel to layer) and y (e.g., across cortical layers) coordinates in tissue}
#' }
#' @param ... additional parameters for passthrough into other functions
#'
#' @return a list with the following entrees:
#' \describe{
#'   \item{mapDat}{mapDat data matrix is passed through}
#'   \item{scaleDat}{scaled mapDat data matrix}
#'   \item{mappingResults}{Results of the mapping and associated confidence values (if any)}
#'   \item{metadata=metadata}{metadata is passed through unchanged}
#'   \item{scaledX/Y}{scaled x and y coordinates (or unscaled if scaling was not performed)}
#' }
#'
fishScaleAndMap <- function(mapDat, refSummaryDat, genesToMap = NULL, mappingFunction = cellToClusterMapping_byCor, 
  transform = function(x) x, noiselevel = 0, scaleFunction = quantileTruncate, scaleXY = TRUE, 
  metadata = data.frame(experiment = rep("all", dim(mapDat)[2])), ...) {
  
  # Setup
  mappingFunction <- match.fun(mappingFunction)
  scaleFunction <- match.fun(scaleFunction)
  transform <- match.fun(transform)
  if (is.null(genesToMap)) 
    genesToMap <- colnames(mapDat)
  genesToMap <- intersect(genesToMap, rownames(refSummaryDat))
  params <- colnames(metadata)
  refSummaryDat <- refSummaryDat[genesToMap, ]
  mapDat <- mapDat[genesToMap, ]
  
  # Transform the data to be mapped
  scaleDat <- as.matrix(mapDat[genesToMap, ])
  scaleDat[scaleDat < noiselevel] = 0  # Set values less than noiselevel to 0
  if (is.element("area", params)) {
    # Account for spot area in gene expression calculation
    scaleDat <- t(t(scaleDat)/metadata$area) * mean(metadata$area)
  }
  scaleDat <- transform(scaleDat)
  
  # Scale to the reference data
  for (ex in unique(metadata$experiment)) {
    isExp = metadata$experiment == ex
    for (g in genesToMap) scaleDat[g, isExp] <- scaleFunction(scaleDat[g, isExp], 
      maxVal = max(refSummaryDat[g, ]), ...)
  }
  
  # Map the map data to the reference data
  mappingResults <- mappingFunction(refSummaryDat, scaleDat, ...)
  
  # Scale x and y coordinates to (0,1) within experiment, if desired
  if (scaleXY) {
    for (ex in unique(metadata$experiment)) {
      isExp = metadata$experiment == ex
      metadata$x[isExp] = metadata$x[isExp] - min(metadata$x[isExp])
      metadata$x[isExp] = metadata$x[isExp]/max(metadata$x[isExp])
      metadata$y[isExp] = metadata$y[isExp] - min(metadata$y[isExp])
      metadata$y[isExp] = metadata$y[isExp]/max(metadata$y[isExp])
    }
  }
  
  # Return the results
  out = list(mapDat = mapDat, scaleDat = scaleDat, mappingResults = mappingResults, 
    metadata = metadata, scaledX = metadata$x, scaledY = metadata$y)
}


#' Filter cells
#' 
#' Remove a select set of cells from all elements of a fishScaleAndMap output
#'
#' @param datIn a fishScaleAndMap output list
#' @param kpSamp a vector of samples to keep (of the same length as the number of cells in datIn)
#'
#' @return filtered fishScaleAndMap object
#'
filterCells <- function(datIn, kpSamp) {
  datOut <- datIn
  datOut$mapDat <- datOut$mapDat[, kpSamp]
  datOut$scaleDat <- datOut$scaleDat[, kpSamp]
  datOut$mappingResults <- datOut$mappingResults[kpSamp, ]
  datOut$metadata <- datOut$metadata[kpSamp, ]
  datOut$scaledX <- datOut$scaledX[kpSamp]
  datOut$scaledY <- datOut$scaledY[kpSamp]
  return(datOut)
}


#' Plot distributions
#' 
#' Plot the distributions of cells across the tissue with overlaying color information.  This is
#'   a wrapper function for plot
#'
#' @param datIn  a fishScaleAndMap output list
#' @param group a character vector (or factor) indicating how to split the data (e.g., cluster 
#'   call) or a metadata/mappingResults column name
#' @param groups a character vector of groups to show (default is levels of group)
#' @param colors a character vector (or factor) indicating how to color the plots (e.g., layer 
#'   or gene expression) or a metadata/mappingResults column name (default is all black)
#' @param colormap function to use for the colormap for the data (default gray.colors)
#' @param xlim,ylim for plot, but will be calculated if not entered
#' @param pch,cex for plot.  Can be single values or vectors
#' @param main,xlab,ylab,... other parameters for plot (must be single values)
#'
#' @return Only returns if there is an error
#'
plotDistributions <- function(datIn, group, groups = NULL, colors = rep("black", dim(datIn$mapDat)[2]), 
  colormap = gray.colors, pch = 19, cex = 1.5, xlim = NULL, ylim = NULL, main = "", 
  xlab = "", ylab = "", ...) {
  
  colormap = match.fun(colormap)
  meta = cbind(datIn$metadata, datIn$mappingResults)
  if (length(group) == 1) {
    if (is.element(group, colnames(meta))) {
      group = as.factor(meta[, group])
      if (is.null(groups)) 
        groups = levels(group)
    } else {
      return(paste(group, "is not an available column name for division."))
    }
  }
  if (length(colors) == 1) {
    if (is.element(colors, colnames(meta))) {
      colors = as.numeric(as.factor(meta[, colors]))
      colors = colormap(length(unique(colors)))[colors]
    } else {
      return(paste(colors, "is not an available column name for coloring."))
    }
  } else {
    colors = as.numeric(as.factor(colors))
    colors = colormap(length(unique(colors)))[colors]
  }
  
  if (is.null(xlim)) 
    xlim = range(datIn$scaledX)
  if (is.null(ylim)) 
    ylim = range(-datIn$scaledY)
  
  # Make the plot!
  par(mfrow = c(1, length(groups)))
  for (gp in groups) {
    kp = group == gp
    pch2 = pch
    if (length(pch) > 1) 
      pch2 = pch[kp]
    cex2 = cex
    if (length(cex) > 1) 
      cex2 = cex[kp]
    
    plot(datIn$scaledX[kp], -datIn$scaledY[kp], pch = pch2, col = colors[kp], xlim = xlim, 
      ylim = ylim, main = paste(main, gp), xlab = xlab, ylab = ylab, cex = cex2, 
      ...)
  }
}


#' Plot heatmap
#' 
#' Plot the heatmap of cells ordering by a specified order.  This is a wrapper for heatmap.2
#'
#' @param datIn  a fishScaleAndMap output list
#' @param group a character vector (or factor) indicating how to order the heatmap (e.g., cluster 
#'   call) or a metadata/mappingResults column name
#' @param groups a character vector of groups to show (default is levels of group)
#' @param grouplab label for the grouping in the heatmap (default is 'Grouping' or the value for group)
#' @param useScaled plot the scaled (TRUE) or unscaled (FALSE; default) values
#' @param capValue values above capValue will be capped at capValue (default is none)
#' @param colormap set of values to use for the colormap for the data (default heat_colors)
#' @param Rowv,Colv,dendrogram,trace,margins,rowsep,key,... other parameters for heatmap.2 (some 
#'   default values are different)
#'
#' @return Only returns if there is an error
#'
plotHeatmap <- function(datIn, group, groups = NULL, grouplab = "Grouping", useScaled = FALSE, 
  capValue = Inf, colormap = grey.colors(1000), pch = 19, xlim = NULL, ylim = NULL, Rowv = FALSE, 
  Colv = FALSE, dendrogram = "none", trace = "none", margins = c(3, 10), rowsep = 9, 
  key = FALSE, ...) {
  
  library(gplots)
  
  if (useScaled) {
    plotDat <- datIn$scaleDat
  } else {
    plotDat <- datIn$mapDat
  }
  plotDat = pmin(plotDat, capValue)
  
  meta = cbind(datIn$metadata, datIn$mappingResults)
  if (length(group) == 1) {
    if (is.element(group, colnames(meta))) {
      if (grouplab == "Grouping") 
        grouplab = group
      group = as.factor(meta[, group])
      if (is.null(groups)) 
        groups = levels(group)
    } else {
      return(paste(group, "is not an available column name for division."))
    }
  }
  # Update the levels if needed
  levels(group) = c(groups, setdiff(levels(group), groups))
  
  # Make the plot!
  plotDat = rbind(plotDat,as.numeric(group)*cap/length(groups))
  plotDat = plotDat[,order(group, -colSums(plotDat)) ]
  rownames(plotDat) = c(rownames(plotDat)[1:(dim(plotDat)[1]-1)],grouplab)
  heatmap.2(plotDat, Rowv = Rowv, Colv = Colv, dendrogram = dendrogram, trace = trace, 
    margins = margins, rowsep = rowsep, key = key, col = colormap, ...)
}


#' Return top mapped correlation-based cluster and confidence
#' 
#' Primary function for doing correlation-based mapping to cluster medians and also reporting the
#'   correlations and confidences.  This is wrapper for getTopMatch and corTreeMapping.
#'
#' @param medianDat representative value for each leaf and node.  If not entered, it is calculated
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param refDat normalized data of the REFERENCE data set.  Ignored if medianDat is passed
#' @param clusters  cluster calls for each cell.  Ignored if medianDat is passed
#' @param genesToMap which genes to include in the correlation mapping
#' @param use additional parameter for cor (use='p' as default)
#' @param method additional parameter for cor (method='p' as default)
#' @param ... not used
#'
#' @return data frame with the top match and associated correlation
#'
cellToClusterMapping_byCor <- function(medianDat, mapDat, refDat = NA, clusters = NA, 
  genesToMap = rownames(mapDat), use = "p", method = "p", ...) {
  corVar <- corTreeMapping(medianDat = medianDat, mapDat = mapDat, refDat = refDat, 
    clusters = clusters, genesToMap = genesToMap, use = use, method = method)
  corMatch <- getTopMatch(corVar)
  colnames(corMatch) <- c("Class", "Correlation")
  
  dex <- apply(corVar, 1, function(x) return(diff(sort(-x)[1:2])))
  corMatch$DifferenceBetweenTopTwoCorrelations = dex
  return(corMatch)
}


#' Quantile normalize, truncate, and scale
#' 
#' Quantile normalize, truncate, and scale a numeric vector (e.g. mFISH data from one gene)
#'
#' @param x input data vector
#' @param qprob probs value to result from quantile (default=0.9)
#' @param maxVal max value for scaling (default=1)
#' @param truncate should data above the qprob threshold be truncated (default=yes)
#' @param ... not used
#'
#' @return scaled vector
#' }
#'
quantileTruncate <- function(x, qprob = 0.9, maxVal = 1, truncate = TRUE, ...) {
  qs = quantile(x[x > 0], probs = qprob, na.rm = TRUE)
  if (is.na(qs)) 
    return(x)
  if (truncate) 
    x[x > qs] = qs
  return(x * maxVal/qs)
}


