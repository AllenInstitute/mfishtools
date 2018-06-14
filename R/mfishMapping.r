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
#' @param binarize should the data be binarized? (default=FALSE)
#' @param binMin minimum ON value for the binarized matrix (ignored if binarize=FALSE)
#' @param summaryFunction function (or function name) to be used for summarization
#' @param ... additional parameters for summaryFunction
#'
#' @return matrix of summarized values
#'
summarizeMatrix <- function(mat, group, binarize = FALSE, binMin = 0.5, summaryFunction = median, 
  ...) {
  
  # Make sure the names match up
  if (is.null(colnames(mat))) 
    colnames(mat) = names(group)
  if (is.null(colnames(mat))) 
    colnames(mat) = 1:length(group)
  names(group) <- colnames(mat)
  
  # Calculate the summary
  summaryFunction <- match.fun(summaryFunction)
  summarizedMat <- do.call("cbind", tapply(names(group), group, function(x) apply(mat[, 
    x], 1, summaryFunction, ...)))
  
  # Binarize the data if desired
  if (binarize) {
    summarizedMat = summarizedMat > binMin
    summarizedMat = summarizedMat + 1 - 1
  }
  
  return(summarizedMat)
}



