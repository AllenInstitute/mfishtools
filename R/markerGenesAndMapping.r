#' mfishtools: Building Gene Sets and Mapping mFISH Data.
#'
#' This repository includes code for gene selection for spatial transcriptomics methods and for
#' mapping of spatial transcriptomics (or RNA-Seq data) onto a RNA-Seq reference. Specific topics include:
#' 1) Correlation-based mapping of cells to reference cell types
#' 2) Iterative building of gene panels a greedy algorithm with pre-defined constraints
#' 3) Visualizations related to gene mapping a gene panel selection
#'
#' @docType package
#' @name mfishtools
NULL


#' Table of confused clusters
#'
#' This function returns a table of the top confused clusters (assigned clusters
#'   incorrectly mapped)
#'
#' @param confusionProp confusion matrix (e.g., output from getConfusionMatrix).
#' @param count number of top confusions to show
#'
#' @return a 3 x count matrix of the top confused pairs of clusters with the three
#'   columns corresponding to mapped cluster, assigned cluster, and fraction of
#'   cells incorrectly mapped, respectively.
#'
#' @export
outputTopConfused <- function(confusionProp,
                              count = 10) {
  topConfused <- NULL
  mfp <- confusionProp
  diag(mfp) <- 0
  rn <- rownames(mfp)
  for (i in 1:count) {
    wt <- which(mfp == max(mfp), arr.ind = TRUE)
    wt <- t(wt[1, ])
    topConfused <- rbind(topConfused, c(rn[wt], max(mfp)))
    mfp[wt] <- 0
  }
  colnames(topConfused) <- c("foundCluster", "realCluster", "proportionOff")
  rownames(topConfused) <- NULL
  topConfused <- as.data.frame(topConfused)
  topConfused[, 3] <- as.numeric(topConfused[, 3])
  return(topConfused)
}

#' Confusion plot vs. confidence
#'
#' Produces line plots showing the percent of correctly mapped cells above a certain confidence value (or score).
#'   This is a wrapper for plot.
#'
#' @param foundClusterAndScore matrix where first column is found cluster and second column is
#'   confidence score (e.g., output from getTopMatch)
#' @param realCluster character vector of assigned clusters
#' @param ... additional parameters for the plot function
#'
#' @export
plotConfusionVsConfidence <- function(foundClusterAndScore,
                                      realCluster,
                                      RI = (31:100) / 100,
                                      main = "% mapping (blue) / correct (orange)",
                                      ylab = "Percent",
                                      xlab = "Fraction correctly mapped to leaf",
                                      type = "l",
                                      xlim = range(RI),
                                      ...) {
  fracMap <- fracRight <- NULL
  for (r in RI) {
    isMap <- foundClusterAndScore[, 2] >= r
    fracMap <- c(fracMap, round(1000 * mean(isMap)) / 10)
    isRight <- (realCluster == foundClusterAndScore[, 1])[isMap]
    fracRight <- c(fracRight, round(1000 * mean(isRight)) / 10)
  }
  plot(RI, fracMap,
    xlim = xlim, col = "blue", main = main,
    ylab = ylab, xlab = xlab, type = "l", ...
  )
  lines(RI, fracRight, col = "orange")
  abline(h = 5 * (0:20), col = "grey", lty = "dotted")
}


#' Confusion matrix
#'
#' Returns a confusion matrix of the found (mapped) vs. real (assigned) clusters.
#'
#' @param realCluster character vector of assigned clusters
#' @param foundCluster character vector of mapped clusters
#' @param proportions FALSE if the counts are to be returned and TRUE if the proportions are to be returned
#'
#' @export
getConfusionMatrix <- function(realCluster,
                               foundCluster,
                               proportions = TRUE) {
  realCluster <- as.character(realCluster)
  foundCluster <- as.character(foundCluster)
  lev <- sort(unique(c(realCluster, foundCluster)))
  realCluster <- factor(realCluster, levels = lev)
  foundCluster <- factor(foundCluster, levels = lev)
  confusion <- table(foundCluster, realCluster)
  if (proportions) {
    cs <- colSums(confusion)
    for (i in 1:dim(confusion)[1])
      confusion[i, ] <- confusion[i, ] / pmax(cs, 1e-08)
  }
  confusion
}




#' Filter genes for spatial transcriptomics panel
#'
#' Returns a set of genes for inclusion in a spatial transcriptomics panel based on a series of 
#'   hard-coded and user-defined constraints 
#'
#' @param summaryExpr Matrix of summarized expression levels for a given cluster.  Typically the median
#'   or mean should be used.  Rows are genes and columns are samples.  ROW NAMES MUST BE GENE SYMBOLS!
#' @param propExpr Proportion of cells expressed in each cluster for use with binary score calculation
#'   (default = summaryExpr, which is not recommended)
#' @param onClusters Vector indicating which clusters should be included in the gene panel (default
#'   is all clusters. Can be logical or numeric, or a character string of cluster names)
#' @param offClusters Vector indidicating from which clusters expression should be avoided
#' @param geneLength Optional vector of gene lengths in same order as summaryExpr.  Default is NULL
#' @param numBinaryGenes Number of genes to include in the final panel.  Genes are sorted by binary
#'   score using `getBetaScore` and this number of genes are chosen (default = 500)
#' @param minOn Minimum summary expression level in most highly expressed "on" cluster (default = 10)
#' @param maxOn Maximum summary expression level in most highly expressed "on" cluster (default = 250) 
#' @param maxOn Maximum summary expression level in most highly expressed "off" cluster (default = 50)
#' @param minLength Minimum gene length for marker gene selection.  Ignored if geneLength is not 
#'   provided (default = 960)
#' @param fractionOnClusters What is the maximum fraction of clusters in which a gene can be expressed
#'   (as defined by propExpr>0.5; default = 0.5).  This prevents nearly ubiquitous genes from selection
#' @param excludeGenes Which genes should be excluded from the analysis (default is none)
#' @param excludeFamilies Which gene classes or families should be excluded from the analysis?  More
#'   specifically, any gene that contain these strings of characters anywhere in the symbol will be
#'   excluded (default is "LOC","LINC","FAM","ORF","KIAA","FLJ","DKFZ","RIK","RPS","RPL","\\-"). 
#' 
#' @return A character vector of genes meeting all constraints
#'
#' @export
filterPanelGenes <- function(summaryExpr,
                             propExpr = summaryExpr,
                             onClusters = 1:dim(summaryExpr)[2],
                             offClusters = NULL,
                             geneLengths = NULL,
                             startingGenes = c("GAD1","SLC17A7"),
                             numBinaryGenes = 500,                             
                             minOn = 10,
                             maxOn = 250,
                             maxOff = 50,
                             minLength = 960,
                             fractionOnClusters = 0.5,
                             excludeGenes = NULL,
                             excludeFamilies = c("LOC","LINC","FAM","ORF","KIAA","FLJ","DKFZ","RIK","RPS","RPL","\\-")) {
  ## Variable check
  if(!is.element("matrix",class(summaryExpr))) {
    summaryExpr <- as.matrix(summaryExpr)
  }
  if(!is.numeric(summaryExpr[1,1])) {
    stop("summaryExpr must be a matrix of numeric values.")
  }
  if(is.null(rownames(summaryExpr))) {
    stop("Please provide summaryExpr with genes as row names.")
  }
  if(!is.numeric(fractionOnClusters)){
    stop("fractionOnClusters needs to be numeric.")
  }
  fractionOnClusters = fractionOnClusters[1]
  if(fractionOnClusters>1) fractionOnClusters = fractionOnClusters/100 # Assume if it is greater than 1 then a percentage was given
  genes  <- rownames(summaryExpr)
  genesU <- toupper(genes)
  
  ## Define excluded genes
  excludeFamilies <- toupper(excludeFamilies)
  excludeGenes <- is.element(genes,excludeGenes)
  if(length(excludeFamilies)>0){
    for(i in 1:length(excludeFamilies)){
      excludeGenes <- excludeGenes|grepl(excludeFamilies[i],genesU)
    }
  }
  
  ## Determine max expression levels in on and off clusters
  if(is.character(onClusters)){
    onClusters <- is.element(colnames(summaryExpr),onClusters)
  }
  if(is.numeric(onClusters)){
    onClusters <- is.element(1:dim(summaryExpr)[2],onClusters)
  }
  if(sum(onClusters)<2){
    stop("Please provide at least two onClusters.  If cluster names were provided, check colnames(summaryExpr).")
  }
  if(is.character(offClusters)){
    offClusters <- is.element(colnames(summaryExpr),offClusters)
  }
  if(is.numeric(offClusters)){
    offClusters <- is.element(1:dim(summaryExpr)[2],offClusters)
  }
  
  ## Now add the actual constraints from above
  maxExprOn <- apply(summaryExpr[,onClusters],1,max)
  if(sum(offClusters)>1){
    maxExprOff <- apply(summaryExpr[,offClusters],1,max)
  } else if(sum(offClusters)==1){
    maxExprOff <- summaryExpr[,offClusters]
  } else {
    maxExprOff <- maxExprOn*Inf # Essentially this is saying there is no off constraint
  }
  
  ## Set the gene lengths, if needed
  if(!is.null(geneLengths)) if(length(geneLengths)!=length(offClusters)){
    stop("geneLengths must be in the same order and of same length as rows of summaryExpr.")
  }
  if(!is.null(geneLengths)) if(!is.numeric(geneLengths)){
    stop("geneLengths must be numeric.")
  }
  if(is.null(geneLengths)){
    geneLengths <- maxExprOn*Inf # Essentially this is saying there is no off constraint
  }
  
  ## Determine the acceptable genes
  keepGenes <- (!excludeGenes)&(maxExprOn>minOn)&(maxExprOn<=maxOn)&(maxExprOff<=maxOff)&
               (geneLengths>=minLength)&(rowMeans(propExpr>0.5)<=fractionOnClusters)
  keepGenes[is.na(keepGenes)] <- FALSE                               
  
  ## Find the top binary genes (if needed) and return gene list
  message(paste(sum(keepGenes),"total genes pass constraints prior to binary score calculation."))
  if(sum(keepGenes)<=numBinaryGenes){
    warning("Fewer genes pass constraints than numBinaryGenes, so binary score was not calculated.")
    return(sort(union(runGenes,startingGenes)))
  }
  
  topBeta     <- getBetaScore(propExpr[keepGenes,onClusters],FALSE)
  runGenes    <- names(topBeta)[topBeta<=numBinaryGenes]
  runGenes    <- sort(union(runGenes,startingGenes))
  runGenes
}






#' Branch list
#'
#' Returns branches of a dendrogram in a specific format
#'
#' @param dend dendrogram for mapping.  Ignored if medianDat is passed
#' @param branches do not change from default
#' @param allTips do not change from default
#'
#' @return a list of branch information for use with leafToNodeMedians
#'
#' @export
getBranchList <- function(dend,
                          branches = list(),
                          allTips = as.character(dend %>% labels())) {
  library(dendextend)
  numBranch <- dend %>% nnodes()
  if (numBranch > 1) {
    lb <- attr(dend, "label")
    lab1 <- paste("BranchInTree___", lb, sep = "")
    branches[[lab1]] <- list()
    cn <- as.character(dend %>% labels())
    for (i in 1:length(dend)) {
      nm <- as.character(dend[[i]] %>% labels())
      branches[[lab1]][[attr(dend[[i]], "label")]] <- list(nm, setdiff(cn, nm))
      if ((length(nm) > 1) & (length(nm) < (length(allTips) - 1))) {
        lab2 <- paste("BranchVsAll___", attr(dend[[i]], "label"), sep = "")
        branches[[lab2]] <- list()
        branches[[lab2]][["branch"]] <- list(nm, setdiff(allTips, nm))
      }
      branches <- getBranchList(dend[[i]], branches, allTips)
    }
  } else {
    leaf <- dend %>% labels()
    lab <- paste("LeafOnly___", leaf, sep = "")
    branches[[lab]][["leaf"]] <- list(leaf, setdiff(allTips, leaf))
  }
  branches <- branches[order(substr(names(branches), 1, 8))]
  branches
}

#' Return mean node expression
#'
#' Define expression at a node as the MEAN expression for each leaf as default (using the
#'   median removes all specific marker genes!)
#'
#' @param dend dendrogram for mapping.  Ignored if medianDat is passed
#' @param medianDat median expression data at each node
#' @param branches a particular format of branch information from the dendrogram structure
#' @param fnIn function to use to wrap up to the node level (default = mean)
#'
#' @return a matrix of mean node expression (rows=genes, columns=nodes)
#'
#' @export
leafToNodeMedians <- function(dend,
                              medianDat,
                              branches = getBranchList(dend),
                              fnIn = mean) {
  library(dendextend)
  if (is.null(rownames(medianDat))) {
    rownames(medianDat) <- 1:dim(medianDat)[1]
  }
  allGenes <- rownames(medianDat)
  brNames <- names(branches)
  brNames <- brNames[grep("BranchInTree", brNames)]
  medianNode <- matrix(0, nrow = length(allGenes), ncol = length(brNames))
  rownames(medianNode) <- allGenes
  colnames(medianNode) <- brNames
  for (n in brNames) medianNode[, n] <-
      apply(medianDat[, c(branches[[n]][[1]][[1]], branches[[n]][[1]][[2]])], 1, fnIn)
  medianNode <- cbind(medianDat, medianNode)
  nameOrd <-
    dend %>%
    get_nodes_attr("label", id = 1:(dend %>% nnodes()))
  nameOrd[substr(nameOrd, 1, 1) == "n"] <-
    paste0("BranchInTree___", nameOrd[substr(nameOrd, 1, 1) == "n"])
  medianNode <- medianNode[, nameOrd]
  colnames(medianNode) <- gsub("BranchInTree___", "", colnames(medianNode))
  medianNode
}


#' Greedy algorithm for building marker gene panel
#'
#' This is the primary function that iteratively builds a marker gene panel, one gene at a
#'   time by iteratively adding the most informative gene to the existing gene panel.
#'
#' @param mapDat normalized data of the mapping (=reference) data set.
#' @param medianDat representative value for each leaf.  If not entered, it is calculated
#' @param clustersF cluster calls for each cell.
#' @param panelSize number of genes to include in the marker gene panel
#' @param subSamp number of random nuclei to select from each cluster (to increase speed);
#'   set as NA to not subsample
#' @param maxFcGene maximum number of genes to consider at each iteration (to increase speed)
#' @param qMin minimum quantile for fold change comparison (between 0 and 1, higher = more
#'   specific marker genes are included)
#' @param seed for reproducibility
#' @param currentPanel starting panel.  Default is NULL.
#' @param panelMin if there are fewer genes than this, the top number of these genes by fc
#'   rank are set as the starting panel.  Cannot be less than 2.
#' @param writeText should gene names and marker scores be output (default TRUE)
#' @param corMapping if TRUE (default) map by correlation; otherwise, map by Euclidean
#'   distance (not recommended)
#' @param optimize if 'FractionCorrect' (default) will seek to maximize the fraction of
#'   cells correctly mapping to final clusters
#'   if 'CorrelationDistance' will seek to minimize the total distance between actual
#'   cluster calls and mapped clusters
#'   if 'DendrogramHeight' will seek to minimize the total dendrogram height between
#'   actual cluster calls and mapped clusters
#' @param clusterDistance only used if optimize='CorrelationDistance'; a matrix (or
#'   vector) of cluster distances.  Will be calculated if NULL and if clusterGenes
#'   provided. (NOTE: order must be the same as medianDat and/or have column and row
#'   names corresponding to clusters in clustersF)
#' @param clusterGenes a vector of genes used to calculate the cluster distance.
#'   Only used if optimize='CorrelationDistance' and clusterDistance=NULL.
#' @param dend only used if optimize='DendrogramHeight' dendrogram; will error out of not provided
#' @param percentSubset for each iteration the function can subset the set of possible
#'   genes to speed up the calculation.
#'
#' @return an ordered character vector corresponding to the marker gene panel
#'
#' @export
buildMappingBasedMarkerPanel <- function(mapDat,
                                         medianDat = NA,
                                         clustersF = NA,
                                         panelSize = 50,
                                         subSamp = 20,
                                         maxFcGene = 1000,
                                         qMin = 0.75,
                                         seed = 10,
                                         currentPanel = NULL,
                                         panelMin = 5,
                                         writeText = TRUE,
                                         corMapping = TRUE,
                                         optimize = "FractionCorrect",
                                         clusterDistance = NULL,
                                         clusterGenes = NULL,
                                         dend = NULL,
                                         percentSubset = 100) {

  # Return an error if optimize='DendrogramHeight'
  # and a dendrogram is not provided
  if ((optimize == "DendrogramHeight") & is.null(dend)) {
    return("Error: dendrogram not provided")
  }

  # CALCULATE THE MEDIAN
  if (is.na(medianDat[1])) {
    names(clustersF) <- colnames(mapDat)
    medianDat <- do.call("cbind", tapply(
      names(clustersF), clustersF, function(x) rowMedians(mapDat[, x])
    ))
    rownames(medianDat) <- rownames(mapDat)
  }
  if (is.null(rownames(medianDat))) {
    rownames(medianDat) <- rownames(mapDat)
  }

  # Convert the dendrogram height into a correlation
  # distance if dendrogram height is entered as the
  # option
  if (optimize == "FractionCorrect") {
    clusterDistance <- NULL
  }
  if (optimize == "CorrelationDistance") {
    if (is.null(clusterDistance)) {
      corDist <- function(x) return(as.dist(1 - cor(x)))
      clusterGenes <- intersect(clusterGenes, rownames(medianDat))
      clusterDistance <- as.matrix(corDist(medianDat[clusterGenes, ]))
    }
    if (is.matrix(clusterDistance)) {
      if (!is.null(rownames(clusterDistance))) {
        clusterDistance <- clusterDistance[colnames(medianDat), colnames(medianDat)]
      }
      clusterDistance <- as.vector(clusterDistance)
    }
  }

  if (optimize == "DendrogramHeight") {
    lcaTable <- makeLCAtable(dend)
    clusterDistance <- 1 - getNodeHeight(dend)[lcaTable]
    optimize <- "clusterDistance"
  }

  # TAKE THE TOP DEX GENES Use fold change (rather
  # than beta) because this function only receives
  # median as input
  fcDiff <- rank(apply(medianDat, 1, function(x) return(diff(quantile(x, c(1, qMin))))))
  if (dim(medianDat)[1] > maxFcGene) {
    kpGene <- names(fcDiff)[fcDiff <= maxFcGene]
    mapDat <- mapDat[kpGene, ]
    medianDat <- medianDat[kpGene, ]
  }

  panelMin <- max(2, panelMin)
  if (length(currentPanel) < panelMin) {
    panelMin <- max(2, panelMin - length(currentPanel))
    currentPanel <- unique(c(currentPanel, names(sort(fcDiff))[1:panelMin]))
    if (writeText) {
      print(paste("Setting starting panel as:", paste(currentPanel, sep = ", ", collapse = ", ")))
    }
  }

  # FIND THE NEXT GENE IN THE PANEL, IF THE DESIRED
  # PANEL SIZE IS NOT REACHED
  if (length(currentPanel) < panelSize) {

    # SUBSAMPLE
    if (!is.na(subSamp)) {
      kpSamp <- subsampleCells(clustersF, subSamp, seed)
      mapDat <- mapDat[, kpSamp]
      clustersF <- clustersF[kpSamp]
      subSamp <- NA
    }

    # CORRELATION MAPPING FOR EACH POSSIBLE ADDITION OF
    # ONE GENE
    otherGenes <- setdiff(rownames(mapDat), currentPanel)
    if (percentSubset < 100) {
      # Only look at a subset of genes if desired
      set.seed(seed + length(currentPanel))
      otherGenes <- otherGenes[sort(sample(
        1:length(otherGenes),
        ceiling(length(otherGenes) * percentSubset / 100)
      ))]
    }
    matchCount <- rep(0, length(otherGenes))
    clustIndex <- match(clustersF, colnames(medianDat))
    for (i in 1:length(otherGenes)) {
      ggnn <- c(currentPanel, otherGenes[i])
      if (corMapping) {
        corMapTmp <- corTreeMapping(mapDat = mapDat, medianDat = medianDat, genesToMap = ggnn)
      }
      if (!corMapping) {
        corMapTmp <- distTreeMapping(mapDat = mapDat, medianDat = medianDat, genesToMap = ggnn)
      }
      corMapTmp[is.na(corMapTmp)] <- -1
      topLeafTmp <- getTopMatch(corMapTmp)
      if (is.null(clusterDistance)) {
        matchCount[i] <- mean(clustersF == topLeafTmp[, 1])
      } else {
        tmpVal <- dim(medianDat)[2] * (match(topLeafTmp[, 1], colnames(medianDat)) - 1) + clustIndex # NEED TO CHECK THIS!!!!!!!
        matchCount[i] <- -mean(clusterDistance[tmpVal])
      }
    }
    wm <- which.max(matchCount)
    addGene <- as.character(otherGenes)[wm]
    if (writeText) {
      if (optimize == "FractionCorrect") {
        print(paste(
          "Added", addGene, "with", signif(matchCount[wm], 3),
          "now matching [", length(currentPanel), "]."
        ))
      } else {
        print(paste(
          "Added", addGene, "with average cluster distance",
          -signif(matchCount[wm], 3), "[", length(currentPanel), "]."
        ))
      }
    }
    currentPanel <- c(currentPanel, addGene)
    currentPanel <- buildMappingBasedMarkerPanel(
      mapDat = mapDat, medianDat = medianDat, clustersF = clustersF,
      panelSize = panelSize, subSamp = subSamp, maxFcGene = maxFcGene,
      qMin = qMin, seed = seed, currentPanel = currentPanel, panelMin = panelMin,
      writeText = writeText, corMapping = corMapping, optimize = optimize,
      clusterDistance = clusterDistance, clusterGenes = clusterGenes, dend = dend,
      percentSubset = percentSubset
    )
  }
  currentPanel
}


#' Correlation between nodes and leafs (deprecated)
#'
#' Returns the correlation between expression of each cell and representative
#'   value for each node and leaf.  NOTE: this function is unstable and will
#'   eventually be merged with corTreeMapping.
#'
#' @param dend dendrogram for mapping.  Ignored if medianDat is passed
#' @param refDat normalized data of the REFERENCE data set.  Ignored if medianExpr
#'   and propExpr are passed
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the
#'   data onto itself.
#' @param medianExpr representative value for each leaf.  If not entered, it is
#'   calculated
#' @param propExpr proportion of cells in each type expressing a given gene.  If not
#'   entered, it is calculated
#' @param filterMatrix a matrix of TRUE/FALSE values to indicate whether a given
#'   cluster is possible
#' @param clusters cluster calls for each cell.  Ignored if medianExpr and propExpr
#'   are passed
#' @param numberOfGenes how many variables genes
#' @param outerLimitGenes choose different numberOfGenes per cell from the top overall
#'   outerLimitGenes (to speed up function)
#' @param genesToMap which genes to include in the correlation mapping
#' @param use,... additional parameters for cor
#'
#' @return a matrix of correlation values with rows as mapped cells and columns as clusters
#'
#' @export
corTreeMapping_withFilter <- function(dend = NA,
                                      refDat = NA,
                                      mapDat = refDat,
                                      medianExpr = NA,
                                      propExpr = NA,
                                      filterMatrix = NA,
                                      clusters = NA,
                                      numberOfGenes = 1200,
                                      outerLimitGenes = 7200,
                                      rankGeneFunction = function(x) getBetaScore(x, returnScore = FALSE),
                                      use = "p",
                                      ...) {

  # -- prepare the data
  if (is.na(medianExpr[1])) {
    names(clusters) <- colnames(refDat)
    medianExpr <- do.call("cbind", tapply(
      names(clusters), clusters, function(x) rowMedians(refDat[, x])
    ))
    rownames(medianExpr) <- rownames(refDat)
  }
  if (is.na(propExpr[1])) {
    names(clusters) <- colnames(refDat)
    medianExpr <- do.call("cbind", tapply(
      names(clusters), clusters, function(x) rowMeans(refDat[, x] > 1)
    ))
    rownames(propExpr) <- rownames(refDat)
  }
  filterMatrix <- filterMatrix[, colnames(medianExpr)]

  # -- take the top outerLimitGenes for the
  # proportion and median
  rankGn <- rankGeneFunction(propExpr)
  kpGn <- rankGn <= outerLimitGenes
  medianDat <- medianExpr[kpGn, ]
  propDat <- propExpr[kpGn, ]

  # -- find all possible filters
  filterVec <- apply(filterMatrix, 1, function(x) paste(x, collapse = "|", sep = "|"))
  vecs <- unique(filterVec)

  # -- find all gene lists based on these filters
  geneLists <- list()
  for (v in vecs) {
    kp <- filterMatrix[which(filterVec == v)[1], ]
    if (sum(kp) > 1) {
      geneLists[[v]] <- rownames(propDat)[rankGeneFunction(propDat[, kp]) <= numberOfGenes]
    } else {
      geneLists[[v]] <- rownames(propDat)[1:2]
    }
  }

  ## -- find the correlations
  kpVar <- intersect(
    names(rankGn)[rankGn <= numberOfGenes],
    intersect(rownames(mapDat), rownames(medianDat))
  )
  corrVar <- cor(mapDat[kpVar, ], medianDat[kpVar, ], use = use, ...)
  for (v in vecs) {
    kpRow <- filterVec == v
    kpCol <- filterMatrix[which(filterVec == v)[1], ]
    if (sum(kpCol) > 1) {
      kpVar <- intersect(geneLists[[v]], intersect(rownames(mapDat), rownames(medianDat)))
      corrVar[kpRow, kpCol] <- cor(mapDat[kpVar, kpRow], medianDat[kpVar, kpCol], use = use) # ,...)
    }
  }
  corrVar <- corrVar * filterMatrix[, colnames(corrVar)]
  corrVar
}


#' (Euclidean) distance mapping
#'
#' Returns the distance between expression of each cell and representative value for each node and
#'   leaf (default is based on euclidean distance).  In our hands this is does not work very well.
#'
#' @param dend dendrogram for mapping.  Ignored if medianDat is passed
#' @param refDat normalized data of the REFERENCE data set.  Ignored if medianDat is passed
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param medianDat representative value for each leaf and node.  If not entered, it is calculated
#' @param clusters cluster calls for each cell.  Ignored if medianDat is passed
#' @param genesToMap which genes to include in the correlation mapping
#' @param returnSimilarity FALSE to return distance, TRUE to return something like a similarity
#' @param use,... additional parameters for dist (for back-compatiblity; doesn't work)
#'
#' @return matrix of Euclidean distances between cells (rows) and clusters (columns)
#'
#' @export
distTreeMapping <- function(dend = NA,
                            refDat = NA,
                            mapDat = refDat,
                            medianDat = NA,
                            clusters = NA,
                            genesToMap = rownames(mapDat),
                            returnSimilarity = TRUE,
                            use = "p",
                            ...) {
  library(pdist)

  if (is.na(medianDat[1])) {
    names(clusters) <- colnames(refDat)
    medianDat <- do.call("cbind", tapply(
      names(clusters), clusters, function(x) rowMedians(refDat[, x])
    ))
    rownames(medianDat) <- rownames(refDat)
    medianDat <- leafToNodeMedians(dend, medianDat)
  }
  kpVar <- intersect(genesToMap, intersect(rownames(mapDat), rownames(medianDat)))
  if (length(kpVar) == 1) {
    kpVar <- c(kpVar, kpVar)
  }
  eucDist <- as.matrix(pdist(t(mapDat[kpVar, ]), t(medianDat[kpVar, ]), ...))
  rownames(eucDist) <- colnames(mapDat)
  colnames(eucDist) <- colnames(medianDat)
  if (!returnSimilarity) {
    return(eucDist)
  }
  eucDist <- sqrt(eucDist / max(eucDist))
  1 - eucDist
}


#' Tree-based mapping
#'
#' Returns the mapping membership of each cell to each node and leaf using a
#'   tree-based method.  This is a wrapper function for map_dend.
#'
#' @param dend dendrogram for mapping
#' @param refDat normalized data of the REFERENCE data set
#' @param clustersF factor indicating which cluster each cell type is actually assigned to
#'   in the reference data set
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#'   to map to. otherwise, a random branch is chosen.
#' @param seed added for reproducibility
#'
#' @return a matrix of confidence scores (from 0 to 100) with rows as cells and columns
#'   as tree node/leafs.  Values indicate the fraction of permutations in which the cell
#'   mapped to that node/leaf using the subset of cells/genes in map_dend
#'
#' @export
rfTreeMapping <- function(dend,
                          refDat,
                          clustersF,
                          mapDat = refDat,
                          p = 0.7,
                          low.th = 0.15,
                          seed = 1) {
  refDat <- as.matrix(refDat)
  mapDat <- as.matrix(mapDat)
  pseq.cells <- colnames(mapDat)
  # isMarker   <- ifelse(is.na(adjustMarkers[1]), FALSE, TRUE) # Doesn't do anything

  pseq.mem <- sapply(1:100, function(i) {
    j <- i
    if (i %% 25 == 0) print(i)
    go <- TRUE
    while (go) {
      # Allow for failures
      j <- j + 1000
      set.seed(j + seed) # Added for reproducibility
      tmp <- try(map_dend(dend, clustersF, refDat,
        mapDat, pseq.cells,
        p = p, low.th = low.th
      ))
      if (length(tmp) > 1) go <- FALSE
    }
    tmp
  }, simplify = F)

  memb <- unlist(pseq.mem)
  memb <- data.frame(cell = names(memb), cl = memb)
  memb$cl <- factor(memb$cl, levels = get_nodes_attr(dend, "label"))
  memb <- table(memb$cell, memb$cl)
  memb <- memb / 100
  return(memb)
}


#' Tree-based mapping
#'
#' Returns the mapping membership of each cell to each node and leaf using a
#'   tree-based method.  This is a wrapper function for map_dend.
#'
#' @param dend dendrogram for mapping
#' @param cl factor indicating which cluster each cell type is actually assigned to
#'   in the reference data set
#' @param dat normalized data of the REFERENCE data set
#' @param map.dat normalized data of the MAPPING data set.  Default is to map the
#'   data onto itself.
#' @param p proportion of marker genes to include in each iteration of the mapping
#'   algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide
#'   on which branch to map to. otherwise, a random branch is chosen.
#' @param default.markers not used
#'
#' @return a matrix of confidence scores (from 0 to 100) with rows as cells and columns
#'   as tree node/leafs.  Values indicate the fraction of permutations in which the cell
#'   mapped to that node/leaf using the subset of cells/genes in map_dend
#'
#' @export
map_dend <- function(dend,
                     cl,
                     dat,
                     map.dat,
                     select.cells,
                     p = 0.8,
                     low.th = 0.2,
                     default.markers = NULL) {
  final.cl <- c(setNames(rep(attr(dend, "label"), length(select.cells)), select.cells))
  if (length(dend) <= 1) {
    return(final.cl)
  }
  markers <- attr(dend, "markers")
  markers <- markers[names(markers) %in% row.names(map.dat)]
  cl.g <- sapply(dend, labels, simplify = F)
  names(cl.g) <- 1:length(cl.g)
  select.cl <- cl[cl %in% unlist(cl.g)]
  ### Sampling the cells from the reference cluster
  cells <- unlist(tapply(names(select.cl), select.cl, function(x) sample(x, round(length(x) * p))))
  genes <- names(markers)
  genes <- union(genes, default.markers)
  ### Compute reference cluster median based on
  ### subsampled cells
  cl.med <- do.call("cbind", tapply(
    cells, droplevels(cl[cells]),
    function(x) rowMedians(dat[genes, x, drop = F])
  ))
  row.names(cl.med) <- genes
  ### determine which branch to take.
  mapped.cl <- resolve_cl(cl.g, cl.med, markers, dat,
    map.dat, select.cells,
    p = p, low.th = low.th
  )
  if (length(mapped.cl) > 0) {
    for (i in unique(mapped.cl)) {
      select.cells <- names(mapped.cl)[mapped.cl == i]
      if (length(select.cells) > 0) {
        final.cl <- c(final.cl, map_dend(dend[[as.integer(i)]],
          cl, dat, map.dat, select.cells,
          p = p, low.th = low.th
        ))
      }
    }
  }
  return(cl = final.cl)
}


#' Tree-based mapping (internal)
#'
#' Returns the mapped cluster call of each cell to each leaf. This function is called by map_dend
#'
#' @param cl.g all clusters
#' @param cl.med cluster medians
#' @param markers gene markers
#' @param dat normalized data of the REFERENCE data set
#' @param map.dat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param select.cells which cells to use?
#' @param p proportion of marker genes to include in each iteration of the mapping algorithm.
#' @param low.th the minimum difference in Pearson correlation required to decide on which branch
#'   to map to. otherwise, a random branch is chosen.
#'
#' @return a vector of the mapped cluster
#'
#' @export
resolve_cl <- function(cl.g,
                       cl.med,
                       markers,
                       dat,
                       map.dat,
                       select.cells,
                       p = 0.7,
                       low.th = 0.2) {
  library(matrixStats)
  ##
  genes <- names(markers)[markers > 0]
  tmp.cl <- unlist(cl.g)

  ### For each branch point, find the highest
  ### expression cluster.
  tmp.med <- sapply(cl.g, function(g) rowMaxs(cl.med[genes, g, drop = F]))
  row.names(tmp.med) <- genes
  ### Make sure the genes are discriminative between
  ### all the branches.
  genes <- genes[rowMaxs(tmp.med) - rowMins(tmp.med) > 1]

  ### Sample the markers based on the weigts. TO DO:
  ### randomforest sometimes give importance value of
  ### 0. adjust for that.
  genes <- sample(genes, round(length(genes) * p), prob = markers[genes])

  ### Compute the correlation with the median cluster
  ### profile. add drop=F
  cl.cor <- cor(map.dat[genes, select.cells, drop = F], cl.med[genes, tmp.cl, drop = F])
  cl.cor[is.na(cl.cor)] <- 0
  ### Compute the best match in each branch.
  tmp.score <- do.call("cbind", sapply(cl.g, function(x) rowMaxs(cl.cor[,
      x,
      drop = F
    ]), simplify = F))
  row.names(tmp.score) <- row.names(cl.cor)
  #### Determine the best match.
  best.score <- setNames(rowMaxs(tmp.score), row.names(tmp.score))
  ### determine the difference from the best match.
  diff.score <- best.score - tmp.score

  #### Give up on cells can't be discriminated,choose
  #### one branch randomly.
  unresolved.cl <- row.names(tmp.score)[rowSums(diff.score < low.th) == ncol(diff.score)]
  mapped.cl <- setNames(sample(colnames(tmp.score), length(unresolved.cl), replace = T), unresolved.cl)

  ### Cells mapped to one or more branches.
  mapped.cells <- setdiff(row.names(cl.cor), unresolved.cl)
  ### For binary branch, done already
  if (length(cl.g) == 2) {
    mapped.cl <- c(mapped.cl, setNames(colnames(diff.score)[apply(diff.score[mapped.cells,
      ,
      drop = F
    ], 1, which.min)], mapped.cells))
    return(mapped.cl)
  }
  ## The remaining options for mapped cells
  tmp.cl <- sapply(mapped.cells, function(x) colnames(diff.score)[which(diff.score[x, ] < low.th)], simplify = F)
  ### cells with multiple options
  resolve.cells <- names(tmp.cl)[sapply(tmp.cl, length) > 1]
  ### cells with only one option. Not further job.
  mapped.cells <- setdiff(mapped.cells, resolve.cells)
  if (length(mapped.cells) > 0) {
    mapped.cl <- c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
  }
  ### Resolve further options.
  if (length(resolve.cells) > 0) {
    tmp.cat <- sapply(tmp.cl[resolve.cells], function(x) paste(x, collapse = " "))
    for (cat in unique(tmp.cat)) {
      tmp.cl <- unlist(strsplit(cat, " "))
      select.cells <- names(tmp.cat)[tmp.cat == cat]
      mapped.cl <- c(mapped.cl, resolve_cl(cl.g[tmp.cl],
        cl.med, markers, dat, map.dat, select.cells,
        p = p, low.th = low.th
      ))
    }
  }
  return(mapped.cl)
}



#' Get top leaf match
#'
#' Returns the top leaf match for each cell and the corresponding fraction mapping there.
#'
#' @param memb.cl membership scores for each leaf
#'
#' @return a matrix where first column is found cluster and second column is confidence score
#'
#' @export
getTopMatch <- function(memb.cl) {
  tmp.cl <- apply(memb.cl, 1, function(x) {
    y <- which.max(x)[1]
    as.character(c(colnames(memb.cl)[y], x[y]))
  })
  rfv <- as.data.frame(t(tmp.cl))
  rfv[, 2] <- as.numeric(rfv[, 2])
  colnames(rfv) <- c("TopLeaf", "Value")
  rfv[is.na(rfv[, 1]), 1] <- "none"
  rfv[is.na(rfv[, 2]), 2] <- 0
  rfv
}


#' Generate reference set of pseudo-cells
#'
#' Creates a new reference set as input for cellToClusterMapping_byRank, where each 'cell' is the
#'   combiniation of several cells and this is run several times using different subsets of cells.
#' @param refDat  normalized data of the REFERENCE data set
#' @param clustersF factor indicating which cluster each cell type is actually assigned to in the reference data set
#' @param genesToMap which genes to include in the correlation mapping
#' @param cellsPerMerge Number of cells to include in each combo cell
#' @param numberOfMerges Number of combo cells to include per cell type
#' @param mergeFunction function for combining cells into combo cells (use rowMeans or rowMedians)
#' @param seed for resproducibility
#'
#' @return list where first element is data matrix of multi-cells by genes and
#'   second element is a vector of corresponding clusters
#'
#' @export
generateMultipleCellReferenceSet <- function(refDat,
                                             clustersF,
                                             genesToUse = rownames(refDat),
                                             cellsPerMerge = 5,
                                             numberOfMerges = 10,
                                             mergeFunction = rowMedians,
                                             seed = 1) {
  if (!is.factor(clustersF)) {
    clustersF <- factor(clustersF)
  }
  names(clustersF) <- colnames(refDat)
  clusts <- levels(clustersF)
  refUse <- refDat[intersect(rownames(refDat), genesToUse), ]
  refOut <- matrix(nrow = dim(refUse)[1], ncol = numberOfMerges * length(clusts))
  rownames(refOut) <- rownames(refUse)
  index <- 0
  for (k in 1:numberOfMerges) {
    i <- NULL
    for (cl in 1:length(clusts)) {
      val <- which(clustersF == clusts[cl])
      set.seed(seed + index + cl + k)
      i <- c(i, sample(val, min(cellsPerMerge, length(val))))
    }
    i <- is.element(1:length(clustersF), i)
    refOut[, (index + 1):(index + length(clusts))] <- do.call(
      "cbind",
      tapply(names(clustersF[i]), clustersF[i], function(x) rowMedians(refUse[, i][, x]))
    )
    index <- index + length(clusts)
  }
  list(data = refOut, clusters = rep(clusts, numberOfMerges))
}


#' Cell-based cluster mapping
#'
#' Maps cells to clusters by correlating every mapped cell with every reference cell,
#'   ranking the cells by correlation, and the reporting the cluster with the lowest average rank.
#'
#' @param mapDat normalized data of the MAPPING data set.
#' @param refDat normalized data of the REFERENCE data set
#' @param clustersF factor indicating which cluster each cell type is actually assigned
#'   to in the reference data set
#' @param genesToMap character vector of which genes to include in the correlation mapping
#' @param mergeFunction function for combining ranks; the tested choices are rowMeans or
#'   rowMedians (default)
#' @param useRank use the rank of the correlation (default) or the correlation itself to
#'   determine the top cluster
#' @param use additional parameter for cor (use='p' as default)
#' @param method additional parameter for cor (method='p' as default)
#'
#' @return a two column data matrix where the first column is the mapped cluster and the second
#'   column is a confidence call indicating how close to the top of the ranked list cells of the
#'   assigned cluster were located relative to their best possible location in the ranked list.
#'   This confidence score seems to be a bit more reliable than correlation at determining how
#'   likely a cell in a training set is to being correctly assigned to the training cluster.
#'
#' @export
cellToClusterMapping_byRank <- function(mapDat,
                                        refDat,
                                        clustersF,
                                        genesToMap = rownames(mapDat),
                                        mergeFunction = rowMedians,
                                        useRank = TRUE,
                                        use = "p",
                                        method = "p") {
  if (is.null(names(clustersF))) names(clustersF) <- colnames(refDat)
  kpVar <- intersect(genesToMap, intersect(rownames(mapDat), rownames(refDat)))
  corrVar <- cor(mapDat[kpVar, ], refDat[kpVar, ], use = use, method = method)
  corrVar[corrVar > 0.999999] <- NA # assume any perfect correlation is either an self-to-self mapping, or a mapping using exactly 1 non-zero gene
  if (useRank) rankVar <- t(apply(-corrVar, 1, rank, na.last = "keep"))
  if (!useRank) rankVar <- -corrVar
  colnames(rankVar) <- names(clustersF) <- paste0("n", 1:length(clustersF))
  clMean <- do.call("cbind", tapply(
    names(clustersF),
    clustersF, function(x) match.fun(mergeFunction)(rankVar[, x], na.rm = TRUE)
  ))
  clMin <- apply(clMean, 1, min, na.rm = TRUE)
  clMin[is.na(clMin)] <- 0
  clMin[clMin == Inf] <- 1e+09
  clBest <- colnames(clMean)[apply(clMean, 1, function(x) return(which.min(x)[1]))]
  clBest[is.na(clBest)] <- colnames(clMean)[1]
  if (useRank) clScore <- (table(clustersF)[clBest] / 2) / pmax(clMin, 1e-11)
  if (!useRank) clScore <- -clMin
  rfv <- data.frame(TopLeaf = clBest, Score = as.numeric(as.character(clScore)))
  rownames(rfv) <- rownames(corrVar)
  rfv
}


#' Correlation-based cluster mapping
#'
#' Primary function for doing correlation-based mapping to cluster medians.  This is wrapper for cor
#'   and returns a correlation matrix.
#'
#' @param mapDat normalized data of the MAPPING data set.  Default is to map the data onto itself.
#' @param medianDat representative value for each leaf and node.  If not entered, it is calculated
#' @param dend dendrogram for mapping.  If provided, correlations to nodes are also returned
#' @param refDat normalized data of the REFERENCE data set.  Ignored if medianDat is passed
#' @param clusters  cluster calls for each cell.  Ignored if medianDat is passed
#' @param genesToMap which genes to include in the correlation mapping
#' @param use additional parameter for cor (use='p' as default)
#' @param method additional parameter for cor (method='p' as default)
#'
#' @return matrix with the correlation between expression of each cell and representative value for
#'   each leaf and node
#'
#' @export
corTreeMapping <- function(mapDat,
                           medianDat,
                           dend = NULL,
                           refDat = NA,
                           clusters = NA,
                           genesToMap = rownames(mapDat),
                           use = "p",
                           method = "p") {
  if (is.na(medianDat[1])) {
    names(clusters) <- colnames(refDat)
    medianDat <- do.call("cbind", tapply(
      names(clusters), clusters, function(x) rowMedians(refDat[, x])
    ))
    rownames(medianDat) <- rownames(refDat)
  }
  if (!is.null(dend)) {
    medianDat <- leafToNodeMedians(dend, medianDat)
  }
  kpVar <- intersect(genesToMap, intersect(rownames(mapDat), rownames(medianDat)))
  corrVar <- cor(mapDat[kpVar, ], medianDat[kpVar, ], use = use, method = method)
  return(corrVar)
}


#' Gets subtree labels for lca function.
#'
#' @param dend a cluster dendrogram
#'
#' @return vector of subtree labels
#'
#' @export
get_subtree_label <- function(dend) {
  library(dendextend)
  l <- attr(dend, "label")
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      l <- c(l, get_subtree_label(dend[[i]]))
    }
  }
  return(l)
}


#' Get lowest common ancestor (defined cluster pairs)
#'
#' Maps a cluster back up the tree to the first node where the mapped and correct clusters agree.
#'
#' @param dend a cluster dendrogram
#' @param l1 a vector of node labels
#' @param l2 a second fector of node labels (of the same length as l1)
#' @param l do not adjust; required for recursive function
#'
#' @return The function will return a vector for lowest common ancestor for every pair of nodes in l1 and l2
#'
#' @export
lca <- function(dend,
                l1,
                l2,
                l = rep(attr(dend, "label"), length(l1))) {
  library(dendextend)
  node.height <- setNames(get_nodes_attr(dend, "height"), get_nodes_attr(dend, "label"))
  if (length(dend) > 1) {
    for (i in 1:length(dend)) {
      tmp.l <- attr(dend[[i]], "label")
      labels <- get_subtree_label(dend[[i]])
      select <- l1 %in% labels & l2 %in% labels
      if (sum(select) > 0) {
        select <- which(select)[node.height[l[select]] > node.height[tmp.l]]
        l[select] <- tmp.l
        l <- lca(dend[[i]], l1, l2, l)
      }
    }
  }
  l
}


#' Get lowest common ancestor (all cluster pairs in tree)
#'
#' Calculates the vector for lowest common ancestor for every pair of leaves in a tree and returns a
#'   vector in a specific format for faster look-up.
#'
#' @param dend a cluster dendrogram
#' @param includeInternalNodes should internal nodes be included in the output?
#' @param verbose if TRUE, status will be printed to the screen, since function is relatively slow
#'   for large trees (default FALSE)
#'
#' @return The function will return a vector for lowest common ancestor for every pair of leaves
#'   in dend.  Vector names are l1|||l2 for string parsing in other functions.
#'
#' @export
makeLCAtable <- function(dend,
                         includeInternalNodes = FALSE,
                         verbose = FALSE) {
  library(dendextend)
  nodes <- get_leaves_attr(dend, "label")
  if (includeInternalNodes) nodes <- get_nodes_attr(dend, "label")
  out <- nm <- rep(0, length(nodes)^2)
  i <- 1
  for (l1 in nodes) {
    if (verbose) print(l1)
    for (l2 in nodes) {
      nm[i] <- paste(l1, l2, sep = "|||")
      out[i] <- lca(dend, l1, l2)
      i <- i + 1
    }
  }
  names(out) <- nm
  out
}


#' Plot dendrogram
#'
#' Plots a dendrogram with set not colors, shapes, sizes and labels.  This is a wrapper for plot.
#'
#' @param tree a dendrogram object
#' @param value numeric vector corresponding to the size of each node
#' @param cexScale a global cex multiplier for node sizes
#' @param margins set the margins using par(mar=margins)
#' @param cols vector of node colors (or a single value)
#' @param pch vector of node pch shapes (or a single value)
#' @param ... additional parameters for the plot function
#'
#' @export
plotNodes <- function(tree,
                      value = rep(1, length(labels(tree))),
                      cexScale = 2,
                      margins = c(10, 5, 2, 2),
                      cols = "black",
                      pch = 19,
                      ...) {
  tree <- set(tree, "nodes_pch", pch)
  tree <- set(tree, "nodes_col", cols)
  tree <- set(tree, "labels_cex", 1)
  treeN <- set(tree, "nodes_cex", cexScale * value)
  par(mar = margins)
  plot(treeN, ylab = "height", ...)
}



#' Build panel for one cluster (beta)
#'
#' This UNTESTED function finds the best small marker panel for marking a single cluster, using
#'   proportion difference as the metric for determining the starting panel.
#'
#' @param mapDat normalized data of the mapping (=reference) data set.
#' @param clustersF cluster calls for each cell.
#' @param medianDat median value for each leaf
#' @param propIn proportions of cells with expression > 1 in each leaf
#' @param clust which cluster to target?
#' @param subSamp number of random nuclei to select from each cluster, EXCEPT the target cluster;
#'   set as NA to not subsample
#' @param seed for reproducibility
#' @param maxSize maximum size of marker gene panel
#' @param dexCutoff criteria for stopping: when improvement in fraction of cells properly mapped
#'   dips below this value
#' @param topGeneCount number of top genes by proportion to consider
#'
#' @return a matrix of the top marker genes for each cluster.  Output matrix includes five columns:
#'   clust = cluster; panel = ordered genes in the panel for that cluster; onCorrect = fraction of
#'   correctly assigned cells in cluster; offCorrect = fraction of cells correctly assigned outside
#'   of cluster; dexTotal = additional dex explained by last gene added.
#'
#' @export
buildPanel_oneCluster <- function(mapDat,
                                  clustersF,
                                  medianDat = NA,
                                  propIn = NA,
                                  clust = as.character(clustersF[1]),
                                  subSamp = NA,
                                  seed = 10,
                                  maxSize = 20,
                                  dexCutoff = 0.001,
                                  topGeneCount = 100) {

  # SUBSAMPLE
  if (!is.na(subSamp)) {
    kpSamp <- subsampleCells(clustersF, subSamp, seed)
    kpSamp[as.character(clustersF) == clust] <- TRUE
    mapDat <- mapDat[, kpSamp]
    clustersF <- clustersF[kpSamp]
    subSamp <- NA
  }

  # REFORMAT THE CLUSTER VARIABLE
  clust <- as.character(clust)
  clustersIn <- clustersF
  clustersF <- rep(clust, length(clustersIn))
  clustersF[as.character(clustersIn) != clust] <- "other"
  clustersF <- factor(clustersF, levels = c(clust, "other"))

  # CALCULATE THE PROPORTION OF CELLS EXPRESSED IN
  # EACH CLUSTERS, AND THE MEDIANS (SEND TO OTHER
  # FUNCTIONS AS MEDIAN)
  names(clustersF) <- colnames(mapDat)
  if (is.na(propIn[1])) {
    propIn <- do.call("cbind", tapply(
      names(clustersF), clustersF, function(x) rowMeans(mapDat[, x] >= 1)
    ))
  }
  rownames(propIn) <- rownames(mapDat)
  if (is.na(medianDat[1])) {
    medianDat <- do.call("cbind", tapply(
      names(clustersF), clustersF, function(x) rowMedians(mapDat[, x])
    ))
  } # switched clustersIn to clustersF
  rownames(medianDat) <- rownames(mapDat)

  # FIND THE BEST GENE IN THE PANEL, UNTIL THE
  # DESIRED PANEL SIZE IS REACHED
  propDat <- cbind(propIn[, clust], rowMeans(propIn[, colnames(propIn) != clust]))
  colnames(propDat) <- c(clust, "Other")
  panel <- onCorrect <- offCorrect <- dexTotal <- NULL
  first <- TRUE
  tt <- dex <- 0
  propDex <- propDat[, 1] - propDat[, 2]
  topMark <- names(-sort(-propDex))[1:topGeneCount] # Some semi-reasonable way to cut down the gene count
  while ((((length(panel) < maxSize) & (dex > dexCutoff))) | first) {

    # EUCLIDEAN MAPPING FOR EACH POSSIBLE ADDITION OF ONE GENE
    first <- FALSE
    otherGenes <- setdiff(topMark, panel)
    matchCount <- offTarget <- totalCount <- rep(0, length(otherGenes))
    for (i in 1:length(otherGenes)) {
      corMapTmp <- distTreeMapping(
        mapDat = mapDat, medianDat = medianDat, genesToMap = c(panel, otherGenes[i])
      )
      corMapTmp[is.na(corMapTmp)] <- -1
      higherOn <- corMapTmp[, clust] == apply(corMapTmp, 1, max)
      matchCount[i] <- sum((clustersF == clust) & (higherOn)) / sum(clustersF == clust)
      offTarget[i] <- sum((clustersF != clust) & (!higherOn)) / sum(clustersF != clust)
      totalCount[i] <- mean(c(matchCount[i], offTarget[i]))
    }
    wm <- which.max(totalCount)
    dex <- totalCount[wm] - tt
    cr <- matchCount[wm]
    tt <- totalCount[wm]
    ot <- offTarget[wm]
    addGene <- as.character(otherGenes)[wm]
    panel <- c(panel, addGene)
    onCorrect <- c(onCorrect, cr)
    offCorrect <- c(offCorrect, ot)
    dexTotal <- c(dexTotal, dex)
  }
  out <- data.frame(clust, panel, onCorrect, offCorrect, dexTotal)
  out <- out[1:dim(out)[1] - 1, ]
  out
}


#' Fraction of cells per layer
#'
#' Determines the expected proportions in each layer based on input
#'
#' @param layerIn a list corresponding to all layers of dissection for a given sample
#' @param layerNm names of all layers.  set to NULL to have this calculated
#' @param scale if TRUE (default), scale to the total number of cells
#'
#' @return vector indicating the fraction of cells in each layerNm layer
#'
#' @export
layerScale <- function(layerIn,
                       layerNm = c("L1", "L2/3", "L4", "L5", "L6"),
                       scale = TRUE) {
  if (is.null(layerNm)) {
    for (l in 1:length(layerIn)) layerNm <- c(layerNm, layerIn[[l]])
    layerNm <- sort(unique(layerNm))
  }
  total <- rep(0, length(layerNm))
  names(total) <- layerNm
  for (l in 1:length(layerIn)) total[layerIn[[l]]] <- total[layerIn[[l]]] +
      1 / length(total[layerIn[[l]]])
  if (scale) total <- total / sum(total)
  total
}


#' Layer weights per cell
#'
#' Returns a numeric vector saying how to weight a particular cell for each layer, using a smart
#'   weighting strategy
#' @param layerIn a list corresponding to all layers of dissection for a given sample
#' @param useLayer target layer
#' @param spillFactor fractional amount of cells in a layer below which it is assumed no cells are
#'   from that layer in multilayer dissection
#' @param weightCutoff anything less than this is set to 0 for convenience and to avoid rare types
#' @param layerNm names of all layers.  set to NULL to have this calculated
#'
#' @return numeric vector saying how to weight a particular cell for each layer, using a smart
#'   weighting strategy
#'
#' @export
smartLayerAllocation <- function(layerIn,
                                 useLayer = "L1",
                                 spillFactor = 0.15,
                                 weightCutoff = 0.02,
                                 layerNm = c("L1", "L2/3", "L4", "L5", "L6")) {
  if (is.null(layerNm)) {
    for (i in 1:length(layer)) layerNm <- c(layerNm, layerNm[[i]])
    layerNm <- sort(unique(layerNm))
  }
  layerMat <- matrix(0, nrow = length(layerIn), ncol = length(layerNm))
  colnames(layerMat) <- layerNm
  for (i in 1:length(layerIn)) layerMat[i, layerIn[[i]]] <- 1
  oneCount <- rowSums(layerMat) == 1
  wgtFrac <- colSums(rbind(layerMat[oneCount, ], layerMat[oneCount, ])) / 2
  wgtFrac <- wgtFrac / max(wgtFrac)
  wgtFrac[is.na(wgtFrac)] <- 0.01
  wgtFrac[wgtFrac < spillFactor] <- 1e-06 * wgtFrac[wgtFrac < spillFactor]
  wgtFrac <- wgtFrac / sum(wgtFrac)
  for (i in which(!oneCount))
    layerMat[i, ] <- (layerMat[i, ] * wgtFrac) / sum(layerMat[i, ] * wgtFrac)
  layerMat[layerMat < weightCutoff] <- 0
  for (i in 1:length(oneCount))
    layerMat[i, ] <- layerMat[i, ] / sum(layerMat[i, ])
  out <- layerMat[, useLayer]
  out
}


#' Layer weights per cell
#'
#' Returns a numeric vector saying how to weight a particular cell for each layer.  This is a
#'   wrapper for smartLayerAllocation
#'
#' @param layerIn a list corresponding to all layers of dissection for a given sample
#' @param useLayer target layer
#' @param cluster if passed the weights are smartly allocated based on laminar distributions by cluster
#' @param ... additional variables for smartLayerAllocation
#'
#' @return numeric vector with weights for cells in input layer
#'
#' @export
layerFraction <- function(layerIn,
                          useLayer = "L1",
                          cluster = NA,
                          ...) {
  weight <- rep(0, length(layerIn))
  if (is.na(cluster[1])) {
    for (l in 1:length(weight)) weight[l] <- sum(layerIn[[l]] == useLayer) / length(layerIn[[l]])
    return(weight)
  }
  for (cli in unique(cluster))
    weight[cli == cluster] <- smartLayerAllocation(layerIn[cli == cluster], useLayer, ...)
  weight[is.na(weight)] <- 0
  weight
}


#' Filter possible cluster calls using priors
#'
#' This function will return a vector of possible clusters for cells that meet a set of priors for each layer
#'
#' @param cluster vector of all clusters
#' @param layer list of layers for each cluster entry (for data sets with only laminar dissections,
#'   each list entry will be of length 1)
#' @param subsetVector a vector of TRUE/FALSE values indicated whether the entry is in the subset of
#'   interest (e.g., Cre lines); default is all
#' @param useClusters a set of clusters to be considered a priori (e.g., GABA vs. glut); default is all
#' @param rareLimit define any values less than this as 0.  The idea is to exclude rare cells
#' @param layerNm names of all layers.  set to NULL to have this calculated
#' @param scaleByLayer if TRUE, scales to the proportion of cells in each layer
#' @param scaleByFn what function should be used for the layer scaling (default=max, ignored
#'   if scaleByLayer=FALSE)
#' @param smartWeight if TRUE, multilayer dissections are weighted smartly by cluster, rather
#'   than evenly by cluster (FALSE)
#' @param spillFactor fractional amount of cells in a layer below which it is assumed no cells
#'   are from that layer in multilayer dissection
#' @param weightCutoff anything less than this is set to 0 for convenience
#'
#' @return a vector of possible clusters for cells that meet a set of priors for each layer
#'
#' @export
possibleClustersByPriors <- function(cluster,
                                     layer,
                                     subsetVector = rep(TRUE, length(cluster)),
                                     useClusters = sort(unique(cluster)),
                                     rareLimit = 0.005,
                                     layerNm = c("L1", "L2/3", "L4", "L5", "L6"),
                                     scaleByLayer = TRUE,
                                     scaleByFn = max,
                                     smartWeight = TRUE,
                                     spillFactor = 0.15,
                                     weightCutoff = 0.02) {
  if (is.null(layerNm)) {
    for (i in 1:length(layer)) layerNm <- c(layerNm, layerNm[[i]])
    layerNm <- sort(unique(layerNm))
  }
  if (!is.factor(cluster)) cluster <- factor(cluster)
  allClusters <- levels(cluster)
  out <- matrix(NA, nrow = length(allClusters), ncol = length(layerNm))
  rownames(out) <- allClusters
  colnames(out) <- layerNm
  isClust <- is.element(cluster, useClusters)
  sb <- subsetVector & isClust
  if (sum(sb) == 0) {
    return(out)
  } # Don't run on cases with no data
  subLayer <- layerScale(layer[sb], layerNm = layerNm) /
    layerScale(layer[isClust], layerNm = layerNm)
  subLayer <- subLayer / sum(subLayer)
  kpLay <- names(subLayer)[subLayer >= rareLimit]
  out[useClusters, kpLay] <- 0

  smartClust <- NA
  if (smartWeight) smartClust <- cluster

  for (lay in kpLay) {
    weight <- layerFraction(layer, lay, smartClust,
      spillFactor = spillFactor,
      weightCutoff = weightCutoff, layerNm = layerNm
    )
    kpCl <- unique(cluster[sb & (weight > 0)])
    for (cli in kpCl) out[cli, lay] <- sum(weight[sb & (cluster == cli)])
  }
  if (scaleByLayer) {
    for (lay in kpLay) out[, lay] <- out[, lay] / scaleByFn(out[, lay], na.rm = TRUE)
    out[out < rareLimit] <- 0
  }
  out
}


#' Get node height
#'
#' Returns the heights of each node, scaled from 0 (top) to 1 (leafs); this is a wrapper for dendextend functions
#'
#' @param tree a dendrogram object
#'
#' @return a vector of node heights
#'
#' @export
getNodeHeight <- function(tree) {
  nodeHeight <- get_nodes_attr(tree, "height")
  nodeHeight <- 1 - nodeHeight / max(nodeHeight)
  names(nodeHeight) <- get_nodes_attr(tree, "label")
  return(nodeHeight)
}


#' Fraction of correct calls per node
#'
#' This function returns the fraction correctly assigned to each node (as defined that the actual
#'   and predicted cluster are both in the same node)
#'
#' @param dendIn dendrogram for mapping.  Ignored if minimizeHeight=FALSE
#' @param clActual character vector of actual cluster assignments
#' @param clPredict character vector of predicted cluster assignments
#' @param minCount set to 0 results from clusters with fewer than this number of cells (default
#'   is to consider all clusters)
#' @param defaultSum value to return in cases where there are fewer than minCount cells in the
#'   actual cluster (e.g., cases that aren't considered at all)
#' @param out required for recursive function.  Do not set!
#'
#' @return matrix of two columns: (1) node name and (2) the fraction of cells in that node that
#'   are correctly assigned
#'
#' @export
fractionCorrectPerNode <- function(dendIn,
                                   clActual,
                                   clPredict,
                                   minCount = 0.1,
                                   defaultSum = -1,
                                   out = NULL) {
  clActual <- as.character(clActual)
  clPredict <- as.character(clPredict)
  if (length(dendIn) > 1) {
    for (i in 1:length(dendIn)) {
      allLabels <- labels(dendIn[[i]])
      nodeName <- attr(dendIn[[i]], "label")
      isActual <- is.element(clActual, allLabels)
      isPredict <- is.element(clPredict, allLabels)
      fractionCorrrect <- signif(sum(isActual & isPredict) / sum(isPredict + 1e-11), 3) # CHECK THIS
      if (sum(isActual) < minCount) fractionCorrrect <- defaultSum
      out <- rbind(out, c(nodeName, fractionCorrrect))
      colnames(out) <- c("nodeName", "fractionCorrrect")
      out <- fractionCorrectPerNode(dendIn[[i]], clActual, clPredict, minCount, defaultSum, out)
    }
  }
  out <- as.data.frame(out)
  out$fractionCorrrect <- as.numeric(as.character(out$fractionCorrrect))
  rownames(out) <- out[, 1]
  out
}


#' Filter by meta-data
#'
#' Return a filter of TRUE/FALSE values for a given piece of meta-data (e.g., broad class).
#'
#' @param classVector vector corresponding to the class information for filtering (e.g., vector
#'   of label calls)
#' @param sampleInfo matrix of sample information with rows corresponding to cells and columns
#'   corresponding to meta-data
#' @param classColumn column name of class information
#' @param clusterColumn column name of cluster information
#' @param threshold minimum fraction of cluster cells from a given class to be considered present
#'
#' @return a matrix of filters with rows as clusters and columns as classes with entries of TRUE or
#'   FALSE indicating whether cells from a given class can assigned to that cluster, given threshold.
#'
#' @export
filterByClass <- function(classVector,
                          sampleInfo,
                          classColumn = "cluster_type_label",
                          clusterColumn = "cluster_label",
                          threshold = 0.1) {
  ##
  out <- table(factor(sampleInfo[, clusterColumn]), factor(sampleInfo[, classColumn]))
  out <- out / rowSums(out)
  out <- out > threshold

  # Allow for names that not present in the
  # sampleInfo file
  out <- as.data.frame(out)
  tmp <- setdiff(classVector, colnames(out))
  if (length(tmp) > 0) {
    for (tm in tmp) out[, tm] <- TRUE
  }
  out <- as.matrix(out)

  # Remove '' names
  classVector[classVector == ""] <- "none"
  rownames(out)[rownames(out) == ""] <- "none"
  colnames(out)[colnames(out) == ""] <- "none"

  # Find and output the filter matrix
  outTable <- t(out[, classVector])
  if (length(classVector) == dim(sampleInfo)[1]) {
    rownames(outTable) <- sampleInfo[, "sample_id"]
  }
  outTable
}

#' Subsample cells
#'
#' Subsets a categorical vector to include up to a maximum number of values for each category.
#'
#' @param clusters vector of cluster labels (or any category) in factor or character format
#' @param subSamp maximum number of values for each category to subsample
#' @param seed for reproducibility
#'
#' @return returns a vector of TRUE / FALSE with a maximum of subSamp TRUE calls per category
#'
#' @export
subsampleCells <- function(clusters,
                           subSamp = 25,
                           seed = 5) {
  kpSamp <- rep(FALSE, length(clusters))
  for (cli in unique(as.character(clusters))) {
    set.seed(seed)
    seed <- seed + 1
    kp <- which(clusters == cli)
    kpSamp[kp[sample(1:length(kp), min(length(kp), subSamp))]] <- TRUE
  }
  kpSamp
}


#' Fraction of cells correctly assigned
#'
#' This function takes as input an ordered set of marker genes (e.g., from at iterative
#'   algorithm), and returns a vector showing the fraction of cells correctly mapped.
#'
#' @param orderedGenes an ordered list of input genes (e.g. from an iterative algorithm)
#' @param mapDat normalized data of the mapping (=reference) data set.
#' @param medianDat median value for each leaf
#' @param clustersF cluster calls for each cell
#' @param verbose whether or not to show progress in the function
#' @param plot if TRUE, plotCorrectWithGenes is run
#' @param ... parameters passed to plotCorrectWithGenes (if plot=TRUE)
#' @param return if TRUE, the value is returned
#'
#' @return a vector showing the fraction of cells correctly mapped to each cluster
#'
#' @export
fractionCorrectWithGenes <- function(orderedGenes,
                                     mapDat,
                                     medianDat,
                                     clustersF,
                                     verbose = FALSE,
                                     plot = TRUE,
                                     return = TRUE,
                                     ...) {
  numGn <- 2:length(orderedGenes)
  frac <- rep(0, length(orderedGenes))
  for (i in numGn) {
    gns <- orderedGenes[1:i]
    corMapTmp <- suppressWarnings(corTreeMapping(
      mapDat = mapDat,
      medianDat = medianDat, genesToMap = gns
    ))
    corMapTmp[is.na(corMapTmp)] <- 0
    topLeafTmp <- getTopMatch(corMapTmp)
    frac[i] <- 100 * mean(topLeafTmp[, 1] == clustersF)
  }
  frac[is.na(frac)] <- 0
  if (plot) {
    plotCorrectWithGenes(frac, genes = orderedGenes, ...)
  }
  if (return) {
    return(frac)
  }
}

#' Correct mapping at different tree heights
#'
#' This function takes as input an ordered set of marker genes (e.g., from at iterative algorithm,
#'   and returns an table showing the fraction of cells correctly mapped to a similar cell type
#'   (as defined by the heights parameter).  A height of 1 indicates correct mapping to the leaf.
#'
#' @param orderedGenes an ordered list of input genes (e.g. from an iterative algorithm)
#' @param dend dendrogram for mapping.
#' @param mapDat normalized data of the mapping (=reference) data set.
#' @param medianDat median value for each leaf
#' @param clustersF cluster calls for each cell
#' @param minVal minimum number of genes to consider from the list in the mapping
#' @param heights height in the tree to look at
#' @param verbose whether or not to show progress in the function
#'
#' @return a matrix of fractions of cells correctly mapped for different tree heights (columns)
#'   and different gene panels (rows)
#'
#' @export
buildQualityTable <- function(orderedGenes,
                              dend,
                              mapDat,
                              medianDat,
                              clustersF,
                              minVal = 2,
                              heights = c((0:100) / 100),
                              verbose = FALSE) {
  minVal <- max(2, round(minVal))
  nodeHeight <- get_nodes_attr(dend, "height")
  nodeHeight <- 1 - nodeHeight / max(nodeHeight)
  names(nodeHeight) <- get_nodes_attr(dend, "label")

  outTable <- matrix(0, nrow = length(orderedGenes), ncol = length(heights))
  rownames(outTable) <- orderedGenes
  colnames(outTable) <- paste0("Height_", heights)

  for (r in minVal:length(orderedGenes)) {
    if (verbose) print(r)
    gns <- orderedGenes[1:r]
    topLeafTmp <- suppressWarnings(getTopMatch(corTreeMapping(
      mapDat = mapDat, medianDat = medianDat, genesToMap = gns
    )))
    lcaVector <- nodeHeight[lca(dend, as.character(topLeafTmp[, 1]), clustersF)]
    for (c in 1:length(heights)) {
      outTable[r, c] <- mean(lcaVector >= (1 - heights[c]))
    }
  }
  outTable
}

#' Plot fraction correct
#'
#' This function is a wrapper for plot designd for plotting the fraction correctly mapped for a given gene set.
#'   If geneN is the Nth gene, the plotted value indicates correct mapping using genes 1:N.
#'
#' @param frac a numeric vector indicating the fraction of cells correctly mapped for a given gene panel
#' @param genes ordered character vector (e.g., of genes) to be plotted; default is names(frac)
#' @param ... additional parameters for plot.
#'
#' @export
plotCorrectWithGenes <- function(frac,
                                 genes = names(frac),
                                 xlab = "Number of genes in panel",
                                 main = "All clusters gene panel",
                                 ylim = c(-10, 100),
                                 lwd = 5,
                                 ylab = "Percent of nuclei correctly mapping",
                                 colLine = "grey",
                                 ...) {
  numGn <- 1:length(frac)
  plot(numGn, frac,
    type = "l", col = "grey", xlab = xlab,
    ylab = ylab, main = main, ylim = ylim, lwd = lwd, ...
  )
  abline(h = (-2:20) * 5, lty = "dotted", col = colLine)
  abline(h = 0, col = "black", lwd = 2)
  text(numGn, frac, genes, srt = 90)
}


#' Get binary (aka beta) score
#'
#' Returns a beta score which indicates the binaryness of a gene across clusters.  High scores
#'   (near 1) indicate that a gene is either on or off in nearly all cells of every cluster.
#'   Scores near 0 indicate a cells is non-binary (e.g., not expressed, ubiquitous, or
#'   randomly expressed).  This value is used for gene filtering prior to defining clustering.
#'
#' @param propExpr a matrix of proportions of cells (rows) in a given cluster (columns) with
#'   CPM/FPKM > 1 (or 0, HCT uses 1)
#' @param returnScore if TRUE returns the score, if FALSE returns the ranks
#' @param spec.exp scaling factor (recommended to leave as default)
#'
#' @return returns a numeric vector of beta score (or ranks)
#'
#' @export
getBetaScore <- function(propExpr,
                         returnScore = TRUE,
                         spec.exp = 2) {
  calc_beta <- function(y, spec.exp = 2) {
    d1 <- as.matrix(dist(y))
    eps1 <- 1e-10
    # Marker score is combination of specificity and sparsity
    score1 <- sum(d1^spec.exp) / (sum(d1) + eps1)
    score1
  }

  betaScore <- apply(propExpr, 1, calc_beta)
  betaScore[is.na(betaScore)] <- 0
  if (returnScore) return(betaScore)
  scoreRank <- rank(-betaScore)
  scoreRank
}
