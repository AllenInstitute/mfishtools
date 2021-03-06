% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markerGenesAndMapping.r
\name{leafToNodeMedians}
\alias{leafToNodeMedians}
\title{Return mean node expression}
\usage{
leafToNodeMedians(dend, medianDat, branches = getBranchList(dend), fnIn = mean)
}
\arguments{
\item{dend}{dendrogram for mapping.  Ignored if medianDat is passed}

\item{medianDat}{median expression data at each node}

\item{branches}{a particular format of branch information from the dendrogram structure}

\item{fnIn}{function to use to wrap up to the node level (default = mean)}
}
\value{
a matrix of mean node expression (rows=genes, columns=nodes)
}
\description{
Define expression at a node as the MEAN expression for each leaf as default (using the
  median removes all specific marker genes!)
}
