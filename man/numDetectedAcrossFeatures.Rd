% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/numDetectedAcrossFeatures.R
\name{numDetectedAcrossFeatures}
\alias{numDetectedAcrossFeatures}
\alias{numDetectedAcrossFeatures,ANY-method}
\alias{numDetectedAcrossFeatures,SummarizedExperiment-method}
\title{Number of detected expression values per group of features}
\usage{
numDetectedAcrossFeatures(x, ...)

\S4method{numDetectedAcrossFeatures}{ANY}(
  x,
  ids,
  detection_limit = 0,
  subset_row = NULL,
  subset_col = NULL,
  average = FALSE,
  BPPARAM = SerialParam()
)

\S4method{numDetectedAcrossFeatures}{SummarizedExperiment}(x, ..., exprs_values = "counts")
}
\arguments{
\item{x}{A numeric matrix of counts containing features in rows and cells in columns.
Alternatively, a \linkS4class{SummarizedExperiment} object containing such a count matrix.}

\item{...}{For the generic, further arguments to pass to specific methods.

For the SummarizedExperiment method, further arguments to pass to the ANY method.

For the ANY method, further arguments to pass to the \code{\link{nexprs}} function.}

\item{ids}{A factor of length \code{nrow(x)}, specifying the set to which each feature in \code{x} belongs.

Alternatively, a list of integer or character vectors, where each vector specifies the indices or names of features in a set.}

\item{detection_limit}{Numeric scalar providing the value above which  observations are deemed to be expressed.}

\item{subset_row}{An integer, logical or character vector specifying the features to use.
Defaults to all features.}

\item{subset_col}{An integer, logical or character vector specifying the cells to use.
Defaults to all cells with non-\code{NA} entries of \code{ids}.}

\item{average}{Logical scalar indicating whether the proportion of non-zero counts in each group should be computed instead.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying whether summation should be parallelized.}

\item{exprs_values}{A string or integer scalar specifying the assay of \code{x} containing the matrix of counts
(or any other expression quantity that can be meaningfully summed).}
}
\value{
An integer or numeric matrix containing the number of detected expression values in each group of features (row) and cell (column).
}
\description{
Computes the number of detected expression values (default defined as non-zero counts) for each cell in each group of features.
}
\examples{
example_sce <- mockSCE()

ids <- sample(paste0("GENE_", 1:100), nrow(example_sce), replace=TRUE)
byrow <- numDetectedAcrossFeatures(example_sce, ids)
head(byrow[,1:10])

}
\seealso{
\code{\link{nexprs}}, on which this function is based.
}
\author{
Aaron Lun
}
