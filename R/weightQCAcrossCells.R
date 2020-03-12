#' @export
#' @rdname sumCountsAcrossCells
#' @importFrom SingleCellExperiment altExp altExps altExp<- altExps<-
#' reducedDimNames reducedDim<- reducedDim reducedDims<- reducedDims
#' @examples 
#' example_sce <- mockSCE()
#' example_sce <- addPerCellQC(example_sce, subsets=list(Mito=1:10))
#' ids <- rep(LETTERS[1:5], each = ncol(example_sce)/5)
#' groups <- rep(1:2, times = ncol(example_sce)/2)
#' 
#' qc_weight <- weightQCAcrossCells(example_sce, ids, groups, 'subsets_Mito_percent')
#' 
setMethod("weightQCAcrossComparison", "SingleCellExperiment", function(x, ids, groups, qc_metric,
                                                                  ..., subset_row=NULL, subset_col=NULL,
                                                                  min.cells=10)
{
  y <- .weight_qc_across_comparison(x, ids=ids, groups=groups, qc_metric=qc_metric, subset_row=subset_row, subset_col=subset_col,
                               min.cells=10)
  
  y
})




#' @importFrom SummarizedExperiment SummarizedExperiment
.weight_qc_across_comparison <- function(x, ids, groups, qc_metric, subset_row=NULL, subset_col=NULL, min.cells=10) 
{
  browser()
  new.ids <- .process_ids_with_groups(x=x, ids=ids, groups=groups, subset_col=subet_col, min.cells=min.cells)
  
  lost <- is.na(ids)
  subset_col <- if (any(lost)) which(!lost)
  
  if (!.noOpSubset(subset_row, nrow(x))) {
    x <- x[subset_row,,drop=FALSE]
  }
  if (!.noOpSubset(subset_col, ncol(x))) {
    x <- x[,subset_col,drop=FALSE]
  }
  
  if (length(unique(groups)) == 1) {
    stop("Need at least two groups")
  }
  
  
  
  new.groups <- .process_groups(x, ids=new.ids, groups=groups)
  
  
  
  # Avoid additional parallelization from DA methods.
  oldBP <- getAutoBPPARAM()
  setAutoBPPARAM(SerialParam()) 
  on.exit(setAutoBPPARAM(oldBP))
  
  sub.ids <- ids[!lost]
  out <- bplapply(by.core, FUN=.colsum, group=sub.ids, BPPARAM=BPPARAM)
  out <- do.call(rbind, out)
  
  freq <- table(sub.ids)
  freq <- as.integer(freq[colnames(out)])
  if (average) {
    out <- t(t(out)/freq)
  }
  
  list(mat=out, freq=freq)
  
  
  #####
  
  
  sum.out <- .sum_across_cells(x, ids=new.ids, subset_row=subset_row,
                               average=average, BPPARAM=BPPARAM, modifier=modifier)
  
  mat <- sum.out$mat 
  mapping <- match(colnames(mat), as.character(new.ids))
  coldata <- .create_coldata(original.ids=ids, mapping=mapping, 
                             freq=sum.out$freq, store_number=store_number)
  
  # non-NULL coldata determines the output SE column names,
  # so we make sure they're sync'd to something sensible.
  if (.has_multi_ids(ids)) {
    rownames(coldata) <- colnames(mat) <- NULL
  } else {
    rownames(coldata) <- colnames(mat)
  }
  
  output <- list(mat)
  names(output) <- if(average) {
    "average" 
  } else { 
    "sum" 
  }
  SummarizedExperiment(output, colData=coldata)
}
