#' @export
#' @rdname weightQCAcrossComparison
#' @examples 
#' ids <- rep(LETTERS[1:2], times = 1000)
#' groups <- rep(1:2, each = 1000)
#' pmito <- c(rnorm(1000, 6.8, 1), rnorm(1000, 7.3, 1))
#' 
#' qc_weight <- weightQCAcrossComparison(pmito, ids, groups)
#' 
weightQCAcrossComparison <- function(x, ids, groups, subset=NULL, min.cells=10) 
{
  .weight_qc_across_comparison(x, ids=ids, groups=groups, subset=subset, min.cells=10)
}



.weight_qc_across_comparison <- function(x, ids, groups, subset=NULL, min.cells=10) 
{
  if (length(unique(groups)) == 1) {
    stop("Need at least two groups")
  }
  
  ids <- .process_ids_with_groups(x, ids=ids, groups=groups, subset=subset, min.cells=min.cells)
  
  lost <- is.na(ids)
  subset <- if (any(lost)) which(!lost)
  
  if (!.noOpSubset(subset, length(x))) {
    x <- x[subset]
  }
  
  sub.ids <- ids[!lost]
  sub.groups <- groups[!lost]
  
  xl <- lapply(unique(sub.ids), function(id) x[sub.ids == id])
  groupsl <- lapply(unique(sub.ids), function(id) sub.groups[sub.ids == id])
  
  weights <- mapply(FUN =.weight_qc_across_comparison_by_id, xl, groupsl)
  
  # unroll weights so that same order as x
  inds <- unlist(lapply(unique(sub.ids), function(id) seq_along(x)[sub.ids == id]))
  unlist(weights)[inds]
}

.weight_qc_across_comparison_by_id <- function(x, groups) {
  
  # devide x into quantiles of 5%
  qs <- quantile(x, probs = seq(0, 1, .05))
  xcut <- cut(x, breaks = qs, include.lowest = TRUE)
  
  # number of observations per group in each quantile
  gnames <- unique(groups)
  tab <- lapply(gnames, function(g) tabulate(xcut[groups == g]))
  
  # weights are min observations in bin divided by actual observations
  mins <- do.call(pmin, tab)
  wbins <- lapply(tab, function(x) mins/x)
  names(wbins) <- gnames
  
  weights <- lapply(gnames, function(g) {
    gweights <- xcut[groups == g]
    levels(gweights) <- wbins[[g]]
    gweights
  })
  
  # unroll weights so that same order as x
  inds <- unlist(lapply(gnames, function(g) seq_along(x)[groups == g]))
  weights <- unlist(weights)[inds]
  as.numeric(as.character(weights))
}


.process_ids_with_groups <- function(x, ids, groups, subset, min.cells = 10) {    
  if (.has_multi_ids(ids)) {
    ids <- .df_to_factor(ids)
  } 
  if (length(x)!=length(ids)) {
    stop("length of 'ids' and length of 'x' are not equal")
  }
  if (length(unique(groups)) == 1) {
    stop("need at least two 'groups'")
  }
  
  tab <- table(ids, groups)
  keep <- apply(tab, 1, function(row) any(row >= min.cells))
  ndisc <- sum(!keep)
  keep <- names(keep)[keep]
  
  if (!length(keep)) {
    stop("need at least 'min.cells' in one 'ids' and 'groups' combination")
  }
  if (ndisc > 0) {
    message(ndisc, " clusters have less than 'min.cells' in one group and will not be used.")
    ids[!ids %in% keep] <- NA_integer_
  }
  
  if (!is.null(subset)) {
    ids[!ids %in% subset] <- NA_integer_
  }
  ids
}
