#' @export
#' @rdname weightQCAcrossComparison
#' @examples 
#' ids <- rep(LETTERS[1:2], times = 1000)
#' groups <- rep(1:2, each = 1000)
#' pmito <- c(rnorm(1000, 6.8, 1), rnorm(1000, 7.3, 1))
#' 
#' qc_weights <- qcWeightsAcrossComparison(pmito, ids, groups)
#' 
qcWeightsAcrossComparison <- function(x, ids, groups, subset=NULL, min.cells=10) 
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
  
  ids <- ids[!lost]
  groups <- groups[!lost]
  
  weights <- numeric(length(x))
  for (id in unique(ids)) {
    is.id <- ids == id
    x.id <- x[is.id]
    groups.id <- groups[is.id]
    
    weights[is.id] <- .qc_weights_across_comparison_by_id(x.id, groups.id)
  }
  return(weights)
}


.qc_weights_across_comparison_by_id <- function(x, groups) {
  
  # devide x into quantiles of 5%
  qs <- quantile(x, probs = seq(0, 1, .05))
  xcut <- cut(x, breaks = qs, include.lowest = TRUE)
  
  # fraction of observations per group in each quantile
  gnames <- unique(groups)
  gfracs <- lapply(gnames, function(g) {
    is.g <- groups == g
    tabulate(xcut[is.g]) / sum(is.g)
  })
  
  # weights are min fraction in quantile divided by observed fractions
  min.fracs <- do.call(pmin, gfracs)
  bin.weights <- lapply(gfracs, function(fracs) min.fracs/fracs)
  names(bin.weights) <- gnames
  
  # expand quantile weights to each observation
  weights <- numeric(length(x))
  for (g in gnames) {
    is.g <- groups == g
    gweights <- xcut[is.g]
    levels(gweights) <- bin.weights[[g]]
    weights[is.g] <- as.numeric(as.character(gweights))
  }
  return(weights)
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
