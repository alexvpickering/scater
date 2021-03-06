.incorporate_common_vis_row <- function(df, se, mode, colour_by, size_by, shape_by, 
    by_exprs_values, by_show_single, other_fields, multiplier=NULL) 
{
    colour_by_out <- retrieveFeatureInfo(se, colour_by, exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    if (!is.null(multiplier)) {
        colour_by_out$value<- colour_by_out$value[multiplier]
    }
    df$colour_by <- colour_by_out$value

    shape_by_out <- retrieveFeatureInfo(se, shape_by, exprs_values = by_exprs_values)
    shape_by <- shape_by_out$name
    if (!is.null(multiplier)) {
        shape_by_out$value<- shape_by_out$value[multiplier]
    }
    df$shape_by <- .coerce_to_factor(shape_by_out$value, 10, "shape_by")

    size_by_out <- retrieveFeatureInfo(se, size_by, exprs_values = by_exprs_values)
    size_by <- size_by_out$name
    if (!is.null(multiplier)) {
        size_by_out$value<- size_by_out$value[multiplier]
    }
    df$size_by <- size_by_out$value

    for (o in other_fields) {
        other <- retrieveFeatureInfo(se, o, exprs_values=by_exprs_values)
        if (!is.null(multiplier)) {
            other$value<- other$value[multiplier]
        }
        df <- .add_other_or_warn(df, other)
    }

    list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by)
}

.incorporate_common_vis_col <- function(df, se, mode, colour_by, size_by, shape_by, 
    by_exprs_values, by_show_single, other_fields, multiplier=NULL) 
{
    colour_by_out <- retrieveCellInfo(se, colour_by, exprs_values = by_exprs_values)
    colour_by <- colour_by_out$name
    if (!is.null(multiplier)) {
        colour_by_out$value<- colour_by_out$value[multiplier]
    }
    df$colour_by <- colour_by_out$value

    shape_by_out <- retrieveCellInfo(se, shape_by, exprs_values = by_exprs_values)
    shape_by <- shape_by_out$name
    if (!is.null(multiplier)) {
        shape_by_out$value<- shape_by_out$value[multiplier]
    }
    df$shape_by <- .coerce_to_factor(shape_by_out$value, 10, "shape_by")

    size_by_out <- retrieveCellInfo(se, size_by, exprs_values = by_exprs_values)
    size_by <- size_by_out$name
    if (!is.null(multiplier)) {
        size_by_out$value <- size_by_out$value[multiplier]
    }
    df$size_by <- size_by_out$value

    for (o in other_fields) {
        other <- retrieveCellInfo(se, o, exprs_values=by_exprs_values)
        if (!is.null(multiplier)) {
            other$value<- other$value[multiplier]
        }
        df <- .add_other_or_warn(df, other)
    }

    list(df=df, colour_by = colour_by, shape_by = shape_by, size_by = size_by)
}

.add_other_or_warn <- function(df, other) {
    if (!is.null(other$name)) {
        if (!other$name %in% colnames(df)) {
            df[[other$name]] <- other$value
        } else {
            warning(sprintf("not adding duplicated '%s' from 'other_fields'", other$name))
        }
    }
    df
}

################################################
## Creating pair plots.

.makePairs <- function(data_matrix) 
# with thanks to Gaston Sanchez, who posted this code online
# https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
{
    if ( is.null(names(data_matrix)) )
        names(data_matrix) <- paste0("row", 1:nrow(data_matrix))
    exp_grid <- expand.grid(x = 1:ncol(data_matrix), y = 1:ncol(data_matrix))
    exp_grid <- exp_grid[exp_grid$x != exp_grid$y,]
    all_panels <- do.call("rbind", lapply(1:nrow(exp_grid), function(i) {
        xcol <- exp_grid[i, "x"]
        ycol <- exp_grid[i, "y"]
        data.frame(xvar = names(data_matrix)[ycol], yvar = names(data_matrix)[xcol],
                   x = data_matrix[, xcol], y = data_matrix[, ycol], data_matrix)
    }))
    all_panels$xvar <- factor(all_panels$xvar, levels = names(data_matrix))
    all_panels$yvar <- factor(all_panels$yvar, levels = names(data_matrix))
    densities <- do.call("rbind", lapply(1:ncol(data_matrix), function(i) {
        data.frame(xvar = names(data_matrix)[i], yvar = names(data_matrix)[i],
                   x = data_matrix[, i])
    }))
    list(all = all_panels, densities = densities)
}
