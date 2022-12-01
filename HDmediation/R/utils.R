gl <- glue::glue

scale_y <- function(y, family, bounds) {
    if (family == "binomial") {
        return(list(y = y, bounds = NULL))
    }

    bounds <- bounds_y(y, bounds)
    list(y = (y - bounds[1]) / (bounds[2] - bounds[1]),
         bounds = bounds)
}

bounds_y <- function(y, bounds) {
    if (is.null(bounds)) {
        return(c(min(y, na.rm = T), max(y, na.rm = T)))
    }
    c(bounds[1], bounds[2])
}

rescale_y <- function(y, bounds) {
    if (is.null(bounds)) {
        return(y)
    }
    (y*(bounds[2] - bounds[1])) + bounds[1]
}
