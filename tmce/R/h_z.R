h_z <- function(data, npsem, pred_m, folds, learners, ...) {
    h_z <- matrix(nrow = nrow(data), ncol = 1)
    for (v in seq_along(folds)) {
        train <- stack_data(origami::training(data, folds[[v]]), npsem, pred_m)
        valid <- origami::validation(data, folds[[v]])
        
        p_zmw <- crossfit(train, list(valid), "tmp_tmce_stack_indicator", 
                          c(npsem$Z, npsem$M, npsem$W), "binomial", learners)[[1]]
        
        p_mw <- crossfit(train, list(valid), "tmp_tmce_stack_indicator", 
                         c(npsem$M, npsem$W), "binomial", learners)[[1]]
        
        h_z[folds[[v]]$validation_set, 1] <- p_zmw * (1 - p_mw) / (1 - p_zmw) * p_mw
    } 
    h_z
}

stack_data <- function(data, npsem, pred_m) {
    delta <- data
    delta[[npsem$M]] <- pred_m
    
    out <- rbind(data, delta)
    out[["tmp_tmce_stack_indicator"]] <- rep(c(0, 1), each = nrow(data))
    out
}
