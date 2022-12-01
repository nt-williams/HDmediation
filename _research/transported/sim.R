suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(tidyverse)
})

# gendata.R is multivariate M and Z
# gendata2.R is binary M and Z
dgp <- 1
if (dgp == 1) {
    source("_research/transported/gendata.R")
} else {
    source("_research/transported/gendata2.R")
}

# source("_research/SL.lightgbm.R")
# source("_research/SL.glmnet3.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

SL.glm.saturated <- function(Y, X, newX, family, obsWeights, ...) {
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    f <- as.formula(paste0("Y ~ .^", ncol(X)))
    fit.glm <- glm(f, data = X, family = family, weights = obsWeights)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}

learners <- c("SL.glm.saturated")
# learners <- c("SL.glm", "SL.glm.interaction", "SL.mean", "SL.lightgbm", "SL.glmnet3", "SL.earth")

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- gendata(n)
    
    # folds <- case_when(n <= 1000 ~ 10,
    #                    n == 5000 ~ 4,
    #                    TRUE ~ 2)
    folds <- 1
    z <- names(dat)[startsWith(names(dat), "Z")]
    m <- names(dat)[startsWith(names(dat), "M")]
    
    mediation(dat, "A", "W1", 
              z, m, "Y", "S",
              family = "binomial", 
              folds = folds, 
              learners_g = "SL.mean", 
              learners_b = learners, 
              learners_e = learners, 
              learners_hz = learners, 
              learners_u = learners, 
              learners_ubar = learners, 
              learners_v = learners, 
              learners_vbar = learners)
}, .id = "n")

res <- mutate(res, 
              n = case_when(
                  n == 1 ~ 500, 
                  n == 2 ~ 1000, 
                  n == 3 ~ 5000, 
                  n == 4 ~ 1e4
              ))

saveRDS(res, glue("_research/data/sim_transported_sat_{dgp}_{id}.rds"))
