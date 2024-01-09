suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(dplyr)
    library(purrr)
    library(foreach)
    library(doFuture)
    library(mlr3extralearners)
})

source("_research/not_transported/gendata4.R")

# Fit with partial TMLE or not
tmle <- T

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

learners <- list("mean", "glm",
                 list("earth", degree = 3, fast.k = 0, fast.beta = 0, id = "earth"),
                 list("lightgbm", min_data_in_leaf = 5),
                 list("ranger", num.trees = 500, id = "ranger1"),
                 list("ranger", num.trees = 1000, id = "ranger2"), 
                 list("nnet", trace = FALSE, decay = 0.05), 
                 list("nnet", trace = FALSE, decay = 0.1), 
                 list("nnet", trace = FALSE, decay = 0.01))

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- sample_data(n)

    z <- names(dat)[startsWith(names(dat), "z")]
    m <- names(dat)[startsWith(names(dat), "m")]
    
    folds <- ifelse(n > 5000, 2, 5)

    mediation(dat, "a", "w", "z", "m", "y",
              family = "binomial",
              folds = folds,
              partial_tmle = tmle,
              learners_g = "glm",
              learners_c = learners,
              learners_b = learners,
              learners_e = learners,
              learners_hz = learners,
              learners_u = learners,
              learners_ubar = learners,
              learners_v = learners,
              learners_vbar = learners)
}, .id = "n")

res <- mutate(res, n = case_when(
    n == 1 ~ 500,
    n == 2 ~ 1000,
    n == 3 ~ 5000,
    n == 4 ~ 1e4
))

saveRDS(res, glue("_research/data/sim_not_transported_{tmle}_cont_{id}.rds"))
