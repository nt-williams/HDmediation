suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(dplyr)
    library(purrr)
    library(foreach)
    library(doFuture)
    library(mlr3extralearners)
})

source("_research/not_transported/gendata5.R")

# Fit with partial TMLE or not
tmle <- FALSE

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

learners <- list("mean", "glm",
                 list("earth", degree = 3, fast.k = 0, fast.beta = 0, id = "earth"),
                 list("lightgbm",
                      num_iterations = 4,
                      learning_rate = 0.40,
                      max_depth = 3,
                      min_data_in_leaf = 25,
                      id = "lgb1"),
                 list("lightgbm",
                      num_iterations = 5,
                      learning_rate = 0.02,
                      max_depth = 17,
                      min_data_in_leaf = 15,
                      id = "lgb2"),
                 list("ranger",
                      mtry.ratio = 0.01405632,
                      replace = FALSE,
                      sample.fraction = 0.3146349,
                      num.trees = 261,
                      id = "ranger1"),
                 list("ranger",
                      mtry.ratio = 0.2329873,
                      replace = TRUE,
                      sample.fraction = 0.309394,
                      num.trees = 711),
                 list("nnet", trace = FALSE, decay = 0.05), 
                 list("nnet", trace = FALSE, decay = 0.1), 
                 list("nnet", trace = FALSE, decay = 0.01))

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- sample_data(n)
    
    z <- names(dat)[startsWith(names(dat), "z")]
    m <- names(dat)[startsWith(names(dat), "m")]
    
    folds <- ifelse(n > 5000, 2, 5)
    
    mediation(dat, "a", "w", z, m, "y",
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

saveRDS(res, glue("_research/data/sim_not_transported_dgpRealistic_{tmle}_{id}.rds"))
