suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(dplyr)
    library(purrr)
    library(foreach)
    library(doFuture)
})

source("_research/not_transported/gendata5.R")

# Fit with partial TMLE or not
tmle <- T

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

folds <- 2

dat <- sample_data(1e4)

z <- names(dat)[startsWith(names(dat), "z")]
m <- names(dat)[startsWith(names(dat), "m")]

learners <- list("mean", "glm", 
                 list("earth", degree = 3, fast.k = 0, fast.beta = 0, id = "earth"),
                 list("lightgbm", min_data_in_leaf = 5),
                 list("ranger", num.trees = 500, id = "ranger1"),
                 list("ranger", num.trees = 1000, id = "ranger2"))

res <- mediation(dat, "a", "w", z, m, "y",
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

saveRDS(res, glue("_research/data/sim_not_transported_dgpRealistic_{tmle}_{id}.rds"))
