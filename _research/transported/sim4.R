# READ ME -----------------------------------------------------------------
#
#       Author: Nick Williams
#       Created: 2023-08-24
#
# -------------------------------------------------------------------------

# Fit with partial TMLE or not
tmle <- TRUE

library(HDmediation)
library(glue)
library(dplyr)
library(purrr)
library(future)

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

source("_research/transported/gendata4.R")

learners <- list("mean", "glm", "earth",
                 list("lightgbm",
                      min_data_in_leaf = 1,
                      num_iterations = 500,
                      learning_rate = 0.05,
                      max_bin = 10,
                      id = "lgbm1"),
                 list("lightgbm",
                      min_data_in_leaf = 1,
                      num_iterations = 1000,
                      learning_rate = 0.025,
                      max_bin = 10,
                      id = "lgbm2"),
                 list("ranger", num.trees = 500, id = "ranger1"),
                 list("ranger", num.trees = 1000, id = "ranger2"))

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- sample_data(n)

    folds <- case_when(n == 500 ~ 10,
                       n == 1000 ~ 10,
                       n == 5000 ~ 4,
                       n == 1e4  ~ 2)

    mediation(dat, "a", "w", "z", "m", "y", "s",
              family = "binomial",
              folds = folds,
              partial_tmle = tmle,
              learners_g = learners,
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

saveRDS(res, glue("_research/data/sim_transported_{tmle}_cont_{id}.rds"))
