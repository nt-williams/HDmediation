suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(dplyr)
    library(purrr)
    library(foreach)
    library(doFuture)
})

source("_research/transported/gendata4.R")

# Fit with partial TMLE or not
tmle <- TRUE

source("_research/SL.lightgbm.R")
# source("_research/SL.glm.saturated.R")
# source("_research/SL.glmnet3.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

learners <- c("SL.earth", "SL.lightgbm", "SL.glm.interaction", "SL.glm", "SL.mean")
folds <- 5

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- sample_data(n)

    z <- names(dat)[startsWith(names(dat), "z")]
    m <- names(dat)[startsWith(names(dat), "m")]

    mediation(dat, "a", "w", "z", "m", "y", "S",
              family = "binomial",
              folds = folds,
              partial_tmle = tmle,
              learners_g = "SL.mean",
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
