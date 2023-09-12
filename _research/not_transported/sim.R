suppressPackageStartupMessages({
    library(HDmediation)
    library(glue)
    library(tidyverse)
})

# gendata.R is multivariate M and Z
# gendata2.R is binary M and Z
dgp <- 1
if (dgp == 1) {
    source("_research/not_transported/gendata3.R")
} else {
    source("_research/not_transported/gendata2.R")
}

# Fit with partial TMLE or not
tmle <- T

# source("_research/SL.lightgbm.R")
source("_research/SL.glm.saturated.R")
source("_research/SL.glmnet3.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

if (dgp == 1) {
    learners <- c("SL.glm.interaction", "SL.glm.saturated", "SL.glmnet3", "SL.glm")
    folds <- 5
} else {
    learners <- "SL.glm.saturated"
    folds <- 1
}

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- gendata(n)
    
    z <- names(dat)[startsWith(names(dat), "Z")]
    m <- names(dat)[startsWith(names(dat), "M")]
    
    mediation(dat, "A", "W1", 
              z, m, "Y", 
              family = "binomial", 
              folds = folds, partial_tmle = tmle,
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

res <- mutate(res, 
              n = case_when(
                  n == 1 ~ 500, 
                  n == 2 ~ 1000, 
                  n == 3 ~ 5000, 
                  n == 4 ~ 1e4
              ))

saveRDS(res, glue("_research/data/sim_not_transported_{tmle}_{dgp}_{id}.rds"))
