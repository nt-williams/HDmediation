library(HDmediation)
library(glue)
library(purrr)
library(dplyr)

# gendata.R is multivariate M and Z
# gendata2.R is binary M and Z
dgp <- 2
if (dgp == 1) {
    source("_research/transported/gendata3.R")
} else {
    source("_research/transported/gendata2.R")
}

# Fit with partial TMLE or not
tmle <- TRUE

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

if (dgp == 1) {
    learners <- c("glm", "sat_glm", "cv_sat_glmnet")
} else {
    learners <- c("sat_glm")
    folds = 1
}

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- gendata(n)
    
    if (!exists("folds")) {
        folds <- case_when(n == 500 ~ 10,
                           n == 1000 ~ 10,
                           n == 5000 ~ 5,
                           n == 1e4  ~ 4)
    }

    z <- names(dat)[startsWith(names(dat), "Z")]
    m <- names(dat)[startsWith(names(dat), "M")]
    
    mediation(dat, "A", "W1", 
              z, m, "Y", "S",
              family = "binomial", 
              folds = folds, partial_tmle = tmle,
              learners_g = c("mean"), 
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

saveRDS(res, glue("_research/data/sim_transported_{tmle}_{dgp}_{id}.rds"))
