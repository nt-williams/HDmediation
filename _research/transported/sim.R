source("_research/gen_data.R")
source("_research/SL.lightgbm.R")

library(furrr)
library(glue)
n <- 500

# id <- Sys.getenv("SGE_TASK_ID")
# 
# if (id == "undefined" || id == "") id <- 1

simulate <- function(n, folds = 1) {
    dat <- gen_data(n)
    
    A <- "A"
    S <- "S"
    W <- "W1"
    Z <- "Z"
    M <- "M"
    Y <- "Y"
    
    psi <- tmce::tmce(dat, A, S, W, Z, M, Y, "binomial", folds = folds)
    
    data.frame(n = n,
               direct = psi$direct,
               var_direct = psi$var_direct,
               indirect = psi$indirect,
               var_indirect = psi$var_indirect)

    # write.csv(
    #     data.frame(
    #         n = n,
    #         direct = psi$direct,
    #         var_direct = psi$var_direct,
    #         indirect = psi$indirect,
    #         var_indirect = psi$var_indirect
    #     ),
    #     glue::glue("_research/sim/data/{id}-{n}.csv"),
    #     row.names = FALSE
    # )
}

# args <- commandArgs(trailingOnly = TRUE)

# args <- list(
#     10000
# )

plan(multisession, workers = 5)

res <- future_map_dfr(1:500, simulate, n = n)
saveRDS(res, glue("_research/data/sim_{n}.rds"))

plan(sequential)
# simulate(as.numeric(args[[1]]), 1)
# 
# quit("no")
