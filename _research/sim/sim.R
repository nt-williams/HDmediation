source("_research/sim/gen_data.R")

id <- Sys.getenv("SGE_TASK_ID")

if (id == "undefined" || id == "") id <- 1

simulate <- function(n, folds) {
    dat <- gen_data(n)
    
    A <- "A"
    S <- "S"
    W <- c("W0", "W1")
    Z <- "Z"
    M <- "M"
    Y <- "Y"
    
    psi <- tmce::tmce(dat, A, S, W, Z, M, Y, "binomial", folds = folds)

    write.csv(
        data.frame(
            seed = seed,
            n = n,
            direct = psi$direct,
            var_direct = psi$var_direct,
            indirect = psi$indirect,
            var_indirect = psi$var_indirect
        ),
        glue::glue("_research/sim/data/{id}-{n}.csv"),
        row.names = FALSE
    )
}

args <- commandArgs(trailingOnly = TRUE)

# args <- list(
#     10000
# )

simulate(as.numeric(args[[1]]), 5)

quit("no")
