source("_research/sim/gen_data.R")

devtools::load_all("tmce")

dat <- gen_data(5000)

A <- "A"
S <- "S"
W <- c("W0", "W1")
Z <- "Z"
M <- "M"
Y <- "Y"

tmce::tmce(dat, A, S, W, Z, M, Y, "binomial", folds = 1)

local({
    source("_research/sim/truth.R", local = TRUE)
    truth <<- truth_mediation(gendata(3e6))
})
