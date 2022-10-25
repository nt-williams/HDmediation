library(HDmediation)
library(glue)

source("_research/SL.lightgbm.R")

g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz <- function(z, a, w) {
    prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"])
    z * prob1 + (1 - z) * (1 - prob1)
}

pm <- function(m, z, a, w) {
    prob1 <- plogis(-log(2) + log(4) * z - log(1.4) * w[, "W1"])
    m * prob1 + (1 - m) * (1 - prob1)
}

pmaw <- function(m, a, w) {
    pm(m, 1, a, w) * pz(1, a, w) + pm(m, 0, a, w) * pz(0, a, w)
}

pmw <- function(m, w) {
    pmaw(m, 1, w) * g(1) + pmaw(m, 0, w) * g(0)
}

r <- function(z, a, m, w) {
    pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

e <- function(a, m, w) {
    pmaw(m, a, w) * g(a) / pmw(m, w)
}

my <- function(m, z, a, w) {
    plogis(-log(5) + log(8) * z  + log(4) * m - log(1.2) * w[, "W1"] + log(1.2) * w[, "W1"] * z)
}

u <- function(z, w, aprime, astar) {
    my(1, z, aprime, w) * pmaw(1, astar, w) + my(0, z, aprime, w) * pmaw(0, astar, w)
}

intu <- function(w, aprime, astar) {
    u(1, w, aprime, astar) * pz(1, aprime, w) +
        u(0, w, aprime, astar) * pz(0, aprime, w)
}

intv <- function(m, w, aprime) {
    my(m, 1, aprime, w) * pz(1, aprime, w) +
        my(m, 0, aprime, w) * pz(0, aprime, w)
}

gendata <- function(N) {
    w0 <- rbinom(N, 1, .5)
    w1 <- rbinom(N, 1, .4)
    w <- data.frame(W1 = w1, W0 = w0)
    a <- rbinom(N, 1, g(1))
    z <- rbinom(N, 1, pz(1, a, w))
    m <- rbinom(N, 1, pm(1, z, a, w))
    y <- rbinom(N, 1, my(m, z, a, w))
    data.frame(W0 = w0, W1 = w1, A = a, Z = z, M = m, Y = y)
}

truth <- function() {
    w <- expand.grid(W1 = c(0, 1))

    prob_w <- vector("numeric", nrow(w))
    for (i in 1:nrow(w)) {
        w1 <- w[i, "W1"] * 0.4 + (1 - w[i, "W1"]) * 0.6
        prob_w[i] <- w1
    }

    aprime <- astar <- 1
    v_11 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)

    astar <- 0
    v_10 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)

    aprime <- 0
    v_00 <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)

    c("indirect" = weighted.mean(v_11 - v_10, prob_w),
      "direct" = weighted.mean(v_10 - v_00, prob_w))
}

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

dat <- gendata(5000)
res <- mediation(dat, "A", c("W0", "W1"), "Z", "M", "Y", S = NULL,
                 family = "binomial", folds = 1, 
                 learners_hz = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.lightgbm", "SL.earth"))

saveRDS(res, glue("_research/data/sim_not_transported_{id}.rds"))
