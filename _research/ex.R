source("../_research/gen_data.R")

dat <- gen_data(1000)

local({
    source("../_research/truth.R", local = TRUE)
    truth <<- truth_mediation(1e5)
})

devtools::load_all()

A <- "A"
S <- "S"
W <- c("W0", "W1")
Z <- "Z"
M <- "M"
Y <- "Y"

aprime <- 1
astar <- 0

npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
folds <- make_folds(dat, 1)

learners <- c("SL.glm", "SL.xgboost")

# works
gg <- g(dat, npsem, folds, learners)
ee <- e(dat, npsem, folds, learners)
cc <- see(dat, npsem, folds, learners)
bb <- b(dat, npsem, "binomial", folds, learners)
hz <- h_z(dat, npsem, folds, learners)
hm <- h_m(hz, gg, ee, aprime, astar)
uu <- u(dat, npsem, bb, hm, aprime, folds, learners)
uubar <- ubar(dat, npsem, uu, aprime, folds, learners)
vv <- v(dat, npsem, bb, hz, aprime, folds, learners)
vvbar <- vbar(dat, npsem, vv, astar, folds, learners)
t <- 1 - mean(dat[[npsem$S]])

# test
D_PY <-
    ((dat[[npsem$S]] == 1) & (dat[[npsem$A]] == aprime)) /
    (t * gg[, gl("g({aprime}|w)")]) *
    (1 - cc[, 1]) /
    cc[, 1] *
    hm * (dat[[npsem$Y]] - bb[, gl("b({aprime},Z,M,W)")])

D_PZ <-
    ((dat[[npsem$S]] == 0) & (dat[[npsem$A]] == aprime)) /
    (t * gg[, gl("g({aprime}|w)")]) *
    (uu[, 1] - uubar[, 1])

D_PM <-
    ((dat[[npsem$S]] == 0) & (dat[[npsem$A]] == astar)) /
    (t * gg[, gl("g({astar}|w)")]) *
    (vv[, 1] - vvbar[, 1])

D_PW <- ((dat[[npsem$S]] == 0) / t) * vvbar[, 1]

D_P <- D_PY + D_PZ + D_PM + D_PW
theta <- mean(D_P)
