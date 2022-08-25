source("_research/sim/gen_data.R")

dat <- gen_data(1e+05, 86746)

local({
    source("_research/sim/truth.R", local = TRUE)

    a <- dat$A
    z <- dat$Z
    m <- dat$M
    y <- dat$Y
    s <- dat$S
    w <- dat[, c("W0", "W1")]
    
    aprime <- astar <- 1
    
    v <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)
    h <- pmaw(m, astar, w, 0) / pm(m, z, aprime, w, 0) * (1 - b(aprime, z, m, w)) / b(aprime, z, m, w)
    
    eifyt <<- (s == 1) / mean(1 - s) * (a == aprime) / g(aprime)  * h * (y - my(m, z, aprime, w))
    eifut <<- (s == 0) / mean(1 - s) * (a == aprime) / g(aprime) * (u(z, w, aprime, astar) - intu(w, aprime, astar))
    eifvt <<- (s == 0) / mean(1 - s) * (a == astar) / g(astar) * (intv(m, w, aprime) - v)
    # eifyt + eifut + eifvt + (s == 0) / mean(1 - s) * v
})

devtools::load_all("tmce")

A <- "A"
S <- "S"
W <- c("W0", "W1")
Z <- "Z"
M <- "M"
Y <- "Y"

npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
folds <- make_folds(dat, 1)

gg <- g(dat, npsem, folds)
ee <- e(dat, npsem, folds)
cc <- see(dat, npsem, folds)
bb <- b(dat, npsem, "binomial", folds)
hz <- h_z(dat, npsem, folds)
t <- 1 - mean(dat[[npsem$S]])

aprime <- 0
astar <- 0

hm <- h_m(hz, gg, ee, aprime, astar)
uu <- u(dat, npsem, bb, hm, aprime, folds)
uubar <- ubar(dat, npsem, uu, aprime, folds)
vv <- v(dat, npsem, bb, hz, aprime, folds)
vvbar <- vbar(dat, npsem, vv, astar, folds)

bb[, gl("b({aprime},Z,M,W)")] 

local({
    source("_research/sim/truth.R", local = TRUE)
    
    a <- dat$A
    z <- dat$Z
    m <- dat$M
    y <- dat$Y
    s <- dat$S
    w <- dat[, c("W0", "W1")]
    
    aprime <- 0
    astar <- 0

    head(pz(z, aprime, w, 0) / r(z, aprime, m, w, 0), 20)
})

# EIF calculation
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

