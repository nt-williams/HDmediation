g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz <- function(z, a, w) {
    prob1 <- plogis(-log(2) + (log(10) * a) - log(2) * w[, "W1"])
    z * prob1 + (1 - z) * (1 - prob1)
}

pm <- function(m, z, a, w) {
    prob1 <- plogis(-log(2) + log(12) * z - log(1.4) * w[, "W1"])
    m * prob1 + (1 - m) * (1 - prob1)
}

pmaw <- function(m, a, w) {
    pm(m, 1, a, w) * pz(1, a, w) +
        pm(m, 0, a, w) * pz(0, a, w)
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
    plogis(-log(5) + log(8) * z + log(10) * m -
               log(1.2) * w[, "W1"] + log(1.2) * w[, "W1"] * z)
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
    w1 <- rbinom(N, 1, .4)
    w <- data.frame(W1 = w1)
    a <- rbinom(N, 1, g(1))
    z <- rbinom(N, 1, pz(1, a, w))
    m <- rbinom(N, 1, pm(1, z, a, w))
    y <- rbinom(N, 1, my(m, z, a, w))
    data.frame(W1 = w1, A = a, Z = z, M = m, Y = y)
}

# h <- pmaw(m, astar, w) / pm(m, z, aprime, w)

If <- function(dat, aprime, astar) {
    w <- dat[, "W1", drop = F]
    
    ipwy <- (dat$A == aprime) / g(aprime)
    hm <- pmaw(dat$M, astar, w) / pm(dat$M, dat$Z, aprime, w)
    eify <- ipwy * hm * (dat$Y - my(dat$M, dat$Z, aprime, w))
    
    ipwz <- ipwy
    eifz <- ipwz * (u(dat$Z, w, aprime, astar) - intu(w, aprime, astar))
    
    ipwm <- (dat$A == astar) / g(astar)
    vbar <- intv(1, w, aprime) * pmaw(1, astar, w) + intv(0, w, aprime) * pmaw(0, astar, w)
    eifm <- ipwm * (intv(dat$M, w, aprime) - vbar)
    
    eify + eifz + eifm + vbar
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

    c("11" = weighted.mean(v_11, prob_w),
      "10" = weighted.mean(v_10, prob_w),
      "00" = weighted.mean(v_00, prob_w),
      "indirect" = weighted.mean(v_11 - v_10, prob_w),
      "direct" = weighted.mean(v_10 - v_00, prob_w))
}
