g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz <- function(z, a, w, s) {
    prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"] + log(1.4)*s)
    z * prob1 + (1 - z) * (1 - prob1)
}

pm <- function(m, z, a, w, s) {
    prob1 <- plogis(-log(2) + log(10) * z - log(1.4) * w[, "W1"] + log(0.3)*s)
    m * prob1 + (1 - m) * (1 - prob1)
}

pmaw <- function(m, a, w, s) {
    pm(m, 1, a, w, s) * pz(1, a, w, s) +
        pm(m, 0, a, w, s) * pz(0, a, w, s)
}

pmw <- function(m, w) {
    pmaw(m, 1, w, 0) * g(1) + pmaw(m, 0, w, 0) * g(0)
}

r <- function(z, a, m, w) {
    pm(m, z, a, w, 0) * pz(z, a, w, 0) / pmaw(m, a, w, 0)
}

e <- function(a, m, w) {
    pmaw(m, a, w, 0) * g(a) / pmw(m, w)
}

my <- function(m, z, a, w) {
    plogis(-log(5) + log(8) * z + log(6) * m -
               log(1.2) * w[, "W1"] + log(1.2) * w[, "W1"] * z)
}

u <- function(z, w, aprime, astar) {
    my(1, z, aprime, w) * pmaw(1, astar, w, 0) + my(0, z, aprime, w) * pmaw(0, astar, w, 0)
}

intu <- function(w, aprime, astar) {
    u(1, w, aprime, astar) * pz(1, aprime, w, 0) +
        u(0, w, aprime, astar) * pz(0, aprime, w, 0)
}

intv <- function(m, w, aprime) {
    my(m, 1, aprime, w) * pz(1, aprime, w, 0) +
        my(m, 0, aprime, w) * pz(0, aprime, w, 0)
}

pw <- function(w) {
    prob1 <- 0.4
    w * prob1 + (1 - w) * (1 - prob1)
}

psw <- function(s, w) {
    prob1 <- plogis(log(1.2) * w)
    s * prob1 + (1 - s) * (1 - prob1)
}

ps <- function(s) {
    psw(s, 1) * pw(1) + psw(s, 0) * pw(0)
}

pzmars <- function(z, a, w) {
    pz(z, a, w, 1) * psw(1, w$W1) + pz(z, a, w, 0) * psw(0, w$W1)
}

psazw <- function(z, a, w) {
    pz(z, a, w, 1) * psw(1, w$W1) / pzmars(z, a, w)
}

pmmars <- function(m, z, a, w) {
    pm(m, z, a, w, 1) * psazw(z, a, w) +
        pm(m, z, a, w, 0) * (1 - psazw(z, a, w))
}

psazmw <- function(a, z, m, w) {
    pm(m, z, a, w, 1) / pmmars(m, z, a, w) *
        pz(z, a, w, 1) / pzmars(z, a, w) * psw(1, w$W1)
}

gendata <- function(N) {
    w1 <- rbinom(N, 1, .4)
    w <- data.frame(W1 = w1)
    s <- rbinom(N, 1, psw(1, w1))
    a <- rbinom(N, 1, g(1))
    z <- rbinom(N, 1, pz(1, a, w, s))
    m <- rbinom(N, 1, pm(1, z, a, w, s))
    y <- rbinom(N, 1, my(m, z, a, w))
    data.frame(S = s, W1 = w1, A = a, Z = z, M = m, Y = y)
}

If <- function(dat, aprime, astar) {
    w <- dat[, "W1", drop = F]
    
    ipwy <- (dat$A == aprime & dat$S == 1) / g(aprime) * ps(0)
    hm <- pmaw(dat$M, astar, w, 0) / pm(dat$M, dat$Z, aprime, w, 0)
    cs <- ((1 - psazmw(aprime, dat$Z, dat$M, w)) / psazmw(aprime, dat$Z, dat$M, w))
    
    eify <- ipwy * hm * cs * (dat$Y - my(dat$M, dat$Z, aprime, w))
    
    ipwz <- (dat$A == aprime & dat$S == 0) / g(aprime) * ps(0)
    eifz <- ipwz * (u(dat$Z, w, aprime, astar) - intu(w, aprime, astar))
    
    ipwm <- (dat$A == astar & dat$S == 0) / g(astar) * ps(0)
    vbar <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)
    eifm <- ipwm * (intv(dat$M, w, aprime) - vbar)
    
    eifw <- (dat$S == 0) / ps(0) * (vbar - mean(vbar[dat$S == 0]))
    
    eify + eifz + eifm + eifw
}

truth <- function() {
    w <- expand.grid(W1 = c(0, 1))

    # prob(W|S=0)
    prob_ws0 <- psw(0, w$W1) * pw(w$W1) / ps(0)

    aprime <- astar <- 1
    v_11 <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)

    astar <- 0
    v_10 <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)

    aprime <- 0
    v_00 <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)

    c("indirect" = weighted.mean(v_11 - v_10, prob_ws0),
      "direct" = weighted.mean(v_10 - v_00, prob_ws0))
}
