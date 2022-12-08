g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz1 <- function(z, a, w, s) {
    prob1 <- 0.25 + 0.1*a + 0.2*w[, "W1"] - 0.05*s
    z * prob1 + (1 - z) * (1 - prob1)
}

pz2 <- function(z, a, w, s) {
    prob1 <- 0.4 + 0.1*a - 0.1*w[, "W1"] + 0.075*s
    z * prob1 + (1 - z) * (1 - prob1)
}

pm1 <- function(m, z1, a, w, s) {
    prob1 <- 0.6 + 0.1*z1 + 0.05*a - 0.3*w[, "W1"]
    m * prob1 + (1 - m) * (1 - prob1)
}

pm2 <- function(m, z2, a, w, s) {
    prob1 <- 0.33 + 0.22*z2 + 0.05*a + 0.15*w[, "W1"] - 0.05*s
    m * prob1 + (1 - m) * (1 - prob1)
}

pz <- function(z1, z2, a, w, s) {
    pz1(z1, a, w, s) * pz2(z2, a, w, s)
}

pm <- function(m1, m2, z1, z2, a, w, s) {
    pm1(m1, z1, a, w, s) * pm2(m2, z2, a, w, s)
}

pm1aw <- function(m1, a, w, s) {
    pm1(m1, 1, a, w, s) * pz1(1, a, w, s) +
        pm1(m1, 0, a, w, s) * pz1(0, a, w, s)
}

pm2aw <- function(m2, a, w, s) {
    pm2(m2, 1, a, w, s) * pz2(1, a, w, s) +
        pm2(m2, 0, a, w, s) * pz2(0, a, w, s)
}

pmaw <- function(m1, m2, a, w, s) {
    pm1aw(m1, a, w, s) * pm2aw(m2, a, w, s)
}

pmw <- function(m1, m2, w) {
    pmaw(m1, m2, 1, w, 0) * g(1) + pmaw(m1, m2, 0, w, 0) * g(0)
}

r <- function(z1, z2, a, m1, m2, w) {
    pm(m1, m2, z1, z2, a, w, 0) * pz(z1, z2, a, w, 0) / pmaw(m1, m2, a, w, 0)
}

e <- function(a, m1, m2, w) {
    pmaw(m1, m2, a, w, 0) * g(a) / pmw(m1, m2, w)
}

my <- function(m1, m2, z1, z2, a, w) {
    plogis(-log(5) + log(8) * z1  + log(4) * m1 -
               log(1.2) * w[, "W1"] - log(2) * z2 + log(1.2) * m2 +
               log(1.2) * w[, "W1"] * z1)
}

u <- function(z1, z2, w, aprime, astar) {
    my(1, 1, z1, z2, aprime, w) * pmaw(1, 1, astar, w, 0) +
        my(0, 1, z1, z2, aprime, w) * pmaw(0, 1, astar, w, 0) +
        my(0, 0, z1, z2, aprime, w) * pmaw(0, 0, astar, w, 0) +
        my(1, 0, z1, z2, aprime, w) * pmaw(1, 0, astar, w, 0)
}

intu <- function(w, aprime, astar) {
    u(1, 1, w, aprime, astar) * pz(1, 1, aprime, w, 0) +
        u(1, 0, w, aprime, astar) * pz(1, 0, aprime, w, 0) +
        u(0, 0, w, aprime, astar) * pz(0, 0, aprime, w, 0) +
        u(0, 1, w, aprime, astar) * pz(0, 1, aprime, w, 0)
}

intv <- function(m1, m2, w, aprime) {
    my(m1, m2, 1, 1, aprime, w) * pz1(1, aprime, w, 0) * pz2(1, aprime, w, 0) +
        my(m1, m2, 1, 0, aprime, w) * pz1(1, aprime, w, 0) * pz2(0, aprime, w, 0) +
        my(m1, m2, 0, 1, aprime, w) * pz1(0, aprime, w, 0) * pz2(1, aprime, w, 0) +
        my(m1, m2, 0, 0, aprime, w) * pz1(0, aprime, w, 0) * pz2(0, aprime, w, 0)
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

pzmars <- function(z1, z2, a, w) {
    pz(z1, z2, a, w, 1) * psw(1, w$W1) + pz(z1, z2, a, w, 0) * psw(0, w$W1)
}

psazw <- function(z1, z2, a, w) {
    pz(z1, z2, a, w, 1) * psw(1, w$W1) / pzmars(z1, z2, a, w)
}

pmmars <- function(m1, m2, z1, z2, a, w) {
    pm(m1, m2, z1, z2, a, w, 1) * psazw(z1, z2, a, w) +
        pm(m1, m2, z1, z2, a, w, 0) * (1 - psazw(z1, z2, a, w))
}

psazmw <- function(a, z1, z2, m1, m2, w) {
    pm(m1, m2, z1, z2, a, w, 1) / pmmars(m1, m2, z1, z2, a, w) *
        pz(z1, z2, a, w, 1) / pzmars(z1, z2, a, w) * psw(1, w$W1)
}

gendata <- function(N) {
    w1 <- rbinom(N, 1, .4)
    w <- data.frame(W1 = w1)
    s <- rbinom(N, 1, psw(1, w1))
    a <- rbinom(N, 1, g(1))
    z1 <- rbinom(N, 1, pz1(1, a, w, s))
    z2 <- rbinom(N, 1, pz2(1, a, w, s))
    m1 <- rbinom(N, 1, pm1(1, z1, a, w, s))
    m2 <- rbinom(N, 1, pm2(1, z2, a, w, s))
    y <- rbinom(N, 1, my(m1, m2, z1, z2, a, w))
    data.frame(S = s, W1 = w1, A = a, Z1 = z1, Z2 = z2, M1 = m1, M2 = m2, Y = y)
}

If <- function(dat, aprime, astar) {
    w <- dat[, "W1", drop = F]
    
    ipwy <- (dat$A == aprime & dat$S == 1) / g(aprime) * ps(0)
    hm <- pmaw(dat$M1, dat$M2, astar, w, 0) / pm(dat$M1, dat$M2, dat$Z1, dat$Z2, aprime, w, 0)
    cs <- ((1 - psazmw(aprime, dat$Z1, dat$Z2, dat$M1, dat$M2, w)) / 
               psazmw(aprime, dat$Z1, dat$Z2, dat$M1, dat$M2, w))
    
    eify <- ipwy * hm * cs * (dat$Y - my(dat$M1, dat$M2, dat$Z1, dat$Z2, aprime, w))
    
    ipwz <- (dat$A == aprime & dat$S == 0) / g(aprime) * ps(0)
    eifz <- ipwz * (u(dat$Z1, dat$Z2, w, aprime, astar) - intu(w, aprime, astar))
    
    ipwm <- (dat$A == astar & dat$S == 0) / g(astar) * ps(0)
    vbar <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w, 0) + 
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w, 0) +
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w, 0) + 
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w, 0)
    eifm <- ipwm * (intv(dat$M1, dat$M2, w, aprime) - vbar)
    
    eifw <- (dat$S == 0) / ps(0) * (vbar - mean(vbar[dat$S == 0]))
    
    eify + eifz + eifm + eifw
}

truth <- function() {
    w <- expand.grid(W1 = c(0, 1))

    # prob(W|S=0)
    prob_ws0 <- psw(0, w$W1) * pw(w$W1) / ps(0)

    aprime <- astar <- 1
    v_11 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w, 0) +
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w, 0) +
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w, 0) +
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w, 0)

    astar <- 0
    v_10 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w, 0) +
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w, 0) +
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w, 0) +
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w, 0)

    aprime <- 0
    v_00 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w, 0) +
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w, 0) +
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w, 0) +
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w, 0)

    c("indirect" = weighted.mean(v_11 - v_10, prob_ws0),
      "direct" = weighted.mean(v_10 - v_00, prob_ws0))
}
