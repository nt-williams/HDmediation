gen_data <- function(n = 1e3) {
    w0 <- rbinom(n, 1, .5)
    w1 <- rbinom(n, 1, 0.5)
    
    probsite <- plogis(log(1.2) * w1)
    site <- rbinom(n, 1, probsite)
    
    a <- rbinom(n, 1, plogis(log(2) * w1))
    z <- rbinom(n, 1, plogis(-log(2) + log(4) * a - log(2) * w1 + log(1.4) * site))
    m <- rbinom(n, 1, plogis(-log(2) + log(4) * a + log(2) * z - log(1.4) * w1 + log(1.4) * site))
    y <- rbinom(n, 1, plogis(-log(5) + log(3) * a + log(1.6) * a * w1 + log(8) * z + log(4) * m - log(1.2) * w1))
    
    data.frame(
        W0 = w0,
        W1 = w1,
        A = a,
        Z = z,
        M = m,
        Y = y,
        S = site
    )
}

ps <- function(w){
    # plogis(log(1.2) * w[, 'W1'] + log(1.2) * w[, 'W0'])
    plogis(log(1.2) * w[, 'W1'])
}

g <- function(a, w) {
    # pscore <- plogis(log(2) * w[, 'W0'])
    pscore <- plogis(log(2) * w[, 'W1'])
    a * pscore + (1 - a) * (1 - pscore)
}

pz <- function(z, a, w, s) {
    prob1 <- plogis(-log(2) + log(4) * a - log(2) * w[, "W1"] + log(1.4) * s)
    z * prob1 + (1 - z) * (1 - prob1)
}

pm <- function(m, z, a, w, s) {
    prob1 <- plogis(-log(2) + log(4) * a + log(2) * z - log(1.4) * w[, "W1"] + log(1.4) * s)
    m * prob1 + (1 - m) * (1 - prob1)
}

pmaw <- function(m, a, w, s) {
    pm(m, 1, a, w, s) * pz(1, a, w, s) + pm(m, 0, a, w, s ) * pz(0, a, w, s)
}

pmw <- function(m, w, s) {
    pmaw(m, 1, w, s) * g(1, w) + pmaw(m, 0, w, s) * g(0, w)
}

r <- function(z, a, m, w, s) {
    pm(m, z, a, w, s) * pz(z, a, w, s) / pmaw(m, a, w, s)
}

e <- function(a, m, w, s) {
    pmaw(m, a, w, s) * g(a, w) / pmw(m, w, s)
}

my <- function(m, z, a, w) {
    plogis(-log(5) + log(3) * a + log(1.6) * a * w[, "W1"] + log(8) * z + log(4) * m - log(1.2) * w[, "W1"])
}

u <- function(z, w, aprime, astar) {
    my(1, z, aprime, w) * pmaw(1, astar, w, 0) + my(0, z, aprime, w) * pmaw(0, astar, w, 0)
}

pzmars <- function(z, a, w) {
    pz(z, a, w, 1) * ps(w) + pz(z, a, w, 0) * (1 - ps(w))
}

psazw <- function(z, a, w) {
    pz(z, a, w, 1) * ps(w) / pzmars(z, a, w)
}

pmmars <- function(m, z, a, w) {
    pm(m, z, a, w, 1) * psazw(z, a, w) +
        pm(m, z, a, w, 0) * (1 - psazw(z, a, w))
}

b <- function(a, z, m, w) {
    pm(m, z, a, w, 1) / pmmars(m, z, a, w) *
        pz(z, a, w, 1) / pzmars(z, a, w) * ps(w)
}

intu <- function(w, aprime, astar) {
    u(1, w, aprime, astar) * pz(1, aprime, w, 0) +
        u(0, w, aprime, astar) * pz(0, aprime, w, 0)
}

intv <- function(m, w, aprime) {
    my(m, 1, aprime, w) * pz(1, aprime, w, 0) +
        my(m, 0, aprime, w) * pz(0, aprime, w, 0)
}

compute_eif <- function(dat, aprime, astar) {
    a <- dat$A
    z <- dat$Z
    m <- dat$M
    y <- dat$Y
    s <- dat$S
    w <- dat[, c("W0", "W1")]
    
    v <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)
    h <- pmaw(m, astar, w, 0) / pm(m, z, aprime, w, 0) * (1 - b(aprime, z, m, w)) / b(aprime, z, m, w)
    
    eifyt <- (s == 1) / mean(1 - s) * (a == aprime) / g(aprime, w)  * h * (y - my(m, z, aprime, w))
    eifut <- (s == 0) / mean(1 - s) * (a == aprime) / g(aprime, w) * (u(z, w, aprime, astar) - intu(w, aprime, astar))
    eifvt <- (s == 0) / mean(1 - s) * (a == astar) / g(astar, w) * (intv(m, w, aprime) - v)
    
    list(theta = mean(v[s == 0]), 
         eif = eifyt + eifut + eifvt + (s == 0) / mean(1 - s) * (v - mean(v[s == 0])))
}

compute_eif2 <- function(dat, aprime, astar) {
    a <- dat$A
    z <- dat$Z
    m <- dat$M
    y <- dat$Y
    s <- dat$S
    w <- dat[, c("W0", "W1")]
    
    v <- intv(1, w, aprime) * pmaw(1, astar, w, 0) + intv(0, w, aprime) * pmaw(0, astar, w, 0)
    h <- pmaw(m, astar, w, 0) / pm(m, z, aprime, w, 0) * (1 - b(aprime, z, m, w)) / b(aprime, z, m, w)
    
    eifyt <- (s == 1) / mean(1 - s) * (a == aprime) / g(aprime, w)  * h * (y - my(m, z, aprime, w))
    eifut <- (s == 0) / mean(1 - s) * (a == aprime) / g(aprime, w) * (u(z, w, aprime, astar) - intu(w, aprime, astar))
    eifvt <- (s == 0) / mean(1 - s) * (a == astar) / g(astar, w) * (intv(m, w, aprime) - v)
    eifwt <- (s == 0) / mean(1 - s) * (v - mean(v[s == 0]))
    list(
        y = eifyt,
        z = eifut, 
        m = eifvt, 
        w = eifwt
    )
}

truth_mediation <- function(dat) {
    # compute influence function with convenience functions
    eif11 <- compute_eif(dat, 1, 1)
    eif10 <- compute_eif(dat, 1, 0)
    eif00 <- compute_eif(dat, 0, 0)
    
    # compute parameter estimate and efficiency bound
    indirect <- (eif11$theta - eif10$theta) + (mean(eif11$eif) - mean(eif10$eif))
    direct   <- (eif10$theta - eif00$theta) + (mean(eif10$eif) - mean(eif00$eif))
    
    eif_indirect <- eif11$eif - eif10$eif
    eif_direct   <- eif10$eif - eif00$eif
    
    var_indirect <- var(eif_indirect)
    var_direct <- var(eif_direct)
    
    data.frame(
        parameter = c("eif11", "eif10", "eif00", "indirect", "direct"),
        truth = c(eif11$theta, eif10$theta, eif00$theta, indirect, direct),
        eff_bound = c(var(eif11$eif), var(eif10$eif), var(eif00$eif), var_indirect, var_direct)
    )
}
