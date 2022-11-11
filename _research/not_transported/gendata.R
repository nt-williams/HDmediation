g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz1 <- function(z, a, w) {
    prob1 <- plogis(-log(2) + (log(10) * a) - log(2) * w[, "W1"])
    z * prob1 + (1 - z) * (1 - prob1)
}

pz2 <- function(z, a, w) {
    prob1 <- plogis(-log(2) + (log(5) * a) - log(0.55) * w[, "W1"])
    z * prob1 + (1 - z) * (1 - prob1)
}

pz <- function(z1, z2, a, w) {
    pz1(z1, a, w) * pz2(z2, a, w)
}

pm1 <- function(m, z1, a, w) {
    prob1 <- plogis(-log(2) + log(7) * z1 - log(1.4) * w[, "W1"])
    m * prob1 + (1 - m) * (1 - prob1)
}

pm2 <- function(m, z2, a, w) {
    prob1 <- plogis(0.05 + log(2.5) * z2 - log(1.4) * w[, "W1"])
    m * prob1 + (1 - m) * (1 - prob1)
}

pm1aw <- function(m1, a, w) {
    pm1(m1, 1, a, w) * pz1(1, a, w) + 
        pm1(m1, 0, a, w) * pz1(0, a, w)
}

pm2aw <- function(m2, a, w) {
    pm2(m2, 1, a, w) * pz2(1, a, w) + 
        pm2(m2, 0, a, w) * pz2(0, a, w)
}

pmaw <- function(m1, m2, a, w) {
    pm1aw(m1, a, w) * pm2aw(m2, a, w)
}

pm <- function(m1, m2, z1, z2, a, w) {
    pm1(m1, z1, a, w) * pm2(m2, z2, a, w)
}

pmw <- function(m1, m2, w) {
    pmaw(m1, m2, 1, w) * g(1) + pmaw(m1, m2, 0, w) * g(0)
}

r <- function(z1, z2, a, m1, m2, w) {
    pm(m1, m2, z1, z2, a, w) * pz(z1, z2, a, w) / pmaw(m1, m2, a, w)
}

e <- function(a, m1, m2, w) {
    pmaw(m1, m2, a, w) * g(a) / pmw(m1, m2, w)
}

my <- function(m1, m2, z1, z2, a, w) {
    plogis(-log(5) + log(8) * z1  + log(4) * m1 - 
               log(1.2) * w[, "W1"] - log(2) * z2 + log(1.2) * m2 + 
               log(1.2) * w[, "W1"] * z1)
}

# u <- function(z, w, aprime, astar) {
#     my(1, z, aprime, w) * pmaw(1, astar, w) + my(0, z, aprime, w) * pmaw(0, astar, w)
# }
# 
# intu <- function(w, aprime, astar) {
#     u(1, w, aprime, astar) * pz(1, aprime, w) +
#         u(0, w, aprime, astar) * pz(0, aprime, w)
# }

intv <- function(m1, m2, w, aprime) {
    my(m1, m2, 1, 1, aprime, w) * pz1(1, aprime, w) * pz2(1, aprime, w) +
        my(m1, m2, 1, 0, aprime, w) * pz1(1, aprime, w) * pz2(0, aprime, w) +
        my(m1, m2, 0, 1, aprime, w) * pz1(0, aprime, w) * pz2(1, aprime, w) +
        my(m1, m2, 0, 0, aprime, w) * pz1(0, aprime, w) * pz2(0, aprime, w)
}

gendata <- function(N) {
    w1 <- rbinom(N, 1, .4)
    w <- data.frame(W1 = w1)
    a <- rbinom(N, 1, g(1))
    z1 <- rbinom(N, 1, pz1(1, a, w))
    z2 <- rbinom(N, 1, pz2(1, a, w))
    m1 <- rbinom(N, 1, pm1(1, z1, a, w))
    m2 <- rbinom(N, 1, pm2(1, z2, a, w))
    y <- rbinom(N, 1, my(m1, m2, z1, z2, a, w))
    data.frame(W1 = w1, A = a, Z1 = z1, Z2 = z2, M1 = m1, M2 = m2, Y = y)
}

# h <- pmaw(m1, m2, astar, w) / pm(m1, m2, z1, z2, aprime, w)

truth <- function() {
    w <- expand.grid(W1 = c(0, 1))
    
    prob_w <- vector("numeric", nrow(w))
    for (i in 1:nrow(w)) {
        w1 <- w[i, "W1"] * 0.4 + (1 - w[i, "W1"]) * 0.6
        prob_w[i] <- w1
    }

    aprime <- astar <- 1
    v_11 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w) + 
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w) + 
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w) + 
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w)
    
    astar <- 0
    v_10 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w) + 
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w) + 
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w) + 
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w)
    
    aprime <- 0
    v_00 <- intv(1, 1, w, aprime) * pmaw(1, 1, astar, w) + 
        intv(1, 0, w, aprime) * pmaw(1, 0, astar, w) + 
        intv(0, 1, w, aprime) * pmaw(0, 1, astar, w) + 
        intv(0, 0, w, aprime) * pmaw(0, 0, astar, w)
    
    c("11" = weighted.mean(v_11, prob_w), 
      "10" = weighted.mean(v_10, prob_w), 
      "00" = weighted.mean(v_00, prob_w), 
      "indirect" = weighted.mean(v_11 - v_10, prob_w),
      "direct" = weighted.mean(v_10 - v_00, prob_w))
}
