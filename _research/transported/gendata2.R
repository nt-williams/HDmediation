g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

pz1 <- function(z, a, w, s) {
    prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"] + log(1.4)*s)
    z * prob1 + (1 - z) * (1 - prob1)
}

pz2 <- function(z, a, w, s) {
    prob1 <- plogis(-(log(1.5) * a) - log(0.55) * w[, "W1"] + log(0.5)*s)
    z * prob1 + (1 - z) * (1 - prob1)
}

pm1 <- function(m, z1, a, w, s) {
    prob1 <- plogis(-log(2) + log(4) * z1 - log(1.4) * w[, "W1"] + log(0.3)*s)
    m * prob1 + (1 - m) * (1 - prob1)
}

pm2 <- function(m, z2, a, w, s) {
    prob1 <- plogis(0.05 + log(1.8) * z2 - log(1.4) * w[, "W1"] + log(0.1)*s)
    m * prob1 + (1 - m) * (1 - prob1)
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

# pmw <- function(m, w) {
#     pmaw(m, 1, w) * g(1) + pmaw(m, 0, w) * g(0)
# }
#
# r <- function(z, a, m, w) {
#     pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
# }
#
# e <- function(a, m, w) {
#     pmaw(m, a, w) * g(a) / pmw(m, w)
# }

my <- function(m1, m2, z1, z2, a, w) {
    plogis(-log(5) + log(8) * z1  + log(4) * m1 -
               log(1.2) * w[, "W1"] - log(2) * z1 + log(1.2) * m2 +
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
