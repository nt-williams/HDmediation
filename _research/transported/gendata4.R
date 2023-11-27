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

g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

dz <- function(z, a, w, s) {
    mu <- 0.25 + 0.1*a + 0.2*w + 0.05*s
    dnorm(z, mean = mu)
}

dm <- function(m, z, a, w, s) {
    mu <- 0.6 + 0.1*z + 0.05*a - 0.3*w + 0.075*s
    dnorm(m, mean = mu)
}

my <- function(m, z, a, w) {
    plogis(-log(5) + log(8) * z + log(6) * m - log(1.2) * w + log(3) * w * z)
}

dmaw <- function(m, a, w, s) {
    f <- function(z) {
        dm(m, z, a, w, s) * dz(z, a, w, s)
    }
    integrate(f, -Inf, Inf)$value
}

dmw <- function(m, w) {
    pma1w <- dmaw(m, 1, w, 0)
    pma0w <- dmaw(m, 0, w, 0)
    pma1w * g(1) + pma0w * g(0)
}

e <- function(a, m, w) {
    g(a) * dmaw(m, a, w, 0) / dmw(m, w)
}

r <- function(z, a, m, w) {
    dz(z, a, w, 0) * dm(m, z, a, w, 0) / dmaw(m, a, w, 0)
}

ipw <- function(A, ap, as, s, z, m, w) {
    ((as.numeric(A == ap) & (s == 1)) / (ps(0) * g(ap))) *
        psw(0, w) / psw(1, w) *
        g(ap) / g(as) *
        dz(z, ap, w, 0) / r(z, ap, m, w) *
        e(as, m, w) / e(ap, m, w)
}

sample_data <- function(n) {
    w <- rbinom(n, 1, pw(1))
    s <- rbinom(n, 1, psw(1, w))
    a <- rbinom(n, 1, g(1))
    z <- (0.25 + 0.1*a + 0.2*w + 0.05*s) + rnorm(n)
    m <- (0.6 + 0.1*z + 0.05*a - 0.3*w + 0.075*s) + rnorm(n)
    y <- rbinom(n, 1, my(m, z, a, w))
    data.table::data.table(S = as.numeric(s),
                           w = as.numeric(w),
                           a = as.numeric(a),
                           z = as.numeric(z),
                           m = as.numeric(m),
                           y = as.numeric(y))
}
