pw <- function(w) {
    prob1 <- 0.4
    w * prob1 + (1 - w) * (1 - prob1)
}

g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

dz <- function(z, a, w) {
    mu <- 0.25 + 0.1*a + 0.2*w
    dnorm(z, mean = mu)
}

dm <- function(m, z, a, w) {
    mu <- 0.6 + 0.1*z + 0.05*a - 0.3*w
    dnorm(m, mean = mu)
}

my <- function(m, z, a, w) {
    plogis(-log(5) + log(8) * z + log(6) * m - log(1.2) * w + log(3) * w * z)
}

dmaw <- function(m, a, w) {
    f <- function(z) {
        dm(m, z, a, w) * dz(z, a, w)
    }
    integrate(f, -Inf, Inf)$value
}

dmw <- function(m, w) {
    pma1w <- dmaw(m, 1, w)
    pma0w <- dmaw(m, 0, w)
    pma1w * g(1) + pma0w * g(0)
}

e <- function(a, m, w) {
    g(a) * dmaw(m, a, w) / dmw(m, w)
}

r <- function(z, a, m, w) {
    dz(z, a, w) * dm(m, z, a, w) / dmaw(m, a, w)
}

ipw <- function(A, ap, as, z, m, w) {
    (as.numeric(A == ap) / g(ap)) *
        g(ap) / g(as) *
        dz(z, ap, w) / r(z, ap, m, w) *
        e(as, m, w) / e(ap, m, w)
}

sample_data <- function(n) {
    w <- rbinom(n, 1, pw(1))
    a <- rbinom(n, 1, g(1))
    z <- (0.25 + 0.1*a + 0.2*w) + rnorm(n)
    m <- (0.6 + 0.1*z + 0.05*a - 0.3*w) + rnorm(n)
    y <- rbinom(n, 1, my(m, z, a, w))
    data.table::data.table(w = as.numeric(w),
                           a = as.numeric(a),
                           z = as.numeric(z),
                           m = as.numeric(m),
                           y = as.numeric(y))
}
