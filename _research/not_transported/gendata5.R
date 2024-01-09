pw <- function(w) {
    prob1 <- 0.4
    w * prob1 + (1 - w) * (1 - prob1)
}

g <- function(a) {
    pscore <- .5
    a * pscore + (1 - a) * (1 - pscore)
}

# binary Z
pz1 <- function(z1, a, w) {
    prob1 <- plogis(0.5 + a - 0.75*w)
    z1 * prob1 + (1 - z1) * (1 - prob1)
}

# ordinal logit model for 2
pz2 <- function(a, w) {
    beta <- matrix(c(0, 0.2, 0.1, 
                     -0.75, 0.2, 0.1, 
                     -2, 0.2, 0.1, 
                     -3, 0.2, 0.1), 
                   ncol = 3, 
                   byrow = TRUE)
    
    x <- matrix(c(rep(1, length(a)), a, w), ncol = 3)
    
    probs <- rbind(matrix(rep(1, length(a)), nrow = 1), 
                   plogis(beta %*% t(x)), 
                   matrix(rep(0, length(a)), nrow = 1))
    
    probs <- diff(probs[nrow(probs):1, , drop = FALSE])
    t(probs[nrow(probs):1, , drop = FALSE])
}

dz2 <- function(z2, a, w) {
    probs <- pz2(a, w)
    purrr::map2_dbl(seq_along(a), z2, function(i, j) probs[i, j])
}

rz2 <- function(a, w) {
    LaplacesDemon::rcat(length(a), pz2(a, w))
}

# ordinal logit model for Z3
pz3 <- function(a, w) {
    beta <- matrix(c(-0.25, 0.2, 0.1, 
                     -0.75, 0.2, 0.1, 
                     -2, 0.2, 0.1, 
                     -2.75, 0.2, 0.1), 
                   ncol = 3, 
                   byrow = TRUE)
    
    x <- matrix(c(rep(1, length(a)), a, w), ncol = 3)
    
    probs <- rbind(matrix(rep(1, length(a)), nrow = 1), 
                   plogis(beta %*% t(x)), 
                   matrix(rep(0, length(a)), nrow = 1))
    
    probs <- diff(probs[nrow(probs):1, , drop = FALSE])
    t(probs[nrow(probs):1, , drop = FALSE])
}

dz3 <- function(z3, a, w) {
    probs <- pz3(a, w)
    purrr::map2_dbl(seq_along(a), z3, function(i, j) probs[i, j])
}

rz3 <- function(a, w) {
    LaplacesDemon::rcat(length(a), pz3(a, w))
}

# continuous Z, bounded [0, 1] (beta distributed)
dz4 <- function(z4, a, w) {
    shape1 <- 1.9 + 0.075*a + 0.02*w
    dbeta(z4, shape1, 2)
}

rz4 <- function(a, w) {
    shape1 <- 1.9 + 0.075*a + 0.02*w
    rbeta(length(a), shape1, 2)
}

# Joint distribution of all Z
dz <- function(z1, z2, z3, z4, a, w) {
    pz1(z1, a, w)*dz2(z2, a, w)*dz3(z3, a, w)*dz4(z4, a, w)
}

# continuous M, bounded [0, 1] (beta distributed)
# dm1 <- function(m1, z1, z2, z3, z4, a, w) {
#     shape1 <- 2 + 0.01*z1 + 0.03*z2 - 0.02*z3 + 0.5*a + 0.03*w
#     dbeta(m1, shape1, 2)
# }

dm1 <- function(m1, z1, z2, z3, z4, a, w) {
    mu <- 1 + 0.01*z1 + 0.03*z2 - 0.02*z3 + 0.5*a + 0.03*w
    dnorm(m1, mu)
}

# rm1 <- function(z1, z2, z3, z4, a, w) {
#     shape1 <- 2 + 0.01*z1 + 0.03*z2 - 0.02*z3 + 0.5*a + 0.03*w
#     rbeta(length(a), shape1, 2)
# }

rm1 <- function(z1, z2, z3, z4, a, w) {
    mu <- 1 + 0.01*z1 + 0.03*z2 - 0.02*z3 + 0.5*a + 0.03*w
    rnorm(length(a), mu)
}

# continuous M, bounded [0, 1] (beta distributed)
# dm2 <- function(m2, z1, z2, z3, z4, a, w) {
#     shape1 <- 1.5 + 0.04*z2 + 0.03*z4 + .75*a + 0.025*w
#     dbeta(m2, shape1, 2)
# }

dm2 <- function(m2, z1, z2, z3, z4, a, w) {
    mu <- 1 + 0.04*z2 + 0.03*z4 + .75*a + 0.025*w
    dnorm(m2, mu)
}

# rm2 <- function(z1, z2, z3, z4, a, w) {
#     shape1 <- 1.5 + 0.04*z2 + 0.03*z4 + .75*a + 0.025*w
#     rbeta(length(a), shape1, 2)
# }

rm2 <- function(z1, z2, z3, z4, a, w) {
    mu <- 1 + 0.04*z2 + 0.03*z4 + .75*a + 0.025*w
    rnorm(length(a), mu)
}

# Joint distribution of all M
dm <- function(m1, m2, z1, z2, z3, z4, a, w) {
    dm1(m1, z1, z2, z3, z4, a, w)*dm2(m2, z1, z2, z3, z4, a, w)
}

dm_noz1 <- function(m1, m2, z2, z3, z4, a, w) {
    dm(m1, m2, z1 = 0, z2, z3, z4, a, w)*pz1(0, a, w) + 
        dm(m1, m2, z1 = 1, z2, z3, z4, a, w)*pz1(1, a, w)
}

dm_noz1z2 <- function(m1, m2, z3, z4, a, w) {
    dm_noz1(m1, m2, z2 = 1, z3, z4, a, w)*dz2(1, a, w) + 
        dm_noz1(m1, m2, z2 = 2, z3, z4, a, w)*dz2(2, a, w) + 
        dm_noz1(m1, m2, z2 = 3, z3, z4, a, w)*dz2(3, a, w) + 
        dm_noz1(m1, m2, z2 = 4, z3, z4, a, w)*dz2(4, a, w) + 
        dm_noz1(m1, m2, z2 = 5, z3, z4, a, w)*dz2(5, a, w)
}

dm_noz1z2z3 <- function(m1, m2, z4, a, w) {
    dm_noz1z2(m1, m2, z3 = 1, z4, a, w)*dz3(1, a, w) +
        dm_noz1z2(m1, m2, z3 = 2, z4, a, w)*dz3(2, a, w) + 
        dm_noz1z2(m1, m2, z3 = 3, z4, a, w)*dz3(3, a, w) + 
        dm_noz1z2(m1, m2, z3 = 4, z4, a, w)*dz3(4, a, w) + 
        dm_noz1z2(m1, m2, z3 = 5, z4, a, w)*dz3(5, a, w)
}

dmaw <- function(m1, m2, a, w) {
    integrate(function(z) dm_noz1z2z3(m1, m2, z4 = z, a, w), 0, 1)$value
}

my <- function(m1, m2, z1, z2, z3, z4, a, w) {
    plogis(-6 + log(5)*a + log(8)*z2 + log(6)*m1 + log(4)*m2 - log(1.2)*w + .25*z1 + log(2)*z3 + log(3)*w + z4)
}

sample_data <- function(n) {
    w <- rbinom(n, 1, pw(1))
    a <- rbinom(n, 1, g(1))
    z1 <- rbinom(n, 1, pz1(1, a, w))
    z2 <- rz2(a, w)
    z3 <- rz3(a, w)
    z4 <- rz4(a, w)
    m1 <- rm1(z1, z2, z3, z4, a, w)
    m2 <- rm2(z1, z2, z3, z4, a, w)
    y <- rbinom(n, 1, my(m1, m2, z1, z2, z3, z4, a, w))
    data.table::data.table(w = as.numeric(w),
                           a = as.numeric(a),
                           z1 = as.numeric(z1),
                           z2 = as.numeric(z2),
                           z3 = as.numeric(z3),
                           z4 = as.numeric(z4),
                           m1 = as.numeric(m1),
                           m2 = as.numeric(m2),
                           y = as.numeric(y))
}
