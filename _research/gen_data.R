gen_data <- function(n = 1e4) {
    w0 <- rbinom(n, 1, .5)
    w1 <- rbinom(n, 1, .4 + (.2 * w0))
    
    probsel <- plogis(-1 + log(4) * w1 + log(4) * w0)
    psel <- rbinom(n, 1, probsel)
    
    probsite <- plogis(log(1.2) * w1 + log(1.2) * w0 + log(1.2) * w0 * w1)
    site <- rbinom(n, 1, probsite)
    
    a <- rbinom(n, 1, .5)
    z <- rbinom(n, 1, plogis(-log(2) + log(4) * a - log(2) * w1 + log(1.4) * site + log(1.43) * site * a))
    m <- rbinom(n, 1, plogis(-log(2) + log(4) * z - log(1.4) * w1 + log(1.4) * site))
    y <- rbinom(n, 1, plogis(-log(5) + log(8) * z + log(4) * m - log(1.2) * w1 + log(1.2) * w1 * z))
    
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
