h_m <- function(h_z, g, e, aprime, astar) {
    `g(a'|w)` <- g[, gl("g({aprime}|w)")]
    `g(a*|w)` <- g[, gl("g({astar}|w)")]
    `e(a'|m,w)` <- e[, gl("e({aprime}|m,w)")]
    `e(a*|m,w)` <- e[, gl("e({astar}|m,w)")]
    h_z[, 1] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}

# true_hm <- function(data, aprime, astar) {
#     pz <- function(z, a, w, s) {
#         prob1 <- plogis(-log(2) + (log(4) * a) - log(2) * w[, "W1"] + log(1.4) * s + log(1.43) * s * a)
#         z * prob1 + (1 - z) * (1 - prob1)
#     }
#     
#     pm <- function(m, z, a, w, s) {
#         prob1 <- plogis(-log(2) + log(4) * z - log(1.4) * w[, "W1"] + log(1.4) * s)
#         m * prob1 + (1 - m) * (1 - prob1)
#     }
#     
#     pmaw <- function(m, a, w, s) {
#         pm(m, 1, a, w, s) * pz(1, a, w, s) + pm(m, 0, a, w, s ) * pz(0, a, w, s)
#     }
#     
#     z <- data$Z
#     m <- data$M
#     w <- data[, c("W0", "W1")]
#     
#     pmaw(m, astar, w, 0) / pm(m, z, aprime, w, 0) 
# }
