h_m <- function(h_z, g, e, aprime, astar) {
    `g(a'|w)` <- g[, gl("g({aprime}|w)")]
    `g(a*|w)` <- g[, gl("g({astar}|w)")]
    `e(a'|m,w)` <- e[, gl("e({aprime}|m,w)")]
    `e(a*|m,w)` <- e[, gl("e({astar}|m,w)")]
    h_z[, gl("h_z({aprime})")] * `g(a'|w)` / `g(a*|w)` * `e(a*|m,w)` / `e(a'|m,w)`
}
