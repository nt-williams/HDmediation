source("_research/transported/gendata.R")

n <- 1e5
set.seed(7856)
dat <- gendata(n)

converged <- function(x, truth) {
    plot(truth, x)
    abline(0, 1)
}

w <- dat[, "W1", drop = F]

# g
converged(inspect$g[, 1], rep(0.5, n))

# e
converged(inspect$e[, 1], e(0, dat$M1, dat$M2, w))
converged(inspect$e[, 2], e(1, dat$M1, dat$M2, w))

# c
converged(inspect$c[, 1], psazmw(0, dat$Z1, dat$Z2, dat$M1, dat$M2, w))
converged(inspect$c[, 2], psazmw(1, dat$Z1, dat$Z2, dat$M1, dat$M2, w))

# b
converged(inspect$b[, 1], my(dat$M1, dat$M2, dat$Z1, dat$Z2, 0, w))
converged(inspect$b[, 2], my(dat$M1, dat$M2, dat$Z1, dat$Z2, 1, w))

# hz
converged(inspect$hz[, 1], pz(dat$Z1, dat$Z2, 0, w, 0) /
              r(dat$Z1, dat$Z2, 0, dat$M1, dat$M2, w))

converged(inspect$hz[, 2], pz(dat$Z1, dat$Z2, 1, w, 0) /
              r(dat$Z1, dat$Z2, 1, dat$M1, dat$M2, w))

# hm, u, ubar, v, vbar
aprime <- astar <- 1
converged(inspect$hm_11,
          pmaw(dat$M1, dat$M2, astar, w, 0) /
              pm(dat$M1, dat$M2, dat$Z1, dat$Z2, aprime, w, 0))

converged(inspect$u_11, u(dat$Z1, dat$Z2, w, aprime, astar))

converged(inspect$ubar_11, intu(w, aprime, astar))

converged(inspect$v_11, intv(dat$M1, dat$M2, w, aprime))

astar <- 0
converged(inspect$hm_10,
          pmaw(dat$M1, dat$M2, astar, w, 0) /
              pm(dat$M1, dat$M2, dat$Z1, dat$Z2, aprime, w, 0))

converged(inspect$u_10, u(dat$Z1, dat$Z2, w, aprime, astar))
converged(inspect$ubar_10, intu(w, aprime, astar))
converged(inspect$v_10, intv(dat$M1, dat$M2, w, aprime))

aprime <- 0
converged(inspect$hm_00,
          pmaw(dat$M1, dat$M2, astar, w, 0) /
              pm(dat$M1, dat$M2, dat$Z1, dat$Z2, aprime, w, 0))

converged(inspect$u_00, u(dat$Z1, dat$Z2, w, aprime, astar))
converged(inspect$ubar_00, intu(w, aprime, astar))
converged(inspect$v_00, intv(dat$M1, dat$M2, w, aprime))
