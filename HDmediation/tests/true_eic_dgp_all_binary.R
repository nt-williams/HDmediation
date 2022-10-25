g <- function(a, w) {
  pscore <- .5
  return(a * pscore + (1 - a) * (1 - pscore))
}

pz <- function(z, a, w) {
  prob1 <- plogis((-log(1.3)*(rowSums(w)) / 3) + 2*a - 1)
  return(z * prob1 + (1 - z) * (1 - prob1))
}

pm <- function(m, z, a, w) {
  prob1 <- plogis(-log(1.1)*w[, 3] + 2*z - 0.9)
  return(m * prob1 + (1 - m) * (1 - prob1))
}

my_inst <- function(m, z, w) {
  plogis((-log(1.3)*(rowSums(w)) / 3) + z + m)
}

intv11 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pz(1, astar, w) * pm(1, 1, astar, w) + 
    my_inst(0, 1, w) * pz(1, astar, w) * pm(0, 1, astar, w)
}

intv10 <- function(w, aprime, astar) {
  my_inst(1, 1, w) * pm(1, 0, astar, w) * (pz(1, aprime, w) - pz(1,astar, w)) +
    my_inst(0, 1, w) * pm(0, 0, astar, w)* (pz(1, aprime, w) - pz(1, astar, w))
}

intv00 <- function(w, aprime, astar) {
  my_inst(1, 0, w) * pm(1, 0, astar, w) * pz(0, aprime, w) + 
    my_inst(0, 0, w) * pm(0, 0, astar, w) * pz(0, aprime, w) 
}

pzmw <- function(m, z, w) {
  pm(m, z, 1, w) * g(1, w) +
    pm(m, z, 0, w) * g(0, w)
}

e <- function(a, z, m, w) {
  pm(m, z, a, w) * g(a, w) / pzmw(m, z,w)
}

pmaw <- function(m, a, w) {
  pm(m, 1, a, w) * pz(1, a, w) +
    pm(m, 0, a, w) * pz(0, a, w)
}

r <- function(z, a, m, w) {
  pm(m, z, a, w) * pz(z, a, w) / pmaw(m, a, w)
}

eic <- function(data, aprime, astar) {
  A <- data$a
  Z <- data$z
  m <- data$m
  Y <- data$y
  Z <- data$z
  w <- data[, paste0("w", 1:3)]
  
  `P(a'|W)` <- g(astar, w)
  `P(a|W)` <- g(aprime, w)
  `P(a|M,1,W)` <- e(aprime, 1, m, w)
  `P(a'|M,1,W)` <- e(astar, 1, m, w)
  `P(Z=0|a',W)` <- pz(0, astar, w)
  `P(Z=0|a,W)` <- pz(0, aprime, w)
  `P(Z=1|a',W)` <- pz(1, astar, w)
  `P(Z=1|a,W)` <- pz(1, aprime, w)
  `P(Z=0|M,a',W)` <- r(0, m, astar, w)
  `P(Z=1|M,a',W)` <- r(1, m, astar, w)
  `P(a|M,0,W)` <- e(aprime, 0, m, w)
  `P(a'|M,0,W)` <- e(astar, 0, m, w)
  `E(Y|a,M,1,W)` <- my_inst(m, 1, w)
  `E(Y|a,M,0,W)` <- my_inst(m, 0, w)
  `E(Y|A,M,Z,W)` <- my_inst(m, Z, w)
  `P(Z=1|A,W)` <- pz(1, A, w)
  
  H_Y11 <- ((Z == 1 & A == aprime) / `P(a'|W)`) * (`P(a'|M,1,W)` / `P(a|M,1,W)`)
  H_Y10 <- ((Z == 1 & A == aprime) / (`P(a'|W)` * `P(Z=0|a',W)`)) * 
    ((`P(a'|M,1,W)` / `P(a|M,1,W)`) * (`P(Z=0|M,a',W)` / `P(Z=1|M,a',W)`)) * 
    (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
  H_Y00 <- ((Z == 0 & A == aprime) / `P(a|W)`) * (`P(a'|M,0,W)` / `P(a|M,0,W)`) * 
    (`P(Z=0|a,W)` / `P(Z=0|a',W)`)
  
  H_M11 <- (Z == 1 & A == astar) / `P(a'|W)`
  H_M10 <- ((Z == 0 & A == astar) / (`P(a'|W)` * `P(Z=0|a',W)`)) * (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
  H_M00 <- ((Z == 0 & A == astar) / (`P(a'|W)` * `P(Z=0|a',W)`)) * `P(Z=0|a,W)`
  
  H_Z11 <- (A == astar) / `P(a'|W)` * intv11(w, aprime, astar)
  H_Z10 <- (((A == aprime) / `P(a|W)`) - ((A == astar) / `P(a'|W)`)) * intv10(w, aprime, astar)
  H_Z00 <- -((A == aprime) / `P(a|W)`) * intv00(w, aprime, astar)
  
  H_W11 <- intv11(w, aprime, astar) * `P(Z=1|a',W)`
  H_W10 <- intv10(w, aprime, astar) * (`P(Z=1|a,W)` - `P(Z=1|a',W)`)
  H_W00 <- intv00(w, aprime, astar) * `P(Z=0|a,W)`
  
  eic_11 <- H_Y11 * (Y - `E(Y|A,M,Z,W)`) + H_Z11 * (Z - `P(Z=1|A,W)`) + 
    H_M11 * (`E(Y|a,M,1,W)` - intv11(w, aprime, astar)) + 
    H_W11
  
  eic_10 <- H_Y10 * (Y - `E(Y|A,M,Z,W)`) + H_Z10 * (Z - `P(Z=1|A,W)`) + 
    H_M10 * (`E(Y|a,M,1,W)` - intv10(w, aprime, astar)) + 
    H_W10
  
  eic_00 <- H_Y00 * (Y - `E(Y|A,M,Z,W)`) + H_Z00 * (Z - `P(Z=1|A,W)`) + 
    H_M00 * (`E(Y|a,M,0,W)` - intv00(w, aprime, astar)) + 
    H_W00
  
  eic_11 + eic_10 + eic_00
}

tmp <- simdata(1e6)
eic_11 <- eic(tmp, 1, 1)
eic_10 <- eic(tmp, 1, 0)
eic_00 <- eic(tmp, 0, 0)



