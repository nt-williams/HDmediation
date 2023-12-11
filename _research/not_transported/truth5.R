dat <- sample_data(1e5)

aprime <- 1
astar <- 1

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

psi_11 <- mean((as.numeric(dat$a == aprime) / g(aprime)) * (prob_maw / prob_mazw)*dat$y)

aprime <- 1
astar <- 0

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

psi_10 <- mean((as.numeric(dat$a == aprime) / g(aprime)) * (prob_maw / prob_mazw)*dat$y)

aprime <- 0
astar <- 0

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

psi_00 <- mean((as.numeric(dat$a == aprime) / g(aprime)) * (prob_maw / prob_mazw)*dat$y)

psi_11 - psi_10
psi_10 - psi_00
