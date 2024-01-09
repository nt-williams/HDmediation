dat <- sample_data(1e5)

aprime <- 1
astar <- 1

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

ipwy <- (as.numeric(dat$a == aprime) / 0.5)
h <- prob_maw / prob_mazw

y11 <- mean((ipwy*h / mean(ipwy*h))*dat$y)
y11 <- mean(ipwy*h*dat$y)

aprime <- 1
astar <- 0

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

ipwy <- (as.numeric(dat$a == aprime) / g(aprime))
h <- prob_maw / prob_mazw

y10 <- mean((ipwy*h / mean(ipwy*h))*dat$y)
y10 <- mean(ipwy*h*dat$y)

aprime <- 0
astar <- 0

prob_mazw <- dm(dat$m1, dat$m2, dat$z1, dat$z2, dat$z3, dat$z4, aprime, dat$w)
prob_maw <- vector("numeric", length = nrow(dat))
for (i in 1:nrow(dat)) {
    prob_maw[i] <- dmaw(dat$m1[i], dat$m2[i], astar, dat$w[i])
}

ipwy <- (as.numeric(dat$a == aprime) / g(aprime))
h <- prob_maw / prob_mazw

y00 <- mean((ipwy*h / mean(ipwy*h))*dat$y)
y00 <- mean(ipwy*h*dat$y)

indirect <- y11 - y10
direct <- y10 - y00
