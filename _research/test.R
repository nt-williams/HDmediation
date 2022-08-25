source("_research/sim/gen_data.R")
# source("_research/sim/truth.R")

library(devtools)

load_all("tmce")

n <- 1e4
dat <- gen_data(n)

A <- "A"
S <- "S"
W <- c("W0", "W1")
Z <- "Z"
M <- "M"
Y <- "Y"

res <- tmce(dat, A, S, W, Z, M, Y, "binomial", folds = 10)

se <- sqrt(var(res$eif$`11` - res$eif$`10`) / n)
(res$theta$`11` - res$theta$`10`) + c(-1, 1) * se * qnorm(0.975)

se <- sqrt(var(res$eif$`10` - res$eif$`00`) / n)
(res$theta$`10` - res$theta$`00`) + c(-1, 1) * se * qnorm(0.975)

res$theta$`11` - res$theta$`10`
res$theta$`10` - res$theta$`00`

library(purrr)

tots <- 
    map(res, \(x) map(x, mean)) |> 
    map(\(x) reduce(x, sum))

# tots$`11` - tots$`10`
# tots$`10` - tots$`00`

map(res, \(x) map(x, mean))

local({
    source("_research/sim/truth.R", local = TRUE)
    truth <<- truth_mediation(gendata(3e6))
})

local({
    source("_research/sim/truth.R", local = TRUE)
    comp <<- list(
        `11` = compute_eif2(gen_data(3e6), 1, 1),
        `10` = compute_eif2(gen_data(3e6), 1, 0),
        `00` = compute_eif2(gen_data(3e6), 0, 0)
    )
})

abs(mean(comp$`11`$y) - mean(res$`11`$y)) * sqrt(n)
abs(mean(comp$`11`$z) - mean(res$`11`$z)) * sqrt(n)
abs(mean(comp$`11`$m) - mean(res$`11`$m)) * sqrt(n)
abs(mean(comp$`11`$w) - mean(res$`11`$w)) * sqrt(n)

# abs(mean(comp$`10`$y) - mean(res$`10`$y)) * sqrt(n)
# abs(mean(comp$`10`$z) - mean(res$`10`$z)) * sqrt(n)
# abs(mean(comp$`10`$m) - mean(res$`10`$m)) * sqrt(n)
# abs(mean(comp$`10`$w) - mean(res$`10`$w)) * sqrt(n)
# 
# abs(mean(comp$`00`$y) - mean(res$`00`$y)) * sqrt(n)
# abs(mean(comp$`00`$z) - mean(res$`00`$z)) * sqrt(n)
# abs(mean(comp$`00`$m) - mean(res$`00`$m)) * sqrt(n)
# abs(mean(comp$`00`$w) - mean(res$`00`$w)) * sqrt(n)

local({
    source("_research/sim/truth.R", local = TRUE)
    comp <<- list(
        `11` = compute_eif2(dat, 1, 1), 
        `10` = compute_eif2(dat, 1, 0), 
        `00` = compute_eif2(dat, 0, 0)
    )
})

plot(comp$`11`$y, res$`11`$y)
plot(comp$`11`$z, res$`11`$z)
plot(comp$`11`$m, res$`11`$m)
plot(comp$`11`$w, res$`11`$w)
