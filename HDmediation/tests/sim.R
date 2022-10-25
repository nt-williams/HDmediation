source("HDmediation/tests/no_transport.R")

library(devtools)
load_all("HDmediation")

id <- Sys.getenv("SGE_TASK_ID")

truth()

n <- 1e4
dat <- simdata(n)

psi <- mediation(dat, "A", paste0("W_", 1:3), "Z", "M", "Y", S = NULL, family = "binomial", folds = 1) 

res <- data.frame(n = n,
                  direct = psi$direct,
                  var_direct = psi$var_direct,
                  indirect = psi$indirect,
                  var_indirect = psi$var_indirect)

write.csv(res, paste0("HDmediation/tests/data/sim_", id, ".csv"))
