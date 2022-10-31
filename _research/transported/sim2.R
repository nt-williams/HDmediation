library(HDmediation)
library(glue)
suppressPackageStartupMessages(library(tidyverse))

source("_research/transported/gendata2.R")
# source("_research/SL.lightgbm.R")
# source("_research/SL.bart.R")

id <- Sys.getenv("SGE_TASK_ID")
if (id == "undefined" || id == "") id <- 1

res <- map_dfr(c(500, 1000, 5000, 1e4), function(n) {
    dat <- gendata(n)
    mediation(dat, "A", "W1",
              c("Z1", "Z2"),
              c("M1", "M2"),
              "Y", "S",
              family = "binomial",
              folds = 1)
}, .id = "n")

res <- mutate(res,
              n = case_when(
                  n == 1 ~ 500,
                  n == 2 ~ 1000,
                  n == 3 ~ 5000,
                  n == 4 ~ 1e4
              ))

saveRDS(res, glue("_research/data/sim_not_transported_{id}.rds"))
