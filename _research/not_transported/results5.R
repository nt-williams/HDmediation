library(tidyverse)
library(glue)
library(doFuture)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

tmle <- TRUE
res <- bind_rows(read_zip_rds(glue("_research/data/sim_not_transported_dgpRealistic_hal_{tmle}_2.zip")))

direct <- c(0.096, 0.097, 0.093, 0.094, 0.099, 0.108, 0.101, 0.102, 0.092, 0.092, 0.109)
indirect <- c(0.043, 0.042, 0.049, 0.048, 0.039, 0.042, 0.043, 0.039, 0.041, 0.051, 0.039)

direct <- mean(direct)
indirect <- mean(indirect)

bind_rows(
    {
        res |>
            filter(between(direct, -1, 1)) |>
            group_by(n) |>
            summarise(estim = mean(direct), 
                      abs_bias = abs(mean(direct) - !!direct),
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, ~ between(!!direct, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "direct")
    }, {
        res |>
            filter(between(indirect, -1, 1)) |>
            group_by(n) |>
            summarise(estim = mean(indirect), 
                      abs_bias = abs(mean(indirect) - !!indirect),
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, ~ between(!!indirect, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "indirect")
    }
) |>
    select(estimand, n, estim, abs_bias, rootn_bias, covr)
