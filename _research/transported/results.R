library(tidyverse)
library(glue)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

tmle <- TRUE
dgp <- 2

res <- bind_rows(read_zip_rds(glue("_research/data/sim_transported2_{tmle}_{dgp}.zip")))

if (dgp == 1) {
    source("_research/transported/gendata3.R")
}

if (dgp == 2) {
    source("_research/transported/gendata2.R")
}

direct <- truth()["direct"]
indirect <- truth()["indirect"]

# Efficiency bounds
tmp <- gendata(5e6)
bindirect <- var(If(tmp, 1, 1) - If(tmp, 1, 0))
bdirect <- var(If(tmp, 1, 0) - If(tmp, 0, 0))

bind_rows(
    {
        res |> 
            # filter(between(direct, -1, 1)) |> 
            group_by(n) |>
            summarise(#psi = median(direct), 
                      abs_bias = abs(median(direct) - !!direct),
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, 
                                           ~ between(!!direct, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "direct")
    }, {
        res |> 
            # filter(between(indirect, -1, 1)) |> 
            group_by(n) |>
            summarise(#psi = median(indirect),
                      abs_bias = abs(median(indirect) - !!indirect),
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, 
                                           ~ between(!!indirect, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "indirect")
    }
) |>
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_transported2_{tmle}_{dgp}.rds"))
