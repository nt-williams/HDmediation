library(tidyverse)
library(glue)
library(doFuture)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

tmle <- TRUE

res <- bind_rows(read_zip_rds(glue("_research/data/sim_not_transported_dgpRealistic_{tmle}.zip")))

bind_rows(
    {
        res |>
            filter(between(direct, -1, 1)) |>
            group_by(n) |>
            summarise(abs_bias = abs(mean(direct) - !!direct),
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, ~ between(!!direct, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "direct")
    }, {
        res |>
            filter(between(indirect, -1, 1)) |>
            group_by(n) |>
            summarise(abs_bias = abs(mean(indirect) - !!indirect),
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, ~ between(!!indirect, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "indirect")
    }
) |>
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_not_transported_{tmle}_cont.rds"))
