library(tidyverse)
library(glue)
library(doFuture)

devtools::source_gist("https://gist.github.com/nt-williams/3afb56f503c7f98077722baf9c7eb644")

source("_research/transported/gendata4.R")

tmp <- sample_data(1e6)

plan(multisession)

ipw00 <- foreach(a = tmp$a, s = tmp$S, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
                     ipw(a, 0, 0, s, z, m, w)
                 }

ipw10 <- foreach(a = tmp$a, s = tmp$S, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
                     ipw(a, 1, 0, s, z, m, w)
                 }

ipw11 <- foreach(a = tmp$a, s = tmp$S, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
                     ipw(a, 1, 1, s, z, m, w)
                 }

plan(sequential)

y11 <- mean(ipw11 * tmp$y)
y10 <- mean(ipw10 * tmp$y)
y00 <- mean(ipw00 * tmp$y)
indirect <- y11 - y10
direct <- y10 - y00

tmle <- FALSE
res <- map_dfr(1:2, function(i) bind_rows(read_zip_rds(glue("_research/data/sim_transported_{tmle}_cont_{i}.zip"))))

bind_rows({
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
}) |>
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_transported_{tmle}_cont.rds"))
