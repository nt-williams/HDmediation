library(tidyverse)
library(glue)

read_zip <- function(tar) {
    files <- unzip(tar, list = TRUE)$Name
    p <- progressr::progressor(along = 1:length(files))
    purrr::map(files, function(file) {
        p()
        con <- gzcon(unz(tar, file))
        x <- readRDS(con)
        close(con)
        x
    })
}

tmle <- FALSE
res <- bind_rows(read_zip(glue("_research/data/sim_transported_{tmle}_cont.zip")))

indirect <- 0.0146
direct <- 0.0266

summary(filter(res, n == 1e4)$indirect)
summary(filter(res, n == 1e4)$direct)

bind_rows(
    {
        res |>
            group_by(n) |>
            summarise(abs_bias = abs(mean(direct) - !!direct),
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, ~ between(!!direct, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "direct")
    }, {
        res |>
            group_by(n) |>
            summarise(abs_bias = abs(mean(indirect) - !!indirect),
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, ~ between(!!indirect, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "indirect")
    }
) |>
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_transported_{tmle}_cont.rds"))
