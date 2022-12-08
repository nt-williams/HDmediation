suppressPackageStartupMessages(library(tidyverse))
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

tmle <- F
dgp <- 2
if (dgp == 1) {
    source("_research/transported/gendata.R")
    res <- bind_rows(read_zip(glue("_research/data/sim_transported_1.zip")))
}

if (dgp == 2) {
    source("_research/transported/gendata2.R")
    res <- bind_rows(read_zip(glue("_research/data/sim_transported_{tmle}_2.zip")))
}

if (dgp == 3) {
    source("_research/transported/gendata3.R")
    res <- bind_rows(read_zip(glue("_research/data/sim_transported_{tmle}_1.zip")))
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
            #filter(between(direct, -1, 1)) |> 
            group_by(n) |>
            summarise(#psi = median(direct), 
                      abs_bias = abs(median(direct) - !!direct),
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, ~ between(!!direct, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "direct")
    }, {
        res |> 
            #filter(between(indirect, -1, 1)) |> 
            group_by(n) |>
            summarise(#psi = median(indirect),
                      abs_bias = abs(median(indirect) - !!indirect),
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, ~ between(!!indirect, .x, .y)))) |>
            mutate(rootn_bias = abs_bias * sqrt(n),
                   estimand = "indirect")
    }
) |>
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_transported_{tmle}_{dgp}.rds"))
