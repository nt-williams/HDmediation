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

dgp <- 1
if (dgp == 1) {
    source("_research/not_transported/gendata.R")
    res <- bind_rows(read_zip(glue("_research/data/sim_not_transported_sat_1.zip")))
} else {
    source("_research/not_transported/gendata2.R")
    res <- bind_rows(read_zip(glue("_research/data/sim_not_transported_sat_2.zip")))
}

direct <- truth()["direct"]
indirect <- truth()["indirect"]

bind_rows(
    {
        filter(res, between(direct, 0, 1)) |> 
            group_by(n) |> 
            summarise(abs_bias = abs(mean(direct) - !!direct), 
                      covr = mean(map2_lgl(ci_direct_low, ci_direct_high, ~ between(!!direct, .x, .y)))) |> 
            mutate(rootn_bias = abs_bias * sqrt(n), 
                   estimand = "direct")
    }, {
        filter(res, between(indirect, 0, 1)) |> 
            group_by(n) |> 
            summarise(abs_bias = abs(mean(indirect) - !!indirect), 
                      covr = mean(map2_lgl(ci_indirect_low, ci_indirect_high, ~ between(!!indirect, .x, .y)))) |> 
            mutate(rootn_bias = abs_bias * sqrt(n), 
                   estimand = "indirect")
    }
) |> 
    select(estimand, n, abs_bias, rootn_bias, covr) |>
    saveRDS(glue("_research/data/summary_not_transported_sat_{dgp}.rds"))
