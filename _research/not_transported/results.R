suppressPackageStartupMessages(library(tidyverse))
library(glue)

source("_research/not_transported/gendata.R")

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

res <- bind_rows(read_zip(glue("_research/data/sim_not_transported.zip")))

direct <- truth()["direct"]
indirect <- truth()["indirect"]

filter(res, estimand == "direct") |> 
    group_by(n) |> 
    summarise(abs_bias = abs(mean(psi) - direct), 
              covr = mean(map2_lgl(conf_low, conf_high, ~ between(direct, .x, .y)))) |> 
    mutate(rootn_bias = abs_bias * sqrt(n))

filter(res, estimand == "indirect") |> 
    group_by(n) |> 
    summarise(abs_bias = abs(mean(psi) - indirect), 
              covr = mean(map2_lgl(conf_low, conf_high, ~ between(indirect, .x, .y)))) |> 
    mutate(rootn_bias = abs_bias * sqrt(n))
