suppressPackageStartupMessages(library(tidyverse))
library(glue)

read_zip <- function(tar) {
    files <- unzip(tar, list = TRUE)$Name
    p <- progressr::progressor(along = 1:length(files))
    purrr::map(files, function(file) {
        p()
        con <- unz(tar, file)
        read.csv(con)
    })
}

n <- 10000
res <- read_zip(glue("_research/data/res_sl_{n}.zip"))
res <- bind_rows(res)

write_csv(res, "_research/results_sl.csv", append = file.exists("_research/results_sl.csv"))
