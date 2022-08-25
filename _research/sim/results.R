suppressPackageStartupMessages(library(tidyverse))

read_zip <- function(tar) {
    files <- unzip(tar, list = TRUE)$Name
    p <- progressr::progressor(along = 1:length(files))
    purrr::map(files, function(file) {
        p()
        con <- unz(tar, file)
        read.csv(con)
    })
}

local({
    source("_research/sim/truth.R", local = TRUE)
    truth <<- truth_mediation(gendata(3e6))
})

indirect <- subset(truth, parameter == "indirect", select = "truth", drop = TRUE)
direct <- subset(truth, parameter == "direct", select = "truth", drop = TRUE)

res <- read_zip("_research/sim/data/res-5000.zip")
res <- bind_rows(res)

# convergence
abs(mean(res$direct) - direct) 
#* sqrt(5000)
abs(mean(res$indirect) - indirect) 
#* sqrt(5000)

mean(res$var_direct)
mean(res$var_indirect)

# coverage
ci <- function(psi, var, n) {
    psi + (c(-1, 1) * (qnorm(.975) * sqrt(var / n)))
}

map2(res$indirect, res$var_indirect, ci, n = 5000) |> 
    map_lgl(\(bounds) between(indirect, bounds[1], bounds[2])) |> 
    mean()

map2(res$direct, res$var_direct, ci, n = 5000) |> 
    map_lgl(\(bounds) between(direct, bounds[1], bounds[2])) |> 
    mean()
