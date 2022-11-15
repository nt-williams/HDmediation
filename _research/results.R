suppressPackageStartupMessages({
    library(tidyverse)
    library(kableExtra)
})

ntrnsprt <- readRDS("_research/data/summary_not_transported.rds")
trnsprt <- readRDS("_research/data/summary_transported.rds")

filter(ntrnsprt, estimand == "direct") |>
    bind_rows(filter(trnsprt, estimand == "direct")) |>
    select(-estimand) |>
    kbl("latex", booktabs = TRUE,
        digits = 3,
        linesep = "",
        col.names = c("$n$", "$\\text{Bias}$", "$\\sqrt{n} \\times \\text{Bias}$", "95\\% CI Covr."),
        escape = FALSE,
        align = "lccc")

filter(ntrnsprt, estimand == "indirect") |>
    bind_rows(filter(trnsprt, estimand == "indirect")) |>
    select(-estimand) |>
    kbl("latex", booktabs = TRUE,
        digits = 3,
        linesep = "",
        col.names = c("$n$", "$\\text{Bias}$", "$\\sqrt{n} \\times \\text{Bias}$", "95\\% CI Covr."),
        escape = FALSE,
        align = "lccc")
