suppressPackageStartupMessages({
    library(tidyverse)
    library(kableExtra)
    
})

transport <- TRUE
if (transport) {
    dgp1t <- readRDS("_research/data/summary_transported_TRUE_3.rds")
    dgp2t <- readRDS("_research/data/summary_transported_TRUE_2.rds")
    dgp3t <- readRDS("_research/data/summary_transported_TRUE_cont.rds")
    dgp1f <- readRDS("_research/data/summary_transported_FALSE_3.rds")
    dgp2f <- readRDS("_research/data/summary_transported_FALSE_2.rds")
    dgp3f <- readRDS("_research/data/summary_transported_FALSE_cont.rds")
} else {
    dgp1t <- readRDS("_research/data/summary_not_transported_TRUE_3.rds")
    dgp2t <- readRDS("_research/data/summary_not_transported_TRUE_2.rds")
    dgp3t <- readRDS("_research/data/summary_not_transported_TRUE_cont.rds")
    dgp1f <- readRDS("_research/data/summary_not_transported_FALSE_3.rds")
    dgp2f <- readRDS("_research/data/summary_not_transported_FALSE_2.rds")
    dgp3f <- readRDS("_research/data/summary_not_transported_FALSE_cont.rds")
}

filter(dgp2f, estimand == "direct") |>
    bind_rows(filter(dgp2t, estimand == "direct"), 
              filter(dgp1f, estimand == "direct"), 
              filter(dgp1t, estimand == "direct"), 
              filter(dgp3f, estimand == "direct"), 
              filter(dgp3t, estimand == "direct")) |>
    select(-estimand) |>
    bind_cols({
        filter(dgp2f, estimand == "indirect") |>
            bind_rows(filter(dgp2t, estimand == "indirect"), 
                      filter(dgp1f, estimand == "indirect"), 
                      filter(dgp1t, estimand == "indirect"), 
                      filter(dgp3f, estimand == "indirect"), 
                      filter(dgp3t, estimand == "indirect")) |>
            select(-estimand, -n) 
    }) |> 
    kbl("latex", booktabs = TRUE,
        digits = 2,
        linesep = "",
        col.names = c("$n$", "$\\text{Bias}$", "$\\sqrt{n} \\times \\text{Bias}$", "95\\% CI Covr.", 
                      "$\\text{Bias}$", "$\\sqrt{n} \\times \\text{Bias}$", "95\\% CI Covr."),
        escape = FALSE,
        align = "rcccccc")
