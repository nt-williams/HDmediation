suppressPackageStartupMessages(library(tidyverse))
library(kableExtra)

res <- read_csv("_research/results_sl.csv")

# local({
#     source("_research/gen_data.R", local = TRUE)
#     truth <<- truth_mediation(gen_data(3e6))
# })

ci <- function(psi, var, n) {
    psi + (c(-1, 1) * (qnorm(.975) * sqrt(var / n)))
}

res <- group_by(res, n) |> 
    nest(data = !n) |> 
    arrange(n) |> 
    mutate(bias_direct = map_dbl(data, ~abs(mean(.x$direct) - 0.368)), 
           bias_indirect = map_dbl(data, ~abs(mean(.x$indirect) - 0.074)), 
           rootnbias_direct = map_dbl(data, ~abs(mean(.x$direct) - 0.368) * sqrt(n)), 
           rootnbias_indirect = map_dbl(data, ~abs(mean(.x$indirect) - 0.074) * sqrt(n)), 
           coverage_direct = map_dbl(data, function(x) {
               map2(x$direct, x$var_direct, ci, n = n) |> 
                   map_lgl(\(bounds) between(0.368, bounds[1], bounds[2])) |> 
                   mean()
           }), 
           coverage_indirect = map_dbl(data, function(x) {
               map2(x$indirect, x$var_indirect, ci, n = n) |> 
                   map_lgl(\(bounds) between(0.074, bounds[1], bounds[2])) |> 
                   mean()
           })) |> 
    select(-data)

pivot_longer(res, !n, 
             names_to = c(".value", "effect"), 
             names_pattern = "(.+)_(.+)") |> 
    mutate(effect = str_to_title(effect)) |> 
    kbl("latex", 
        booktabs = TRUE, 
        digits = 2, 
        col.names = c("n", "Estimand", "$|Bias|$", "$\\sqrt{n} \\times |Bias|$", "95\\% CI Coverage"), 
        escape = FALSE, 
        align = c("l", "l", "c", "c", "c"),
        linesep = "")
