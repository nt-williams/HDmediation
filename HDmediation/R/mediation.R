#' Interventional (in)direct effects
#'
#' Efficient, doubly-robust estimation of interventional (in)direct effects with
#' possibly high-dimensional mediators and confounders.
#'
#' @param data [\code{data.frame}]\cr
#'  A \code{data.frame} containing all necessary variables for the estimation problem.
#' @param A [\code{character(1)}]\cr
#'  A vector containing the column name of the treatment variable. \code{A} must be binary, code with 0 and 1.
#' @param W [\code{character}]\cr
#'  A vector containing the column names of baseline confounders to be
#'  included for adjustment.
#' @param Z [\code{character}]\cr
#'  A vector containing the column names of the mediator-outcome confounders.
#' @param M [\code{character}]\cr
#'  A vector containing the column names of the mediators.
#' @param Y [\code{character(1)}]\cr
#'  A vector containing the column name of the outcome variable.
#' @param S [\code{character(1)}]\cr
#'  An optional vector containing the column name for the variable indicating
#'  what population an observation belongs to. If not \code{NULL}, transported
#'  in(direct) will be estimated. \code{S} must be binary, code with 0 and 1.
#' @param family [\code{character(1)}]\cr
#'  Outcome variable type (i.e., "gaussian", "binomial").
#' @param folds [\code{integer(1)}]\cr
#'  The number of folds to be used for cross-fitting.
#'
#' @return A list of the estimates
#' @export
#'
#' @examples
#' n <- 1000
#' w0 <- rbinom(n, 1, 0.5)
#' w1 <- rbinom(n, 1, 0.5)
#'
#' probsite <- plogis(log(1.2) * w1)
#' site <- rbinom(n, 1, probsite)
#'
#' a <- rbinom(n, 1, plogis(log(2) * w1))
#' z <- rbinom(n, 1, plogis(-log(2) + log(4) * a - log(2) * w1 + log(1.4) * site))
#' m <- rbinom(n, 1, plogis(-log(2) + log(4) * a + log(2) * z - log(1.4) * w1 + log(1.4) * site))
#' y <- rbinom(n, 1, plogis(-log(5) + log(3) * a + log(1.6) * a * w1 + log(8) * z + log(4) * m - log(1.2) * w1))
#'
#' tmp <- data.frame(
#'     W0 = w0,
#'     W1 = w1,
#'     A = a,
#'     Z = z,
#'     M = m,
#'     Y = y,
#'     S = site
#' )
#'
#' mediation(tmp, "A", c("W0", "W1"), "Z", "M", "Y", "S", "binomial", 1)
mediation <- function(data, A, W, Z, M, Y, S = NULL, family = c("binomial", "gaussian"), folds = 1) {
    checkmate::assertDataFrame(data[, c(A, S, W, Z, M, Y)])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    if (!is.null(S)) {
        ans <- transported(data, A, S, W, Z, M, Y, family, folds)
    }

    if (is.null(S)) {
        ans <- not_transported(data, A, W, Z, M, Y, family, folds)
    }

    ans
}
