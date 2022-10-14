Npsem <- R6::R6Class(
    "Npsem",
    public = list(
        S = NULL,
        W = NULL,
        A = NULL,
        Z = NULL,
        M = NULL,
        Y = NULL,
        initialize = function(S = NULL, W, A, Z, M, Y) {
            checkmate::assertCharacter(A)
            checkmate::assertCharacter(W)
            checkmate::assertCharacter(Y)
            checkmate::assertCharacter(Z)
            checkmate::assertCharacter(M)
            checkmate::assertCharacter(S, null.ok = TRUE)

            self$S <- S
            self$W <- W
            self$Z <- Z
            self$M <- M
            self$A <- A
            self$Y <- Y

            invisible()
        },
        history = function(var = c("A", "Z", "M", "Y")) {
            switch(
                match.arg(var),
                A = private$parents_A(),
                Z = private$parents_Z(),
                M = private$parents_M(),
                Y = private$parents_Y()
            )
        }
    ),
    private = list(
        parents_A = function() {
            c(self$S, self$W)
        },
        parents_Z = function() {
            c(private$parents_A(), self$A)
        },
        parents_M = function() {
            c(private$parents_Z(), self$Z)
        },
        parents_Y = function() {
            c(private$parents_M(), self$M)
        }
    )
)
