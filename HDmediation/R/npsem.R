Npsem <- R6::R6Class(
    "Npsem",
    public = list(
        S = NULL,
        W = NULL,
        A = NULL,
        Z = NULL,
        M = NULL,
        Y = NULL,
        cens = NULL,
        initialize = function(S = NULL, W, A, Z, M, Y, cens = NULL) {
            checkmate::assertCharacter(A)
            checkmate::assertCharacter(W)
            checkmate::assertCharacter(Y)
            checkmate::assertCharacter(Z, null.ok = TRUE)
            checkmate::assertCharacter(M)
            checkmate::assertCharacter(S, null.ok = TRUE)
            checkmate::assertCharacter(cens, null.ok = TRUE)

            self$S <- S
            self$W <- W
            self$Z <- Z
            self$M <- M
            self$A <- A
            self$Y <- Y
            self$cens <- cens

            invisible()
        },
        history = function(var = c("A", "Z", "M", "Y", "cens")) {
            switch(
                match.arg(var),
                A = private$parents_A(),
                Z = private$parents_Z(),
                M = private$parents_M(),
                Y = private$parents_Y(),
                cens = private$parents_Y()
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
