tmce <- function(data, A, S, W, Z, M, Y, aprime, astar, learners, family, folds = 1) {
    checkmate::assertDataFrame(data[, c(A, S, W, Z, M, Y)])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)

    tmp <- data
    # src <- tmp[tmp[[np$S]] == 0, ]

    folds <- make_folds(tmp, folds)

    # is this limited to S = 0 or predicting with S = 0!?
    gg <- g(tmp, npsem, folds, learners)
    ee <- e(tmp, npsem, folds, learners)

    hz <- h_z(data, npsem, folds, learners)
    hm <- h_m(hz, gg, ee, aprime, astar)

    bb <- b(data, npsem, family, folds, learners)

}
