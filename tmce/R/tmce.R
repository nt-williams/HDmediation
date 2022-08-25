tmce <- function(data, A, S, W, Z, M, Y, family, folds = 1) {
    checkmate::assertDataFrame(data[, c(A, S, W, Z, M, Y)])
    checkmate::assertNumber(folds, lower = 1, upper = nrow(data) - 1)

    npsem <- Npsem$new(A = A, S = S, W = W, Z = Z, M = M, Y = Y)
    folds <- make_folds(data, folds)

    gg <- g(data, npsem, folds)
    ee <- e(data, npsem, folds)
    cc <- see(data, npsem, folds)
    bb <- b(data, npsem, "binomial", folds)
    t <- 1 - mean(data[[npsem$S]])

    thetas <- eifs <- list()
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hz <- h_z(data, npsem, aprime, folds)
        hm <- h_m(hz, gg, ee, aprime, astar)
        uu <- u(data, npsem, bb, hm, aprime, folds)
        uubar <- ubar(data, npsem, uu, aprime, folds)
        vv <- v(data, npsem, bb, hz, aprime, folds)
        # vv <- true_v(data, aprime)
        vvbar <- vbar(data, npsem, vv, astar, folds)
        # vvbar <- matrix(true_vbar(data, aprime, astar), ncol = 1)

        # EIF calculation
        S <- data[[npsem$S]]
        A <- data[[npsem$A]]
        Y <- data[[npsem$Y]]
            
        eify <-
            ((S == 1) & (A == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (1 - cc[, gl("c({aprime},z,m,w)")]) / cc[, gl("c({aprime},z,m,w)")] *
            hm * (Y - bb[, gl("b({aprime},Z,M,W)")])

        eifz <-
            ((S == 0) & (A == aprime)) /
            (t * gg[, gl("g({aprime}|w)")]) *
            (uu[, 1] - uubar[, 1])

        eifm <-
            ((S == 0) & (A == astar)) /
            (t * gg[, gl("g({astar}|w)")]) *
            (vv[, 1] - vvbar[, 1])
        
        # eif <- eify + eifz + eifm + (1 - S) / t * (vvbar[, 1] - mean(vvbar[, 1])) + mean(vvbar[, 1])
        # 
        # theta <- mean(eif)
        # eif <- eif - (1 - S) / t * theta
        
        eifw <- (S == 0) / t * (vvbar[, 1] - mean(vvbar[S == 0, 1]))
        # # D_PW <- (data[[npsem$S]] == 0) / t * (vvbar[, 1] - mean(vvbar[, 1]))
        # 
        theta <- mean(eify + eifz + eifm + eifw) + mean(vvbar[S == 0, 1])
        # # theta <- mean(D_PY + D_PZ + D_PM + D_PW) + mean(vvbar[, 1])
        eif <- eify + eifz + eifm + eifw
        
        thetas <- c(thetas, list(theta))
        eifs <- c(eifs, list(eif))
        # components <- list(y = D_PY, 
        #      z = D_PZ, 
        #      m = D_PM, 
        #      w = (data[[npsem$S]] == 0) / t * (vvbar[, 1] - mean(vvbar[(data[[npsem$S]] == 0), 1])))
        # 
        # eifs <- c(eifs, list(components))
    }

    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")
    
    ans <- list(indirect = thetas$`11` - thetas$`10`, 
                direct = thetas$`10` - thetas$`00`)
    
    ans$var_indirect <- var(eifs$`11` - eifs$`10`)
    ans$var_direct <- var(eifs$`10` - eifs$`00`)

    # ans <- list(
    #     indirect = mean(eifs$`11` - eifs$`10`),
    #     direct = mean(eifs$`10` - eifs$`00`)
    # )
    # 
    # ans$var_indirect <- var(eifs$`11` - eifs$`10` - (data[[npsem$S]] == 0) / t * ans$indirect)
    # ans$var_direct <- var(eifs$`10` - eifs$`00` - (data[[npsem$S]] == 0) / t * ans$direct)
    
    ans
}
