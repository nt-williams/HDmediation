not_transported <- function(data, A, W, Z, M, Y, cens,
                            family, folds = 1, partial_tmle, bounds,
                            learners_g = "glm",
                            learners_e = "glm",
                            learners_b = "glm",
                            learners_hz = "glm",
                            learners_u = "glm",
                            learners_ubar = "glm",
                            learners_v = "glm",
                            learners_vbar = "glm",
                            learners_cens = "glm") {
    npsem <- Npsem$new(A = A, W = W, Z = Z, M = M, Y = Y, cens = cens)
    folds <- make_folds(data, folds)

    bounds <- scale_y(data[[npsem$Y]], family, bounds)
    data[[npsem$Y]] <- Y <- bounds$y
    Y <- ifelse(is.na(Y), -999, Y)
    A <- data[[npsem$A]]

    gg <- g(data, npsem, folds, learners_g)
    ee <- e(data, npsem, folds, learners_e)
    bb <- b(data, npsem, family, folds, learners_b)
    qq <- mY(data, npsem, family, folds, learners_b)
    hz <- h_z(data, npsem, folds, learners_hz)

    if (!is.null(cens)) {
        prob_obs <- pObs(data, npsem, folds, learners_cens)
    }

    thetas <- ipws <- eifs <- list()
    vvbar <- matrix(nrow = nrow(data), ncol = 3)
    colnames(vvbar) <- c("00", "10", "11")
    for (param in list(c(1, 1), c(1, 0), c(0, 0))) {
        aprime <- param[1]
        astar <- param[2]

        hm <- h_m(hz, gg, ee, aprime, astar)

        if (!is.null(cens)) {
            obs <- data[[npsem$cens]]
            ipcw_ap <- obs / prob_obs[, gl("P(delta=1|A={aprime},Z,M,W)")]
            ipcw_as <- obs / prob_obs[, gl("P(delta=1|A={astar},Z,M,W)")]
        } else {
            ipcw_as <- ipcw_ap <- 1
        }

        ipwy <- ((A == aprime) / gg[, gl("g({aprime}|w)")])*ipcw_ap
        if (partial_tmle) {
            fit <- glm(Y ~ 1, offset = qlogis(bb[, gl("b({aprime},Z,M,W)")]), family = "binomial",
                       subset = A == aprime, weights = ipwy * hm / mean(ipwy * hm))
            bb[, gl("b({aprime},Z,M,W)")] <- plogis(coef(fit) + qlogis(bb[, gl("b({aprime},Z,M,W)")]))
        }

        uu <- u(data, npsem, bb, hm, aprime, folds, learners_u)
        uubar <- ubar(data, npsem, uu, aprime, folds, learners_ubar)
        vv <- v(data, npsem, bb, hz, aprime, folds, learners_v)
        vvbar[, paste(param, collapse = "")] <- vbar(data, npsem, vv, astar, folds, learners_vbar)

        # EIF calculation
        eify <- ipwy * hm / mean(ipwy * hm) * (Y - bb[, gl("b({aprime},Z,M,W)")])

        ipwz <- ((A == aprime) / gg[, gl("g({aprime}|w)")])*ipcw_ap
        eifz <- ipwz / mean(ipwz) * (uu[, 1] - uubar[, 1])

        ipwm <- ((A == astar) / gg[, gl("g({astar}|w)")])*ipcw_as
        eifm <- ipwm / mean(ipwm) * (vv[, 1] - vvbar[, paste(param, collapse = "")])

        eif <- rescale_y(eify + eifz + eifm + vvbar[, paste(param, collapse = "")], bounds$bounds)
        theta <- mean(eif)

        thetas <- c(thetas, list(theta))
        ipws <- c(ipws, list(mean(ipwy * hm / mean(ipwy * hm) * Y)))
        eifs <- c(eifs, list(eif))
    }
    
    if (!is.null(cens)) {
        obs <- data[[npsem$cens]]
        ipcw_a1 <- obs / prob_obs[, "P(delta=1|A=1,Z,M,W)"]
        ipcw_a0 <- obs / prob_obs[, "P(delta=1|A=0,Z,M,W)"]
    } else {
        ipcw_a1 <- ipcw_a0 <- 1
    }
    
    Y <- ifelse(is.na(data[[npsem$Y]]), -999, data[[npsem$Y]])
    mYa <- data[[npsem$A]]*qq[, "Q(1,W)"] + (1 - data[[npsem$A]])*qq[, "Q(0,W)"]
    H_a <- ((data[[npsem$A]] / gg[, "g(1|w)"])*ipcw_a1) - 
        (((1 - data[[npsem$A]]) / gg[, "g(0|w)"])*ipcw_a0)
    eif_ate <- (Y - mYa)*H_a + qq[, "Q(1,W)"] - qq[, "Q(0,W)"]

    names(eifs) <- c("11", "10", "00")
    names(thetas) <- c("11", "10", "00")
    names(ipws) <- c("11", "10", "00")

    ans <- data.frame(ate = mean(eif_ate), 
                      total = thetas$`11` - thetas$`00`,
                      indirect = thetas$`11` - thetas$`10`,
                      direct = thetas$`10` - thetas$`00`,
                      gcomp_indirect = mean(vvbar[, "11"] - vvbar[, "10"]),
                      gcomp_direct = mean(vvbar[, "10"] - vvbar[, "00"]), 
                      ipw_indirect = ipws$`11` - ipws$`10`, 
                      ipw_direct = ipws$`10` - ipws$`00`)
    
    eif_total <- (eifs$`11` - eifs$`00`) - ans$total
    eif_ate <- eif_ate - ans$ate

    ans$var_ate <- var(eif_ate)
    ans$var_total <- var(eif_total)
    ans$var_indirect <- var(eifs$`11` - eifs$`10`)
    ans$var_direct <- var(eifs$`10` - eifs$`00`)
    
    p.value <- pnorm(abs(ans$ate - ans$total) / (sqrt(var(eif_ate - eif_total) / nrow(data))), lower.tail = FALSE) * 2
    cat("Test for difference between total effect and ATE, p-value: ", round(p.value, 4), "\n")

    ci_ate <- ans$ate + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_ate / nrow(data))
    ci_total <- ans$total + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_total / nrow(data))
    ci_indirect <- ans$indirect + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_indirect / nrow(data))
    ci_direct <- ans$direct + c(-1, 1) * qnorm(0.975) * sqrt(ans$var_direct / nrow(data))

    ans$ci_ate_low <- ci_ate[1]
    ans$ci_ate_high <- ci_ate[2]
    ans$ci_total_low <- ci_total[1]
    ans$ci_total_high <- ci_total[2]
    ans$ci_indirect_low <- ci_indirect[1]
    ans$ci_indirect_high <- ci_indirect[2]
    ans$ci_direct_low <- ci_direct[1]
    ans$ci_direct_high <- ci_direct[2]

    ans
}
