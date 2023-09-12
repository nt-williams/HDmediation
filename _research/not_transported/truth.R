# READ ME -----------------------------------------------------------------
#
#       Author: Nick Williams
#       Created: 2023-08-22
#
# Calculates truth for gendata.R DGP
#
# -------------------------------------------------------------------------

library(doFuture)

tmp <- sample_data(1000000)

plan(multisession)

ipw00 <- foreach(a = tmp$a, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
    ipw(a, 0, 0, z, m, w)
}

ipw10 <- foreach(a = tmp$a, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
    ipw(a, 1, 0, z, m, w)
}

ipw11 <- foreach(a = tmp$a, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4)) %dofuture% {
    ipw(a, 1, 1, z, m, w)
}

plan(sequential)

y11 <- mean(ipw11 * tmp$y)
y10 <- mean(ipw10 * tmp$y)
y00 <- mean(ipw00 * tmp$y)
