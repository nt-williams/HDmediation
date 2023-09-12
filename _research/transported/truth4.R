# READ ME -----------------------------------------------------------------
#
#       Author: Nick Williams
#       Created: 2023-08-24
#
# -------------------------------------------------------------------------

library(doFuture)

tmp <- sample_data(1e6)

plan(multisession)

ipw00 <- foreach(a = tmp$a, s = tmp$s, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4,
                                                      seed = TRUE)) %dofuture% {
                                                          ipw(a, 0, 0, s, z, m, w)
                                                      }

ipw10 <- foreach(a = tmp$a, s = tmp$s, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4,
                                                      seed = TRUE)) %dofuture% {
                                                          ipw(a, 1, 0, s, z, m, w)
                                                      }

ipw11 <- foreach(a = tmp$a, s = tmp$s, z = tmp$z, m = tmp$m, w = tmp$w,
                 .combine = c, .options.future = list(chunk.size = 1e4,
                                                      seed = TRUE)) %dofuture% {
                                                          ipw(a, 1, 1, s, z, m, w)
                                                      }

plan(sequential)

y11 <- mean(ipw11 * tmp$y)
y10 <- mean(ipw10 * tmp$y)
y00 <- mean(ipw00 * tmp$y)

# indirect <- 0.0146
# direct <- 0.0266
