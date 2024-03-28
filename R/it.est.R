it.est <- function(y, p, d, wx = 0.1, space = "mean",std = "FALSE", density = "normal", method = "FM")
{
  sw2 <- wx
  y <- as.matrix(y)
  n <- nrow(y)
  aa <- nrow(y)
  l <- array(0, dim = c((aa - p), p))
  for (i in p:1) {
    l[, (p + 1 - i)] <- y[i:(n - p - 1 + i)]
  }
  x <- l
  yy <- y[(p + 1):aa]
  if (space == "mean") {
    if (method == "FM") {
      if (density == "normal") {
        xy.dr <- FMTSN(yy, x, sw2, std, 0)
      } else if (density == "kernel") {
        xy.dr <- FMTSN(yy, x, sw2, std, 1)
      } else {
        stop("Error! Wrong density provided. Plese check the density and try agin")
      }
    } else {
      stop("Error! Wrong method provided. Plese check the method and try again")
    }
  }
  eta_hat <- xy.dr$evectors[, c(1:d)]
  list(eta_hat = eta_hat, M = xy.dr$M)
}
