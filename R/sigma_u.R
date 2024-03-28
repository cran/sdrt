#' The tuning parameter for the estimation of the time series central mean subspace
#'
#' `sigma_u()' estimates the turning parameter needed to estimate time series central mean subspace in Fourier Method.
#'
#' @usage sigma_u(y, p, d, w1_list=seq(0.1,0.5,by=0.1),space="mean",
#'                                  std=FALSE,density="kernel",method="FM",B=20)
#' @param y A univariate time series observations.
#' @param p Integer value. The lag of the time series.
#' @param d Integer value. The dimension of the time series central mean subspace.
#' @param w1_list (default \{0.1, 0.2,0.3,0.4,0.5\}). The sequence of candidate list for the tuning parameter.
#' @param space (default ``mean''). Specify the SDR subspace needed to be estimated.
#' @param std (default FALSE). If TRUE, then standardizing the time series observations.
#' @param density (default ``kernel''). Specify the density function for the estimation (``kernel'' or ``normal'').
#' @param method (default ``FM''). Specify the estimation method. (``FM'' or ``NW'').
#' @param B (default 20). Number of block bootstrap samples.
#' @return The output is a length(sw2_seq) dimensional vector.
#' \item{dis_sw2}{The average block boostrap distances for each candidate list of values.}
#' @export
#' @references
#' Samadi S. Y. and De Alwis T. P. (2023). Fourier Method of Estimating Time Series Central Mean Subspace.
#' \emph{https://arxiv.org/pdf/2312.02110}.
#' @examples
#' \donttest{
#' data("lynx")
#' y <- log10(lynx)
#' p <- 3
#' d <- 1
#' w1_list=seq(0.1,0.5,by=0.1)
#' Tuning.model=sigma_u(y, p, d, w1_list=w1_list, std=FALSE, B=10)
#' Tuning.model$sigma_u_hat
#' }
#' @importFrom tseries tsbootstrap
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom psych tr
sigma_u <- function(y, p, d, w1_list=seq(0.1,0.5,by=0.1),space="mean",std=FALSE,
                    density="kernel",method="FM",B=20) {
  sigmaw2 <- w1_list
  dj <- matrix(0, nrow = B, ncol = 1)
  dist.r <- matrix(0, nrow = length(sigmaw2), ncol = 1)
  pb <- txtProgressBar(min = 0, max = length(sigmaw2), style = 3)
  y <- as.matrix(y)
  aa <- nrow(y)

  l <- array(0, dim = c((aa - p), p))
  for (i in p:1) {
    l[, (p + 1 - i)] <- y[i:(aa - p - 1 + i)]
  }
  x <- l
  yy <- y[(p + 1):aa]
  y <- yy
  for (j in 1:length(sigmaw2)) {
    sw2 <- sigmaw2[j]
    n <- nrow(y)
    #xy.dr <- FMTS(y, x, sw2, std)
    if (method == "FM") {
      if (density == "normal") {
        xy.dr <- FMTSN(y, x, sw2, std, 0)
      } else if (density == "kernel") {
        xy.dr <- FMTSN(y, x, sw2, std, 1)
      } else {
        stop("Error! Wrong density provided. Plese check the density and try agin")
      }
    }

    s_d <- xy.dr$evectors[, c(1:d)]
    # bootstrapping the original time seires with block size p*p
    y.boost <- tsbootstrap(y, b = round(aa / 2), type = "block", nb = B)

    for (jj in 1:B) {
      y.boostsmpl <- y.boost[, jj]
      aa <- length(y.boostsmpl)
      l <- array(0, dim = c(aa - p, p))
      for (i in p:1) {
        l[, p + 1 - i] <- y.boostsmpl[i:(aa - p - 1 + i)]
      }
      x.boostrap <- l
      y.boostrap <- y.boostsmpl[(p + 1):aa]

      #xy.dr <- FMTS(y.boostrap, x.boostrap, sw2, std)
      if (method == "FM") {
        if (density == "normal") {
          xy.dr <- FMTSN(y.boostrap,x.boostrap, sw2, std, 0)
        } else if (density == "kernel") {
          xy.dr <- FMTSN(y.boostrap, x.boostrap, sw2, std, 1)
        } else {
          stop("Error! Wrong density provided. Plese check the density and try agin")
        }
      }

      s_dj <- xy.dr$evectors[, c(1:d)]
      dist.dj <- dist.space(s_dj, s_d)
      dj[jj] <- dist.dj$r
    }
    dist.r[j, 1] <- mean(dj)
    setTxtProgressBar(pb, j)
  }
  close(pb)
  disttab <- cbind(w1 = w1_list, dbar = dist.r)

  list(dis_sw2 = disttab,sigma_u_hat=disttab[which.min(dist.r),1])
}
