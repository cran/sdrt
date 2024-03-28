#' Select the model parameters using Fourier transformation method.
#'
#' `pd.boots()' estimates the number of lags in the model and the dimension of the time series central mean subspace.
#'
#' @usage pd.boots(y, p_list=seq(2,6,by=1), w1=0.1,  space = "mean",std = FALSE,
#'                                      density = "kernel", method = "FM", B=50)
#' @param y A univariate time series observations.
#' @param p_list (default \{2,3,4,5,6\}). The candidate list of the number of lags, p.
#' @param w1 (default 0.1). The tuning parameter of the estimation.
#' @param space (default ``mean''). Specify the SDR subspace needed to be estimated.
#' @param std (default FALSE). If TRUE, then standardizing the time series observations.
#' @param density (default ``kernel''). Density function for the estimation (``kernel'' or ``normal'').
#' @param method (default ``FM''). Estimation method (``FM'' or ``NW'').
#' @param B (default 50). Number of block bootstrap sample.
#'
#' @return The output is a p-by-p matrix, estimated p and d.
#' \item{dis_dp}{The average block bootsrap distances.}
#'
#' \item{p_hat}{The estimator for p.}
#'
#' \item{d_hat}{The estimator for d.}
#'
#' @export
#' @examples
#' \donttest{
#' data("lynx")
#' y <- log10(lynx)
#' p_list=seq(2,5,by=1)
#' fit.model=pd.boots(y,p_list,w1=0.1,B=10)
#' fit.model$dis_pd
#' fit.model$p_hat
#' fit.model$d_hat
#' }
#' @references
#' Samadi S. Y. and  De Alwis T. P. (2023). Fourier Method of Estimating Time Series Central Mean Subspace.
#' \emph{https://arxiv.org/pdf/2312.02110}.
#'
#' @importFrom tseries tsbootstrap
#' @importFrom utils setTxtProgressBar txtProgressBar
pd.boots <- function(y, p_list=seq(2,6,by=1), w1=0.1,  space = "mean",
                  std = FALSE, density = "kernel", method = "FM", B=50) {
  sw2 <- w1
  ps=length(p_list)
  p_st=p_list[1]
  p_end=p_list[ps]
  if (p_st <= 1 || p_end <= 1 || p_st >= p_end) {
    stop("Error: Check the starting and ending candidate lags (p)")
  } else {
    yy <- y
    yy <- as.matrix(yy)
    aa <- length(yy)
    P <- seq(p_st, p_end, by = 1)
    dj <- matrix(0, nrow = B, ncol = 1)
    dist.r <- matrix(0, nrow = length(P) + 1, ncol = length(P) + 1)

    # create progress bar
    pb <- txtProgressBar(min = 0, max = length(P), style = 3)
    for (j in 1:(length(P))) {
      pp <- P[j]
      for (k in 1:pp) {
        d <- k
        aa <- nrow(yy)
        l <- array(0, dim = c(aa - pp, pp))
        for (i in pp:1) {
          l[, (pp + 1 - i)] <- yy[i:(aa - pp - 1 + i)]
        }
        x <- l
        y <- yy[(pp + 1):aa]

        if (method == "FM") {
          if (density == "normal") {
            xy.dr <- FMTSN(y, x, sw2, std, 0)
          } else if (density == "kernel") {
            xy.dr <- FMTSN(y, x, sw2, std, 1)
          } else {
            stop("Error! Wrong density provided. Plese check the density and try agin")
          }
        }
        s_d = xy.dr$evectors[ , c(1:d)]

        y.boost <- tsbootstrap(yy, b = round(aa / 2), type = "block", nb = B)

        for (jj in 1:B) {
          y.boostsmpl <- y.boost[, jj]
          aa <- length(y.boostsmpl)
          l <- array(0, dim = c(aa - pp, pp))
          for (i in pp:1) {
            l[, (pp + 1 - i)] <- y.boostsmpl[i:(aa - pp - 1 + i)]
          }
          x.boostrap <- l
          y.boostrap <- y.boostsmpl[(pp + 1):aa]


          if (method == "FM") {
            if (density == "normal") {
              xy.dr <- FMTSN(y.boostrap, x.boostrap, sw2, std, 0)
            } else if (density == "kernel") {
              xy.dr <- FMTSN(y.boostrap,x.boostrap, sw2, std, 1)
            } else {
              stop("Error! Wrong density provided. Plese check the density and try agin")
            }
          }
          s_dj <- xy.dr$evectors[, c(1:d)]
          dist.dj <- dist(s_dj, s_d)
          dj[jj] <- dist.dj$r
        }
        dist.r[j, k] <- mean(dj)
        setTxtProgressBar(pb, j)
      }
    }
    close(pb)
    # dist.r=dist.r[,-max(P)]
    dist.r <- dist.r[-max(P), ]

    disttab <- cbind(2:max(P), round(dist.r, 4))
    disttab <- rbind(c("p", 1:(max(P))), disttab)
  }
  dims=dim(dist.r)
  disp=10
  for(k in 1:dims[1]){
    for(r in 1:k){
      if(dist.r[k,r]<disp){
        disp=dist.r[k,r]
        p_hat=k+1
        d_hat=r
      }else{
        disp=disp
      }
    }
  }
  list(dis_pd = noquote(disttab),p_hat=p_hat,d_hat=d_hat)
}
