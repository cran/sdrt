#' Estimate the SDR subspaces for univariate time series data.
#'
#' `sdrt()' is the main function to estimate the SDR subspaces in time series.
#' @usage sdrt(y, p, d, w1 = 0.1, space = "mean", std = FALSE,
#'                                  density = "normal", method = "FM",n.grid=10)
#' @param y A univariate time series observations.
#' @param p Integer value. The lag of the time series.
#' @param d Integer value (<p). The dimension of the time series central mean subspace.
#' @param space (default ``mean''). Specify the SDR subspace needed to be estimated.
#' @param method (default ``FM''). Specify the estimation method (``FM'' or ``NW'').
#' @param w1 (default 0.1). The tuning parameter of the ``FM'' estimation method.
#' @param density (default ``kernel''). Specify the density function for the estimation (``kernel'' or ``normal'').
#' @param std (default FALSE). If TRUE, then standardize the data.
#' @param n.grid (default 10). Number of searches for the initial value in ``NW'' method
#' @export
#'
#' @return The output is a p-by-d basis matrix for the TS-CMS.
#'
#' @seealso
#' \code{\link[sdrt]{pd.boots}}, \code{\link[sdrt]{sigma_u}}
#' @references
#' Park J. H., Sriram T. N. and Yin X. (2010). Dimension Reduction in Time Series.
#' \emph{Statistica Sinica}. 20, 747-770.
#'
#' Samadi S. Y. and De Alwis T. P. (2023). Fourier Method of Estimating Time Series Central Mean Subspace.
#' \emph{https://arxiv.org/pdf/2312.02110}.


#' @examples
#' data("lynx")
#' y <- log10(lynx)
#' p <- 3
#' d <- 1
#' fit.model <- sdrt(y, p, d=1,method="FM",density = "kernel")
#'fit.model$eta_hat
#'
#' @importFrom stats cov sd quantile var rnorm
#' @importFrom pracma fmincon eye rand repmat
#' @importFrom tseries tsbootstrap
#' @importFrom psych tr
#'
#'
sdrt <- function(y, p, d, w1 = 0.1, space = "mean",
    std = FALSE, density = "normal", method = "FM",n.grid=10) {

  if (method == "FM") {
    fit <- it.est(y, p, d, wx = w1, space = space, std = std, density = density)
    eta_hat <- fit$eta_hat
  } else if (method == "NW") {
    fit <- nw.est(y, p, d, n.grid = n.grid)
    eta_hat <- fit$eta_hat
  } else {
    stop("Error!, Please check the method")
  }
  list(eta_hat = eta_hat)
}
