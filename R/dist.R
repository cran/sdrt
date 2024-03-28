#' Return the distance between two subspaces spanning by column space of matrices.
#'
#' The function calculates three metrics for measuring the distance between two subspaces spaning by the columns of two matrices.
#' @usage dist(A, B)
#'
#' @param A Matrix 1 with dimension p-by-d.
#' @param B Matrix 2 with dimension p-by-d.
#'
#' @return The outputs include three scales and one d-dimensional vector.
#' \item{r}{One minuse the summation of eiegenvalues of the matrix B^TAA^TB.}
#'
#' \item{q}{One minues the product of eiegenvalues of the matrix B^TAA^TB.}
#'
#' \item{rho}{rho=sqrt(A^TBB^TA)}
#'
#' \item{m^2}{A d-variate vector giving the colum-wise distance between A and B.}
#' @export
#'
#' @references
#' Samadi S. Y. and De Alwis T. P. (2023). Fourier Method of Estimating Time Series Central Mean Subspace.
#' \emph{https://arxiv.org/pdf/2312.02110}.
#'
#' Ye Z. and Weiss R.E. (2003). Using the Bootstrap to Select One of a New Class of Dimension Reduction Methods,
#' \emph{Journal of the American Statistical Association}, 98,968-978.
#'
#'
dist <- function(A, B) {
  A.orth <- qr.Q(qr(A))
  B.orth <- qr.Q(qr(B))
  p <- nrow(B.orth)
  d <- ncol(B.orth)

  BAAB <- t(B.orth) %*% A.orth %*% t(A.orth) %*% B.orth
  BAAB.eig <- eigen(BAAB, only.values = T)$values

  rohat <- sqrt(abs(det(t(A.orth) %*% B %*% t(B) %*% A.orth)))
  rohat <- sqrt(abs(det(t(B) %*% A.orth %*% t(A.orth) %*% B)))

  h <- (diag(p) - B.orth %*% t(B.orth)) %*% A.orth
  msq <- apply(h, 2, function(x) sum(x^2))

  list(r = 1 - sqrt(mean((BAAB.eig))), q = 1 - sqrt(prod((BAAB.eig))), ro = rohat, msq = msq)
}
