## ======================================================================##
#   File name: FourierMethod.R                                           #
#                                                                        #
#  R functions that are used for estimating the central subspace         #
#                                                                        #
#  need "FMTS-clibs.dll" (windows)                                       #
#    or "FMTS-clibs.so" (unix/linux)                                     #
#                                                                        #
#  Reference:                                                            #
#      Samadi and Priyan (2020).                                         #
#       Fourier Methods for Estimating the Central Subspace              #
#       in Time series                                                   #
#                                                                        #
# 10/29/2020                                                             #
## ======================================================================##
#' @useDynLib sdrt, .registration = TRUE, .fixes = "C_"
FMN <- function(y, x, sw2, std) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)

  if (std == T) {
    sx.raw <- standardize.x(x)
    sigmahalfinv <- sx.raw$tran.M
    x.std <- sx.raw$x.std
    y <- scale(y)
  } else {
    x.std <- x # scale(x,scale=F)
  }

  M <- .C("vecMfmn", as.double(x.std), as.double(y),
    as.integer(n), as.integer(p), as.double(sw2),
    M = double(p * p)
  )$M

  M <- matrix(M, ncol = p, byrow = F)
  M.raw <- eigen(M, symmetric = T)
  if (std == T) {
    M.evectors <- sigmahalfinv %*% M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    M.evectors <- M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  }
  list(evalues = M.raw$values, evectors = M.evectors, M = M)
}
FMTSN <- function(y, x, sw2, std, den_est) {
  outden <- 0
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  yful <- matrix(0, nrow = n + p, ncol = 1)
  yful[1:n] <- y
  yful[(n + 1):(n + p)] <- x[n, ]

  if (std == T) {
    sx.raw <- standardize.x(x)
    sigmahalfinv <- sx.raw$tran.M
    x.std <- sx.raw$x.std
    varx <- apply(x.std, 2, sd)
    y <- scale(y)
  } else {
    x.std <- scale(x, scale = F)
  }
  varx <- apply(x, 2, sd)

  if (den_est == 1) {
    f0gk <- .C("dlogfmarg", as.double(x.std), as.integer(n),
      as.integer(p), as.double(varx), as.integer(outden),
      f0 = double(n), dlogf = double(n * p)
    )
    cutpoint <- quantile(f0gk$f0, 1e-6) ## trim 10% points.
    index <- (f0gk$f0 >= cutpoint)
    gk <- matrix(f0gk$dlogf, ncol = p)[index, ]
    x.std <- x.std[index, ]
    y <- y[index]
    n <- sum(index)

    M <- .C("vecMfmtscms", as.integer(den_est), as.double(x.std), as.double(y), as.double(yful),
      as.integer(n), as.integer(p), as.double(sw2), as.double(gk),
      M = double(p * p), as.double(varx)
    )$M
  } else {
    gk <- x %*% (solve(t(x) %*% x / n))

    M <- .C("vecMfmtscms", as.integer(den_est), as.double(x.std), as.double(y), as.double(yful),
      as.integer(n), as.integer(p), as.double(sw2), as.double(gk),
      M = double(p * p), as.double(varx)
    )$M
    # print(M)
  }
  M <- matrix(M, ncol = p, byrow = F)
  M.raw <- eigen(M, symmetric = T)
  if (std == T) {
    M.evectors <- sigmahalfinv %*% M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    M.evectors <- M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  }
  list(evalues = M.raw$values, evectors = M.evectors, M = M)
}
# =======================================================================##
#            Estimating the candidate matrix (M_{FMTm})                  #
#                                                                        #
#                                                                        #
# Inputs: y- n-p vector                                                  #
#         x- n-p by p matrix                                             #
#         sw2- tuning parameter 1                                        #
#         std- If TRUE, then data standarized, otherwise not             #
#                                                                        #
# Output: M- candidate matrix to estiamte TS-CMS                         #
#         evalues- Eigenvalues of M.                                     #
#         evectors- Eigenvectors of M.                                   #
## ======================================================================##
FMTS <- function(y, x, sw2, std) {
  x <- as.matrix(x)


  n <- nrow(x)
  p <- ncol(x)

  if (std == T) {
    sx.raw <- standardize.x(x)
    sigmahalfinv <- sx.raw$tran.M
    x.std <- sx.raw$x.std
    y <- scale(y)
  } else {
    x.std <- scale(x, scale = F)
  }
  # sqrt(varx(x))
  varx <- apply(x, 2, sd)
  M <- .C("vecMftsm", as.double(x.std), as.double(y),
    as.integer(n), as.integer(p), as.double(sw2),
    M = double(p * p), as.double(varx)
  )$M

  M <- matrix(M, ncol = p, byrow = F)
  M.raw <- eigen(M, symmetric = T)
  if (std == T) {
    M.evectors <- sigmahalfinv %*% M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  } else {
    M.evectors <- M.raw$vectors
    M.evectors <- apply(M.evectors, 2, function(x) x / sqrt(sum(x^2)))
  }
  list(evalues = M.raw$values, evectors = M.evectors, M = M)
}


standardize.x <- function(x) {
  xe <- eigen(cov(x), symmetric = T)
  Mhalf <- xe$vectors %*% diag(1 / sqrt(xe$values)) %*% t(xe$vectors)
  x.std <- scale(x, scale = F) %*% Mhalf
  list(x.std = x.std, tran.M = Mhalf, Sigmax = cov(x))
}
# =======================================================================##
#            Calculate the distance between two column spaces            #
#            of two matrix A and B.                                      #
#                                                                        #
#          Suppose A and B are two orthonormal matries,                  #
#          then distance define as average of eigenvalues (r),           #
#          or the geometry average (q) of module of matrix B'AA'B.       #
## ======================================================================##
dist.space <- function(A, B) {
  A.orth <- qr.Q(qr(A))
  B.orth <- qr.Q(qr(B))
  p <- nrow(B.orth)
  d <- ncol(B.orth)

  BAAB <- t(B.orth) %*% A.orth %*% t(A.orth) %*% B.orth
  BAAB.eig <- eigen(BAAB, only.values = T)$values

  Dhat <- sqrt(1 - tr(A.orth %*% t(A.orth) %*% B.orth %*% t(B.orth)) / p)
  rohat <- sqrt(abs(det(t(A.orth) %*% B %*% t(B) %*% A.orth)))
  rohat <- sqrt(abs(det(t(B) %*% A.orth %*% t(A.orth) %*% B)))

  h <- (diag(p) - B.orth %*% t(B.orth)) %*% A.orth
  msq <- apply(h, 2, function(x) sum(x^2))

  list(r = 1 - sqrt(mean((BAAB.eig))), q = 1 - sqrt(prod((BAAB.eig))), ro = rohat, eigenval = BAAB.eig, Dhat = Dhat, msq = msq)
}

nw <- function(yt, xt, pp, d, nos) {
  nos <- 10
  rss <- matrix(0, nrow = nos, 1)
  n <- length(yt) # 100#input here: sample size
  X <- xt
  x0 <- yt
  lambda <- 1
  a0 <- matrix(rnorm(pp * d), nrow = pp, ncol = d)
  l <- knd3w(X, x0, a0)
  fmincon.fit <- fmincon(a0,
    fn = l, heq = nlcon11,
    maxfeval = 100, maxiter = 100, tol = 0.01
  )
  a <- matrix(fmincon.fit$par, nrow = pp)
  return(list(eta_hat = a))
}


knd3w <- function(X, x0, a) {
  atx <- X %*% t(a)
  n <- nrow(atx)
  p <- ncol(atx)
  s1 <- sqrt(apply(atx, 2, var))
  ak <- (4 / 5)^(1 / 7)
  an <- matrix(0, ncol = p, nrow = 1)
  gk <- array(0, dim = c(n, n, p))
  WW <- 1
  for (j in 1:p) {
    an[j] <- ak * s1[j] * n^(-1 / 7)
    e <- matrix(1, 1, n)
    x1 <- atx[, j] %*% e
    z <- x1 - t(x1)
    gk[, , j] <- 1 / an[j] / sqrt(2 * pi) * exp(-.5 * (z^2) / an[j])
    WW <- WW * gk[, , j]
  }
  W3 <- matrix(apply(WW, 2, sum), n, 1)
  denow <- as.numeric(t(W3) %*% matrix(1, n, 1))
  W4 <- WW / denow
  function(par) {
    if (sum(par < 0) > 0) {
      fun <- .1000
    } else {
      fun <- as.numeric((t((diag(n) - W4) %*% x0)) %*% ((diag(n) - W4) %*% x0) / n)
    }
    return(fun)
  }
}
nlcon11 <- function(a) {
  b <- size(t(a) %*% a)
  k <- b[1]
  Ceq <- t(a) %*% a - eye(k)
  return(Ceq)
}



## ========================== END OF FILE =================================##
