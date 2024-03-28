
nw.est <- function(y, p, d, n.grid = 10) {
  pp <- p
  y <- as.vector(y)
  n <- length(y)
  aa <- length(y)
  l <- array(0, dim = c((aa - pp), pp))
  for (i in pp:1) {
    l[, (pp + 1 - i)] <- y[i:(n - pp - 1 + i)]
  }
  x <- l
  yy <- y[(p + 1):aa]
  nos <- 10
  rss <- matrix(0, nrow = nos, 1)
  n <- length(yy)
  X <- x
  x0 <- yy
  lambda <- 1
  a <- array(0, dim = c(p, d, n.grid))
  for (jj in 1:n.grid) {
    a0 <- matrix(rnorm(pp * d), nrow = pp, ncol = d)
    l <- knd3w(X, x0, a0)
    fmincon.fit <- fmincon(a0,
      fn = l, heq = nlcon11,
      maxfeval = 100, maxiter = 100, tol = 0.1
    )
    a[, , jj] <- matrix(fmincon.fit$par, nrow = pp)
  }
  eta <- apply(a, c(1, 2), mean)
  return(list(eta_hat = eta))
}


knd3w <- function(X, x0, a) {
  atx <- X %*% a
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
  b <- dim(t(a) %*% a)
  k <- b[1]
  Ceq <- t(a) %*% a - eye(k)
  return(Ceq)
}
