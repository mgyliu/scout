gridge <- function(x, rho = 0, v = NULL, thetas = NULL, u = NULL) {
  x <- x / sqrt(nrow(x) - 1)
  p <- 2 * rho
  if (is.null(v) || is.null(thetas) || is.null(u)) {
    covarswap <- x %*% t(x)
    eigenswap <- eigen(covarswap)
    keep <- ((eigenswap$values) > 1e-11)
    v <- t(diag(1 / sqrt(eigenswap$values[keep])) %*% t(eigenswap$vectors[, keep]) %*% x)
    thetas <- eigenswap$values[keep]
    u <- eigenswap$vectors[, keep]
  }
  lambda <- sqrt(4 * p) / 2
  diagmat <- (-lambda + (-thetas + sqrt(thetas^2 + 4 * p)) / 2)
  dbar <- .5 * (thetas + sqrt(thetas^2 + 4 * p)) - sqrt(p)
  return(list(svdstuff = list(u = u, v = v, thetas = thetas), wistuff = list(
    v = v, firstdiag = (1 / sqrt(p)),
    diagsandwich = (-((1 / p) * (1 / (1 / sqrt(p) + 1 / dbar))))
  ), wstuff = list(
    v = v, firstdiag = lambda,
    diagsandwich = (thetas + diagmat)
  )))
}
