predict.scoutobject <- function(object, newx, ...) {
  if (class(object) != "scoutobject") stop("Object must be of class 'scoutobject', created by a call to 'scout' function.")
  scout.obj <- object
  if (!is.matrix(newx) && length(newx) != length(scout.obj$coef[1, 1, ])) stop("newx should be a single observation in vector form, or a matrix of observations. If in matrix form, then there should be one observation per row, with the features on the columns, just like the x matrix.")
  lam1s <- scout.obj$lam1s
  lam2s <- scout.obj$lam2s
  interceptmat <- scout.obj$intercepts
  betamat <- scout.obj$coefficients
  if (is.matrix(newx)) yhat <- array(dim = c(length(lam1s), length(lam2s), nrow(newx)))
  if (!is.matrix(newx)) yhat <- matrix(0, nrow = length(lam1s), ncol = length(lam2s))
  for (i in 1:length(lam1s)) {
    for (j in 1:length(lam2s)) {
      if (is.matrix(newx)) yhat[i, j, ] <- interceptmat[i, j] + newx %*% matrix(betamat[i, j, ], ncol = 1)
      if (!is.matrix(newx)) yhat[i, j] <- interceptmat[i, j] + matrix(newx, nrow = 1) %*% matrix(betamat[i, j, ], ncol = 1)
    }
  }
  return(yhat)
}
