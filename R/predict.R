#' predict.scoutobject
#' @description
#' A function to perform prediction, using an x matrix and the output of
#' the "scout" function.
#' @param object The results of a call to the "scout" function. The
#' coefficients that are part of this object will be used for
#' making predictions
#' @param newx The new x at which predictions should be made. Can be a
#' vector of length ncol(x), where x is the data on which scout.obj was
#' created, or a matrix with ncol(x) columns
#' @return
#' `yhat`: If newx was a vector, then a  matrix will be returned,
#' with dimension length(lam1s)xlength(lam2s) (where lam1s and lam2s
#' are attributes of scout.obj). The (i,j) element of this matrix will
#' correspond to tuning parameter values (lam1s[i], lam2s[j]). If newx
#' is a matrix, then an array of dimension
#' nrow(newx) x length(lam1s) x length(lam2s) will be returned.
#' @param ... Additional arguments to predict
#' @export
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
