#' This predict function is for a single beta, single new_x matrix
#' @param beta
#' @param meany (maybe robust) mean of y data used to fit the model
#' @param new_x X matrix on which to predict
get_yhat <- function(betahat, meanx, meany, sdx, sdy, newx) {
  # Recovered unscaled betahat to compute the intercept
  unscaled_betahat <- betahat / (sdy / sdx)
  # Compute intercept
  icp <- meany - sum((sdy * meanx / sdx) * unscaled_betahat)
  # Compute yhat
  icp + newx %*% betahat
}

predict.scoutobject <- function(object, newx, use_intercept = TRUE, ...) {
  if (class(object) != "scoutobject") stop("Object must be of class 'scoutobject', created by a call to 'scout' function.")
  scout.obj <- object

  # newx has to be a vector or matrix
  if (!is.matrix(newx) & !is.vector(newx)) {
    stop("newx has to be a vector or matrix")
  }

  if (is.vector(newx)) {
    newx <- as.matrix(newx, ncol = 1)
  }

  if (dim(newx)[2] != length(scout.obj$coef[1, 1, ])) {
    stop("newx does not have the correct number of features")
  }

  lam1s <- scout.obj$lam1s
  lam2s <- scout.obj$lam2s

  if (use_intercept) {
    interceptmat <- scout.obj$intercepts
  } else {
    interceptmat <- array(0, dim = dim(scout.obj$intercepts))
  }

  betamat <- scout.obj$coefficients

  yhat <- array(dim = c(length(lam1s), length(lam2s), nrow(newx)))
  for (i in 1:length(lam1s)) {
    for (j in 1:length(lam2s)) {
      yhat[i, j, ] <- interceptmat[i, j] + newx %*% matrix(betamat[i, j, ], ncol = 1)
    }
  }
  return(yhat)
}
