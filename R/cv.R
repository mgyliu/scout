cv.folds <- function(n, folds = 10) {
  split(sample(1:n), rep(1:folds, length = n))
}

#' cv.scout
#' @title
#' Perform cross-validation for covariance-regularized regression, aka the Scout.
#' @description
#' This function returns cross-validation error rates for a range
#' of lambda1 and lambda2 values, and also makes beautiful CV plots if
#' plot=TRUE.
#' @param x A matrix of predictors, where the rows are the samples and
#' the columns are the predictors
#' @param y A matrix of observations, where length(y) should equal nrow(x)
#' @param KNumber of cross-validation folds to be performed; default is 10
#' @param lam1s The (vector of) tuning parameters for regularization of the
#' covariance matrix. Can be NULL if p1=NULL, since then no covariance
#' regularization is taking place. If p1=1 and nrow(x)<ncol(x), then the no value in lam1s
#' should be smaller than 1e-3, because this will cause graphical lasso
#' to take too long. Also, if ncol(x)>500 then we really do not
#' recommend using p1=1, as graphical lasso can be uncomfortably slow.
#' @param lam2s The (vector of) tuning parameters for the $L_1$ regularization of
#' the regression coefficients, using the regularized covariacne
#' matrix. Can be NULL if p2=NULL. (If p2=NULL, then non-zero lam2s
#' have no effect). A value of 0 will result in no
#' regularization.
#' @param p1 The $L_p$ penalty for the covariance regularization. Must be
#' one of 1, 2, or NULL. NULL corresponds to no covariance
#' regularization.
#' @param p2 The $L_p$ penalty for the estimation of the regression
#' coefficients based on the regularized covariance matrix. Must be one
#' of 1 (for $L_1$ regularization) or NULL (for no regularization).
#' @param trace Print out progress as we go? Default is TRUE.
#' @param plot If TRUE (by default), makes beautiful CV plots.
#' @param plotSE Should those beautiful CV plots also display std error
#'     bars for the CV? Default is FALSE
#' @param rescale Scout rescales coefficients, by default, in order to
#'     avoid over-shrinkage
#' @param ... Additional parameters
#' @return
#' - `folds`: The indices of the members of the K test sets are
#'   returned.
#' - `cv`: A matrix of average cross-validation errors is returned.
#' - `cv.error`: A matrix containing the standard errors of the
#'   elements in "cv", the matrix of average cross-validation errors.
#' - `bestlam1`: Best value of lam1 found via cross-validation.
#' - `bestlam2`: Best value fo lam2 found via cross-validation.
#' - `lam1s`: Values of lam1 considered.
#' - `lam2s`: Values of lam2 considered.
#' @export
cv.scout <- function(x, y, K = 10,
                     nlambda1 = 100, lambda1_min_ratio = 0.01,
                     nlambda2 = 100, lambda2_min_ratio = 0.001,
                     lam1s = NULL, # seq(0.001, .2, len = 10),
                     lam2s = NULL, # seq(0.001, .2, len = 10),
                     p1 = 2, p2 = 1,
                     trace = FALSE, plot = FALSE, plotSE = FALSE, rescale = TRUE, ...) {
  call <- match.call()
  if (K == 1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
  if (K > length(y) / 2) stop("Please choose a value of K between 2 and length(y)/2.")
  if (p1 == 0 && p2 == 0) stop("Why would you want to cross-validate least squares?")

  x_std <- scale(x, TRUE, TRUE)
  y_std <- (y - mean(y)) / sd(y)
  if (is.null(lam1s)) {
    lam1s <- get_lambda1_path(
      x_std, p1,
      nlambda = nlambda1,
      lambda_min_ratio = lambda1_min_ratio
    )
  }
  if (is.null(lam2s)) {
    lam2s <- get_lambda2_path(
      x_std, y_std, p2,
      nlambda = nlambda2,
      lambda_min_ratio = lambda2_min_ratio
    )
  }

  if (length(lam1s) < 2 && length(lam2s) < 2) stop("Not a reasonable range of lambdas over which to be cross-validating")
  all.folds <- cv.folds(length(y), K)
  if (length(lam1s) > 1 && length(lam2s) > 1) {
    residmat <- array(0, dim = c(length(lam1s), length(lam2s), K))
    for (i in seq(K)) {
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, rescale = rescale, trace = trace)
      residmat[, , i] <- apply((sweep(fit$yhat, 3, y[omit], "-"))^2, c(1, 2), mean)
    }
    cv <- apply(residmat, c(1, 2), mean)
    cv.error <- sqrt(apply(residmat, c(1, 2), var) / K)
    object <- list(p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1s[which.min(apply(cv, 1, min))], bestlam2 = lam2s[which.min(apply(cv, 2, min))])
    if (plot) {
      myrainbow <- rainbow(length(lam1s) * 1.15)
      if (!plotSE) plot(0, 0, xlim = range(lam2s), ylim = range(cv), col = "white", main = "CV Error", xlab = "lam2", ylab = "CV Error")
      if (plotSE) plot(0, 0, xlim = range(lam2s), ylim = range(c(cv, cv + cv.error, cv - cv.error)), col = "white", main = "CV Error", xlab = "lam2", ylab = "CV Error")
      for (i in 1:length(lam1s)) {
        lines(lam2s, cv[i, ], col = myrainbow[i])
        points(lam2s, cv[i, ], col = myrainbow[i])
        if (plotSE) {
          points(lam2s, cv[i, ] + cv.error[i, ], col = myrainbow[i], type = "l", lty = "dashed")
          points(lam2s, cv[i, ] - cv.error[i, ], col = myrainbow[i], type = "l", lty = "dashed")
        }
      }
      legend("topright", pch = 15, col = myrainbow[1:length(lam1s)], paste("lam1=", sep = "", round(lam1s, 4)))
    }
  } else if (length(lam1s) == 1) {
    lam1 <- lam1s[1]
    residmat <- matrix(0, nrow = length(lam2s), ncol = K)
    for (i in seq(K)) {
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1, lam2s = lam2s)
      residmat[, i] <- apply(sweep(fit$yhat[1, , ], 2, y[omit], "-")^2, 1, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var) / K)
    object <- list(p1 = p1, p2 = p2, lam1s = lam1, lam2s = lam2s, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1, bestlam2 = lam2s[which.min(cv)])
    if (plot) {
      if (!plotSE) plot(0, 0, xlim = range(lam2s), ylim = range(cv), col = "white", main = "CV Error", xlab = "lam2", ylab = "CV Error")
      if (plotSE) plot(0, 0, xlim = range(lam2s), ylim = range(c(cv, cv + cv.error, cv - cv.error)), col = "white", main = "CV Error", xlab = "lam2", ylab = "CV Error")
      lines(lam2s, cv, col = "red")
      points(lam2s, cv, col = "red")
      legend("topright", pch = 15, col = "red", paste("lam1=", sep = "", lam1))
      if (plotSE) {
        points(lam2s, cv + cv.error, col = "red", type = "l", lty = "dashed")
        points(lam2s, cv - cv.error, col = "red", type = "l", lty = "dashed")
      }
    }
  } else if (length(lam2s) == 1) {
    lam2 <- lam2s[1]
    residmat <- matrix(0, nrow = length(lam1s), ncol = K)
    for (i in seq(K)) {
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2)
      residmat[, i] <- apply(sweep(fit$yhat[, 1, ], 2, y[omit], "-")^2, 1, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var) / K)
    object <- list(p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1s[which.min(cv)], bestlam2 = lam2)
    if (plot) {
      if (!plotSE) plot(0, 0, xlim = range(lam1s), ylim = range(cv), col = "white", main = "CV Error", xlab = "lam1", ylab = "CV Error")
      if (plotSE) plot(0, 0, xlim = range(lam1s), ylim = range(c(cv, cv + cv.error, cv - cv.error)), col = "white", main = "CV Error", xlab = "lam1", ylab = "CV Error")
      lines(lam1s, cv, col = "red")
      points(lam1s, cv, col = "red")
      if (plotSE) {
        points(lam1s, cv + cv.error, col = "red", type = "l", lty = "dashed")
        points(lam1s, cv - cv.error, col = "red", type = "l", lty = "dashed")
      }
      legend("topright", pch = 15, col = "red", paste("lam2=", sep = "", lam2))
    }
  }
  if (trace) cat("\n")
  object$call <- call
  object$folds <- all.folds
  class(object) <- "cvobject"
  invisible(object)
}

print.cvobject <- function(x, ...) {
  if (class(x) != "cvobject") stop("Class of x must be 'cvobject', created by call to cv.scout.")
  cat("Call:\t", fill = F)
  dput(x$call)
  cat("\n Cross-validation MSE for each lambda1/lambda2 pair: \n")
  mat <- matrix(round(x$cv, 2), nrow = length(x$lam1s), ncol = length(x$lam2s))
  dimnames(mat) <- list(round(x$lam1s, 3), round(x$lam2s, 3))
  print(mat, quote = F)
  invisible()
}
