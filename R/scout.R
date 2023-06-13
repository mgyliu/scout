# yhat: array of dim (n_lam1 x n_lam2 x n_obs)
# ytrue: vector of length n_obs
# returns: (nlam1 x nlam2) matrix of computed metrics
#          corresponds to the metric for each lam1, lam2 combination
compute_cv_metric <- function(yhat, ytrue, metric = "mse") {
  lengths_match <- all(apply(yhat, c(1, 2), length) == length(ytrue))
  stopifnot("yhat and ytrue dimension mismatch" = lengths_match)

  # For each (lam1, lam2), compute yhat - ytrue
  # resids is (n_lam1 x n_lam2 x n_omit)
  resids <- sweep(yhat, MARGIN = 3, STATS = ytrue, FUN = "-")

  if (metric == "mse") {
    sq_resids <- resids^2
    # For each (lam1, lam2), compute the MSE
    return(apply(sq_resids, c(1, 2), mean))
  } else if (metric == "tausize") {
    return(apply(resids, c(1, 2), pense::tau_size))
  }

  # median absolute prediction error
  else if (metric == "mape") {
    return(apply(resids, c(1, 2), function(r) median(abs(r))))
  }

  # mean absolute prediction error
  else if (metric == "meanape") {
    return(apply(resids, c(1, 2), function(r) mean(abs(r))))
  }

  stop(paste("cv metric", metric, "not yet implemented"))
}

scout <- function(x, y, newx = NULL, p1 = 2, p2 = 1,
                  lam1s = seq(.001, .2, len = 10), lam2s = seq(.001, .2, len = 10),
                  rescale = TRUE, trace = TRUE, standardize = TRUE,
                  rescale_betas = TRUE, alternateCov = "default") {
  call <- match.call()
  if (!is.null(p1) && p1 != 1 && p1 != 2) stop("p1 must be 1, 2, or NULL.")
  if (!is.null(p2) && p2 != 1) stop("p1 must be 1 or NULL.")
  if ((sum(is.na(x)) + sum(is.na(y))) > 0) stop("Please fix the NAs in your data set first. Missing values can be imputed using library 'impute'.")
  x <- as.matrix(x)
  if (min(apply(x, 2, sd)) == 0) stop("Please do not enter an x matrix with variables that are constant.")

  # do standardization if needed
  std_result <- rob_standardize(x, y, standardize, alternateCov)
  meanx <- std_result$meanx
  meany <- std_result$meany
  sdx <- std_result$sdx
  sdy <- std_result$sdy
  x <- std_result$x
  y <- std_result$y

  # Want to re-order lam2s in increasing order
  if (length(lam2s) > 0) {
    lam2s.orig <- lam2s
    s2.out <- sort(lam2s.orig, index.return = TRUE)
    lam2s <- lam2s.orig[s2.out$ix]
  }
  # Done re-ordering lam2s
  # Want to re-order lam1s in increasing order
  if (length(lam1s) > 0) {
    lam1s.orig <- lam1s
    s1.out <- sort(lam1s.orig, index.return = TRUE)
    lam1s <- lam1s.orig[s1.out$ix]
  }
  # Done re-ordering lam1s
  if (is.null(lam1s) || is.null(p1) || sum(lam1s) == 0) {
    p1 <- 0
    lam1s <- c(0)
    lam1s.orig <- c(0)
  }
  if (is.null(lam2s) || is.null(p2) || sum(lam2s) == 0) {
    p2 <- 0
    lam2s <- c(0)
    lam2s.orig <- c(0)
  }
  if (p1 != 0 && p1 != 1 && p1 != 2) stop("p1 must be 1, 2, or 0")
  if (p2 != 0 && p2 != 1) stop("p2 must be 1 or 0")
  if (!is.null(lam1s) && min(lam1s) < 0) stop("lam1s cannot be negative")
  if (!is.null(lam2s) && min(lam2s) < 0) stop("lam2s cannot be negative")
  if (ncol(x) >= nrow(x) && p1 == 0 && p2 == 0) {
    stop("p1 and p2 cannot both be zero when p>n.")
  }

  g.out <- NULL

  # Conditions on p1 start here
  if (p1 == 0) {
    betamat <- array(NA, dim = c(1, length(lam2s), ncol(x)))
    for (j in 1:length(lam2s)) {
      if (trace) cat(j, fill = F)
      if (lam2s[j] == 0) {
        if (ncol(x) >= nrow(x)) stop("Cannot do Least Squares when ncol(x)>=nrow(x)")
        beta <- lsfit(x, y, intercept = FALSE)$coef
        betamat[1, j, ] <- beta
      } else {
        if (j == 1) {
          l_one_res <- lasso_one(rob_cov(x, alternateCov = alternateCov), rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j])
          beta <- l_one_res$beta
          the_cov <- l_one_res$the_cov
        }
        if (j != 1) {
          if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
            l_one_res <- lasso_one(rob_cov(x, alternateCov = alternateCov), rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j], beta.init = beta)
            beta <- l_one_res$beta
            the_cov <- l_one_res$the_cov
            # if got zero for smaller value of lambda 2,
            # then no need to keep computing!!!
          }
        }
        if (rescale && sum(abs(beta)) != 0) {
          beta <- beta * lsfit(x %*% beta, y, intercept = FALSE)$coef
        }
        betamat[1, j, ] <- beta
      }
    }
  } else if (p1 == 1) {
    scout1something_res <- scout1something_huge(x, y, p2, lam1s, lam2s, rescale, trace, alternateCov = alternateCov)
    betamat <- scout1something_res$betamat
    g.out <- scout1something_res$g.out
  } else if (p1 == 2) {
    scout2something_res <- scout2something(x, y, p2, lam1s, lam2s, rescale, trace, alternateCov = alternateCov)
    betamat <- scout2something_res$betamat
    g.out <- scout2something_res$g.out
  }
  interceptmat <- matrix(meany, nrow = length(lam1s), ncol = length(lam2s))
  for (i in 1:length(lam1s)) {
    for (j in 1:length(lam2s)) {
      interceptmat[i, j] <- interceptmat[i, j] - sum((sdy * meanx / sdx) * betamat[i, j, ])
    }
  }

  if (rescale_betas) {
    betamat <- sweep(betamat, 3, sdy / sdx, "*")
  }

  betamat <- array(
    betamat[rank(lam1s.orig), rank(lam2s.orig), ],
    dim = c(length(lam1s.orig), length(lam2s.orig), ncol(x))
  )
  interceptmat <- matrix(interceptmat[rank(lam1s.orig), rank(lam2s.orig)], nrow = length(lam1s.orig), ncol = length(lam2s.orig))
  scout.obj <- (list(
    intercepts = interceptmat,
    coefficients = betamat,
    p1 = p1,
    p2 = p2,
    lam1s = lam1s.orig,
    lam2s = lam2s.orig,
    gout = g.out,
    yhat = NULL,
    call = call
  ))
  class(scout.obj) <- "scoutobject"
  if (!is.null(newx)) {
    yhat <- predict.scoutobject(scout.obj, newx)
    scout.obj$yhat <- yhat
  }
  return(scout.obj)
}

print.scoutobject <- function(x, ...) {
  if (class(x) != "scoutobject") stop("Class of x must be 'scoutobject', created by call to function 'scout'.")
  scout.obj <- x
  cat("Call:\n", fill = F)
  dput(scout.obj$call)
  p1 <- scout.obj$p1
  p2 <- scout.obj$p2
  if (p1 == 0) p1 <- "NULL"
  if (p2 == 0) p2 <- "NULL"
  cat("\n Scout(", p1, ",", p2, ") was performed with lambda1 = (", scout.obj$lam1s, ") and lambda2 = (", scout.obj$lam2s, " ).\n")
  cat("\n Number of non-zero coefficients for each (lambda1, lambda2) pair:\n")
  mat <- apply(x$coef != 0, c(1, 2), sum)
  dimnames(mat) <- list(round(x$lam1s, 3), round(x$lam2s, 3))
  print(mat, quote = F)
  invisible()
}


cv.folds <- function(n, folds = 10) {
  split(sample(1:n), rep(1:folds, length = n))
}

cv.scout <- function(x, y,
                     K = 10,
                     lam1s = seq(0.001, 0.2, len = 10),
                     lam2s = seq(0.001, 0.2, len = 10),
                     p1 = 2, p2 = 1,
                     trace = TRUE, plot = TRUE, plotSE = FALSE,
                     rescale = TRUE, alternateCov = "default", ...) {
  call <- match.call()
  if (K == 1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
  if (K > length(y) / 2) stop("Please choose a value of K between 2 and length(y)/2.")
  if (p1 == 0 && p2 == 0) stop("Why would you want to cross-validate least squares?")
  if (is.null(p1)) lam1s <- 0
  if (is.null(p2)) lam2s <- 0
  lam1s <- c(lam1s)
  lam2s <- c(lam2s)
  if (length(lam1s) < 2 && length(lam2s) < 2) stop("Not a reasonable range of lambdas over which to be cross-validating")
  all.folds <- cv.folds(length(y), K)

  if (length(lam1s) > 1 && length(lam2s) > 1) {
    residmat <- array(0, dim = c(length(lam1s), length(lam2s), K))

    for (i in seq(K)) {
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, rescale = rescale, trace = trace, alternateCov = alternateCov)
      # For each (lam1, lam2), compute yhat - ytrue
      # resids is (n_lam1 x n_lam2 x n_omit)
      resids <- sweep(fit$yhat, MARGIN = 3, STATS = y[omit], FUN = "-")
      sq_resids <- resids^2
      # For each (lam1, lam2), compute the MSE
      residmat[, , i] <- apply(sq_resids, c(1, 2), mean)
    }

    # Now residmat is a (n_lam1 x n_lam2 x K) matrix
    # For each lam1, lam2, and fold, residmat stores the MSPE from evaluating
    # the model on that fold
    #
    # Next line averages out the MSPE over the folds
    cv <- apply(residmat, c(1, 2), mean)
    # TODO: cv.error using `var` which is not robust.
    # use robust measures inside CV here
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
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1, lam2s = lam2s, rescale = rescale, trace = trace, alternateCov = alternateCov)
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
      fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2, trace = trace, rescale = rescale, alternateCov = alternateCov)
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

# Same as cv.scout but without plotting functions
# Shorter because I got rid of the cases depending on number of lambda1 and lambda2
# The result is exactly the same
cv2.scout <- function(x, y, p1 = 2, p2 = 1,
                      K = 10,
                      # Lambda 1 arguments
                      nlambda1 = 100, lambda1_min_ratio = 0.001,
                      lamda1_max_desc_frac = 0.8, lambda1_max_iter = 100,
                      lam1s = NULL,
                      # Lambda 2 arguments
                      nlambda2 = 100, lambda2_min_ratio = 0.001,
                      lamda2_max_desc_frac = 0.8, lambda2_max_iter = 100,
                      lam2s = NULL,
                      trace = FALSE,
                      rescale = TRUE,
                      alternateCov = "default",
                      cvmetric = "mse",
                      ncores = 1,
                      ...) {
  call <- match.call()
  if (K == 1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
  if (K > length(y) / 2) stop("Please choose a value of K between 2 and length(y)/2.")
  if (p1 == 0 && p2 == 0) stop("Why would you want to cross-validate least squares?")

  std_result <- rob_standardize(x, y, alternateCov = alternateCov)

  if (is.null(lam1s)) {
    lam1s <- get_lambda1_path(
      std_result$x, p1,
      alternateCov = alternateCov,
      nlambda = nlambda1,
      lambda_min_ratio = lambda1_min_ratio
    )
  }
  if (is.null(lam2s)) {
    lam2s <- get_lambda2_path(
      std_result$x, std_result$y, p2,
      alternateCov = alternateCov,
      nlambda = nlambda2,
      lambda_min_ratio = lambda2_min_ratio
    )
  }

  if (length(lam1s) < 2 && length(lam2s) < 2) stop("Not a reasonable range of lambdas over which to be cross-validating")
  all.folds <- cv.folds(length(y), K)

  residmat <- array(0, dim = c(length(lam1s), length(lam2s), K))

  for (i in seq(K)) {
    if (trace) cat("\n CV Fold", i, "\t")
    omit <- all.folds[[i]]
    fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, rescale = rescale, trace = trace, alternateCov = alternateCov)

    residmat[, , i] <- compute_cv_metric(fit$yhat, y[omit], cvmetric)
  }

  # Now residmat is a (n_lam1 x n_lam2 x K) matrix
  # For each lam1, lam2, and fold, residmat stores the MSPE from evaluating
  # the model on that fold

  # Next line averages out the MSPE over the folds
  cv <- apply(residmat, c(1, 2), mean)

  cv.error <- sqrt(apply(residmat, c(1, 2), var) / K)

  object <- list(
    p1 = p1,
    p2 = p2,
    lam1s = lam1s,
    lam2s = lam2s,
    cv = cv,
    cv.error = cv.error,
    call = call,
    folds = all.folds,
    bestlam1 = lam1s[which.min(apply(cv, 1, min))],
    bestlam2 = lam2s[which.min(apply(cv, 2, min))]
  )

  if (trace) cat("\n")
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
