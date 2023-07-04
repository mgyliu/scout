#' Computes covariance matrix estimate using either L1 or L2 penalty
#' @param X n x p data matrix. It is assumed to be standardized
#' @param p1 numeric value, either 1 or 2, corresponding to L1 and L2
#' regularization respectively. If p1 = 1, huge::huge.glasso is used to estimate
#' a sparse covariance matrix. If p1 = 2, scout::gridge2 is used to estimate
#' a covariance matrix.
#' @param lam1s numeric vector of penalty parameters for the regularization
#' @return a list of length length(lam1s) where each entry is a covariance
#' matrix estimate. List indices correspond to lam1s indices.
get_xtxs <- function(X, p1, lam1s) {
  if (p1 == 1) {
    g.out <- huge::huge.glasso(cov(X), lambda = lam1s, verbose = FALSE, cov.output = TRUE)
    return(g.out$cov)
  } else if (p1 == 2) {
    g.out <- gridge2(X, rho = lam1s)
    return(g.out$ws)
  } else {
    stop(glue::glue("get_xtxs is not implemented for p1 = {p1}"))
  }
}

#' betas_to_original_scale
#' @description rescales betas to original scale of training data using the
#' standard deviation estimates of the X and Y training data.
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
#' @param sdx a numeric vector. column SDs of original training X data
#' @param sdy a numeric value. SD of original training Y data
#' @return betamat where each column is scaled by sdy/sdx
betas_to_original_scale <- function(betamat, sdx, sdy) {
  sweep(betamat, 1, sdy / sdx, "*")
}

#' compute_intercepts
#' @description computes the intercept term for each beta_hat estimate.
#' The beta matrix should be such that each column contains a different beta_hat
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
#' @param meanx a numeric vector. column means of original training X data
#' @param sdx a numeric vector. column SDs of original training X data
#' @param meany a numeric value. mean of original training Y data
#' @param sdy a numeric value. SD of original training Y data
compute_intercepts <- function(betamat, meanx, sdx, meany, sdy) {
  apply(betamat, MARGIN = 2, FUN = function(beta_hat) {
    meany + sum((sdy * meanx / sdx) * beta_hat)
  })
}

#' compute_rmspe_betamat
#' @description computes the rmspe between the provided Y_test and the predicted
#' X_test %*% beta_hat for each beta_hat in the provided beta matrix.
#' The beta matrix should be such that each column contains a different beta_hat
#' @param X_test n_test x p data matrix (test data, unstandardized)
#' @param Y_test n_test x 1 response vector (test data, unstandardized)
#' @param intercepts numeric vector of length nlambda
#' @param betamat numeric matrix of dimensions p x nlambda
compute_rmspe_betamat <- function(X_test, Y_test, intercepts, betamat) {
  errors <- apply(rbind(intercepts, betamat), MARGIN = 2, FUN = function(beta_hat) {
    yhat <- beta_hat[1] + X_test %*% beta_hat[-1]
    perry::rmspe(yhat, Y_test)
  })
}

#' Rescales beta estimates as in Step 4 of the Scout algo (Witten & Tibshirani)
#' @param X_train standardized training data
#' @param Y_train standardized training response
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
rescale_betas <- function(X_train, Y_train, betamat) {
  apply(betamat, MARGIN = 2, FUN = function(beta_hat) {
    beta_hat * lsfit(X_train %*% beta_hat, Y_train, intercept = FALSE)$coef
  })
}

# assume: p2 = 1
# get_lam2_sequence(cov_x_est)
get_best_lam2_lasso <- function(X_train, Y_train, X_test, Y_test,
                                meanx, meany, sdx, sdy,
                                cov_x_est, lam2s, rescale) {
  # Initialize beta_hat estimates. The beta_hat for a given lambda2 goes into
  # that lambda2's column index.
  betamat <- matrix(NA, nrow = ncol(X_train), ncol = length(lam2s))
  beta <- NA
  for (i in 1:length(lam2s)) {
    if (i == 1) {
      beta <- crossProdLasso(cov_x_est, cov(X_train, Y_train), rho = lam2s[i])$beta
    } else {
      beta <- crossProdLasso(cov_x_est, cov(X_train, Y_train), rho = lam2s[i], beta.init = beta)$beta
    }
    betamat[, i] <- beta
  }

  # Rescale betas if needed (Scout procedure Step 4)
  if (rescale) {
    betamat <- rescale_betas(X_train, Y_train, betamat)
  }

  # Estimate intercepts
  intercepts <- compute_intercepts(betamat, meanx, sdx, meany, sdy)

  # Rescale betas to the original scale of the training data
  betamat <- betas_to_original_scale(betamat, sdx, sdy)

  # Compute yhats based on estimated intercepts, beta matrix, and test data
  errors <- compute_rmspe_betamat(X_test, Y_test, intercepts, betamat)

  # Get best lam2 and corresponding coefficients
  list(
    best_lam2_idx = which.min(errors),
    best_lam2 = lam2s[which.min(errors)],
    intercept = intercepts[which.min(errors)],
    beta_hat = betamat[, which.min(errors)],
    error = min(errors)
  )
}

get_best_lam1 <- function(X_train, Y_train, X_test, Y_test,
                          meanx, meany, sdx, sdy,
                          p1, lam1s, lam2, rescale) {
  xtxs <- get_xtxs(X_train, p1, lam1s)
  # Initialize beta_hat estimates. The beta_hat for a given lambda2 goes into
  # that lambda2's column index.
  betamat <- matrix(NA, nrow = ncol(X_train), ncol = length(lam1s))
  beta <- NA
  for (i in 1:length(lam1s)) {
    if (i == 1) {
      beta <- crossProdLasso(xtxs[[i]], cov(X_train, Y_train), rho = lam2)$beta
    } else {
      beta <- crossProdLasso(xtxs[[i]], cov(X_train, Y_train), rho = lam2, beta.init = beta)$beta
    }
    betamat[, i] <- beta
  }

  # Rescale betas if needed (Scout procedure Step 4)
  if (rescale) {
    betamat <- rescale_betas(X_train, Y_train, betamat)
  }

  # Estimate intercepts
  intercepts <- compute_intercepts(betamat, meanx, sdx, meany, sdy)

  # Rescale betas to the original scale of the training data
  betamat <- betas_to_original_scale(betamat, sdx, sdy)

  # Compute yhats based on estimated intercepts, beta matrix, and test data
  errors <- compute_rmspe_betamat(X_test, Y_test, intercepts, betamat)

  # Get best lam1 and corresponding coefficients
  list(
    best_lam1_idx = which.min(errors),
    best_lam1 = lam1s[which.min(errors)],
    best_xtx = xtxs[[which.min(errors)]],
    intercept = intercepts[which.min(errors)],
    beta_hat = betamat[, which.min(errors)],
    error = min(errors)
  )
}

#' @title scout_alternating_lasso
#' @description implements an alternating search algorithm to efficiently
#' search over a grid of lambda1 x lambda2. L_p penalty for coefficient
#' regularization is assumed to be 1 (i.e., p2 = 1)
#' @param X_train A matrix of predictors. Rows are observations and columns
#' are variables
#' @param Y_train A vector of outcomes
#' @param X_test A matrix of predictors (test set)
#' @param Y_test A vector of outcomes (test set)
#' @param p1 L_p penalty for covariance regularization. must be 1 or 2
#' @param nlambda1 number of regularization terms to use in covariance
#' regularization step
#' @param nlambda2 number of regularization terms to use in coefficient
#' regularization step
#' @param lambda1_min_ratio smallest fraction of lambda1_max to include in
#' lambda1 sequence
#' @param lambda2_min_ratio smallest fraction of lambda2_max to include in
#' lambda2 sequence
#' @param tol convergence tolerance for difference between previous estimated
#' rmspe and current estimated rmspe. if the difference becomes smaller than
#' this value, the algorithm stops and returns the current lam1 and lam2
#' combination
#' @param max_iter maximum number of alternating iterations to compute before
#' stopping. the loop will stop when either the number of iterations exceeds
#' max_iter, or when the errors converge to have a diff below the tolerance
#' @param rescale Should coefficients beta obtained by
#' covariance-regularized regression be re-scaled by a constant, given
#' by regressing $y$ onto $x beta$? This is done in Witten and
#' Tibshirani (2008) and is important for good performance. Default is
#' TRUE.
#' @param standardize whether or not to scale the training X and Y data
#' before performing estimation
#' @param lam1_init one of "random" or "max"
#' @return list with items:
#' * errors: vector of numeric RMSPE values
#' * betas: list of vectors of estimated beta_hats
#' * intercepts: vector of numeric intercept values
#' * lambda_pairs: list of vectors of length 2, each vector is (lam1, lam2)
#' * lambda2_paths: list of vectors of lambda2 paths. user usually doesn't need
#' this but it's provided as a convenience
#' Incides of each list or vector (other than lambda2_paths) correspond to
#' the lambda pair at that index in lambda_pairs
#' @export
scout_alternating_lasso <- function(X_train, Y_train, X_test, Y_test, p1,
                                    nlambda1 = 100, nlambda2 = 100,
                                    lambda1_min_ratio = 0.1, lambda2_min_ratio = 0.001,
                                    tol = 1e-4, max_iter = 10,
                                    rescale = TRUE, standardize = TRUE,
                                    lam1_init = "random") {
  # Can only pass in one lam1_init value
  stopifnot(lam1_init %in% c("random", "max"))

  # Standardize training data if needed
  meany <- mean(Y_train)
  meanx <- apply(X_train, 2, mean)
  if (standardize) {
    sdy <- sd(Y_train)
    sdx <- apply(X_train, 2, sd)
    X_train <- scale(X_train, TRUE, TRUE)
  } else {
    X_train <- scale(X_train, TRUE, FALSE)
    sdx <- rep(1, ncol(X_train))
    sdy <- 1
  }
  Y_train <- (Y_train - meany) / sdy

  # diff keeps track of the difference between current RMSE and previous RMSE
  # niter keeps track of how many loops we did so far
  diff <- tol + 1
  niter <- 0

  # Get lambda 1 sequence
  lam1s <- get_lambda1_path(X_train, p1, nlambda1, lambda1_min_ratio)
  # Compute inital lambda1
  if (lam1_init == "random") {
    lam1 <- lam1s[sample(1:nlambda1, 1)]
  } else {
    lam1 <- lam1s[1]
  }

  # Compute initial covariance estimate
  cov_x_est <- get_xtxs(X_train, p1, lam1)[[1]]

  # More to keep track of: errors, betas, intercepts, lambda pairs
  errors <- c()
  betas <- list()
  intercepts <- c()
  lambda_pairs <- list()
  lam2_paths <- list()

  # Iterate until either:
  # (a) difference in curr and prev errors is <= tol
  # (b) we ran for more than max_iter iterations
  # In each iteration:
  # - get best lam2 for this lam1
  # - get best lam1 for the lam2 in prev step
  # - update covariance estimate for next iteration
  while (diff > tol | niter <= max_iter) {
    niter <- niter + 1

    # Compute best lambda1 for the current cov_x_est
    lam2s <- get_lambda2_path(X_train, Y_train, cov_x_est, 1, nlambda2, lambda2_min_ratio)
    best_lam2_res <- get_best_lam2_lasso(
      X_train, Y_train, X_test, Y_test,
      meanx, meany, sdx, sdy,
      cov_x_est, lam2s, rescale
    )
    lambda_pairs <- c(lambda_pairs, list(c(lam1, best_lam2_res$best_lam2)))
    errors <- c(errors, best_lam2_res$error)
    intercepts <- c(intercepts, best_lam2_res$intercept)
    betas <- c(betas, list(best_lam2_res$beta_hat))
    lam2_paths <- c(lam2_paths, list(lam2s))

    # Update rmspe diff
    if (niter > 1) {
      diff <- abs(diff(tail(errors, 2)))
    }

    if (diff <= tol | niter > max_iter) {
      break
    }

    niter <- niter + 1

    # Update covariance est for next iteration
    best_lam1_res <- get_best_lam1(
      X_train, Y_train, X_test, Y_test,
      meanx, meany, sdx, sdy,
      p1, lam1s, best_lam2_res$best_lam2, rescale
    )
    lambda_pairs <- c(lambda_pairs, list(c(best_lam1_res$best_lam1, best_lam2_res$best_lam2)))
    errors <- c(errors, best_lam1_res$error)
    intercepts <- c(intercepts, best_lam1_res$intercept)
    betas <- c(betas, list(best_lam1_res$beta_hat))
    lam2_paths <- c(lam2_paths, list(lam2s))

    # Update rmspe diff
    if (niter > 1) {
      diff <- abs(diff(tail(errors, 2)))
    }

    cov_x_est <- best_lam1_res$best_xtx
    lam1 <- best_lam1_res$best_lam1
  }

  return(list(
    errors = errors,
    betas = betas,
    intercepts = intercepts,
    lambda_pairs = lambda_pairs,
    lam2_paths = lam2_paths,
    lam1s = lam1s
  ))
}
