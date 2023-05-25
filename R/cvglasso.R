#' Computes the lambda sequence for huge::huge.glasso
#' @param cov_X covariance matrix estimate
#' @param nlambda number of lambdas to return
#' @param lambda.min.ratio smallest value of lambda as a fraction of lambda_max
#' @return numeric vector of length nlambda, in decreasing order, of log-spaced
#' lambda values
huge_glasso_lambda_seq <- function(cov_X, nlambda, lambda.min.ratio) {
  p <- ncol(cov_X)
  lambda.max <- max(abs(cov_X - diag(p)))
  lambda.min <- lambda.min.ratio * lambda.max
  exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
}

#' Uses huge::huge.glasso to estimate an inverse covariance matrix given
#' a data matrix.
#' @param X feature matrix, n \times p
#' @param nlambda number of lambdas to optimize over
#' @param lambda.min.ratio smallest value of lambda as a fraction of lambda_max
#' @param alternateCov either "default" or "cellwise". if "default", passes
#' cov(X) to huge.glasso. if "cellwise", passes the cellwise robust covariance
#' estimate of X to huge.glasso.
#' @param standardize whether or not to standardize the columns of X before
#' performing glasso
#' @param crit criteria to select the optimal lambda. one of "bic" or "loglik"
#' @param scr whether to use lossy screening in huge.glasso
#' @param verbose whether to let huge.glasso print progress messages
#' @return list of:
#' * `icov`: matrix - inverse covariance estimate based on best lambda,
#' * `best_lambda`: numeric - best lambda selected based on `crit`
#' * `lambda`: numeric vector - sequence of lambdas that was used for selection
#' * `errors`: numeric vector - `crit` values corresponding to each value in
#' `lambda`
glasso_select <- function(X, nlambda = 10, lambda.min.ratio = 0.1,
                          standardize = TRUE, alternateCov = "default",
                          crit = "bic",
                          scr = F, verbose = F) {
  p <- ncol(X)
  n <- nrow(X)
  if (standardize) {
    std_result_x <- rob_standardize_matrix(X, alternateCov = alternateCov)
    std_x <- std_result_x$x
  } else {
    std_x <- X
  }

  cov_X <- rob_cov(std_x, alternateCov = alternateCov)
  lambda <- huge_glasso_lambda_seq(cov_X, nlambda, lambda.min.ratio)

  hg.out <- huge::huge.glasso(cov_X, lambda = lambda, scr = scr, verbose = verbose)

  # Compute error criteria for each lambda
  errors <- unlist(lapply(hg.out$icov, function(icov) {
    icov_error(icov, cov_X, n, method = crit)
  }))

  best_idx <- which.min(errors)
  best_lambda <- lambda[best_idx]
  icovx <- hg.out$icov[[best_idx]]

  list(
    icovx = icovx, best_lambda = best_lambda,
    lambda = lambda, errors = errors
  )
}

# Assume X is standardized
#' @param X feature matrix
#' @param cvfolds list of CV folds. each list item is a vector of indices to
#' omit in the training data
cvglasso <- function(X, cvfolds, nlambda = 10, lambda.min.ratio = 0.1,
                     alternateCov = "default", crit = "bic",
                     scr = F, verbose = F) {
  stopifnot(crit %in% c("loglik", "bic"))
  p <- ncol(X)
  std_result_x <- rob_standardize_matrix(X, alternateCov = alternateCov)
  cov_X <- rob_cov(std_result_x$x, alternateCov = alternateCov)
  lambda <- huge_glasso_lambda_seq(cov_X, nlambda, lambda.min.ratio)

  cv_errors <- foreach::foreach(i = 1:length(cvfolds), .combine = "rbind") %do% {
    omit <- cvfolds[[i]]
    X_train <- X[-omit, ]
    X_test <- X[omit, ]
    glasso_traintest(X_train, X_test, lambda, alternateCov, crit, scr, verbose)
  }

  best_idx <- which.min(colMeans(cv_errors))
  best_lambda <- lambda[best_idx]
  icovx <- huge::huge.glasso(cov_X, lambda = best_lambda, scr = scr, verbose = verbose)

  list(
    cv_errors = cv_errors,
    lambda = lambda,
    best_lambda = best_lambda,
    icovx = icovx$icov[[1]]
  )
}

#' Compute glasso validation error for a particular train/test split
#' Assume X_train and X_test are NOT already standardized
#' @param X_train unstandardized train data
#' @param X_test unstandardized test data
#' @param lambda lambda squence (numeric vector)
glasso_traintest <- function(X_train, X_test, lambda, alternateCov, crit,
                             scr, verbose) {
  # Train data dimensions
  n <- nrow(X_train)
  p <- ncol(X_train)

  # Sample covariances
  std_result_X_train <- rob_standardize_matrix(X_train, alternateCov = alternateCov)
  cov_X_train <- rob_cov(std_result_X_train$x, alternateCov = alternateCov)
  std_result_X_test <- rob_standardize_matrix(X_test, alternateCov = alternateCov)
  cov_X_test <- rob_cov(X_test, alternateCov = alternateCov)

  # Inverse covariance estimate of train data
  gh.out <- huge::huge.glasso(cov_X_train, lambda = lambda, scr = scr, verbose = verbose)
  # CV error computation
  unlist(lapply(gh.out$icov, function(icov) {
    icov_error(icov, cov_X_test, n, method = crit)
  }))
}

#' @param icov inverse covariance estimate
#' @param cov covariance estimate
#' @param n number of rows in original data matrix
#' @param method one of "loglik", "bic"
icov_error <- function(icov, cov, n, method = "loglik") {
  stopifnot(method %in% c("loglik", "bic"))

  # negative log likelihood = - log |Theta| + tr(Theta * Sigmabrowser()
  neg_loglik <- -determinant(icov, logarithm = TRUE)$modulus[[1]] + sum(diag(icov %*% cov))

  if (method == "bic") {
    # esum is \sum_{i \leq j} \hat{e}_{ij}
    # where \hat{e}_ij = 0 if \Theta_{ij} = 0 and 1 otherwise
    # i.e., count how many unique pairs of variables have non-zero
    # partial correlation with each other (and include the diagonal)
    esum <- sum(abs(icov[lower.tri(icov, diag = T)]) > 1e-8)
    return(neg_loglik + (log(n) / n) * esum)
  }

  return(neg_loglik)
}
