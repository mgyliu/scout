huge_glasso_lambda_seq <- function(cov_X, p, nlambda, lambda.min.ratio) {
  lambda.max <- max(abs(cov_X - diag(p)))
  lambda.min <- lambda.min.ratio * lambda.max
  exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
}

glasso_select <- function(X, nlambda = 10, lambda.min.ratio = 0.1,
                          alternateCov = "default", crit = "bic",
                          scr = F, verbose = F) {
  p <- ncol(X)
  n <- nrow(X)
  std_result_x <- rob_standardize_matrix(X, alternateCov = alternateCov)
  cov_X <- rob_cov(std_result_x$x, alternateCov = alternateCov)
  lambda <- huge_glasso_lambda_seq(cov_X, p, nlambda, lambda.min.ratio)

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
#' @param X unstandardized feature matrix
#' @param cvfolds list of CV folds. each list item is a vector of indices to
#' omit in the training data
cvglasso <- function(X, cvfolds, nlambda = 10, lambda.min.ratio = 0.1,
                     alternateCov = "default", crit = "bic",
                     scr = F, verbose = F) {
  stopifnot(crit %in% c("loglik", "bic"))
  p <- ncol(X)
  std_result_x <- rob_standardize_matrix(X, alternateCov = alternateCov)
  cov_X <- rob_cov(std_result_x$x, alternateCov = alternateCov)
  lambda <- huge_glasso_lambda_seq(cov_X, p, nlambda, lambda.min.ratio)

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
