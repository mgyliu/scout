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
#' @param crit criteria to select the optimal lambda. one of "bic" or "loglik"
#' @param scr whether to use lossy screening in huge.glasso
#' @param verbose whether to let huge.glasso print progress messages
#' @return list of:
#' * `icov`: matrix - inverse covariance estimate based on best lambda,
#' * `best_lambda`: numeric - best lambda selected based on `crit`
#' * `lambda`: numeric vector - sequence of lambdas that was used for selection
#' * `errors`: numeric vector - `crit` values corresponding to each value in
#' `lambda`
glasso_select <- function(X, nlambda = 100, lambda.min.ratio = 0.1,
                          crit = "bic", standardize = T, scr = F, verbose = F) {
  p <- ncol(X)
  n <- nrow(X)

  # Need to center and scale X
  meanx <- apply(X, 2, mean)
  if (standardize) {
    sdx <- apply(X, 2, sd)
    X <- scale(X, T, T)
  } else {
    X <- scale(X, T, F)
    sdx <- rep(1, ncol(X))
  }

  cov_X <- cov(X)
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

#' @param icov inverse covariance estimate
#' @param cov covariance estimate
#' @param n number of rows in original data matrix
#' @param method one of "loglik", "bic"
icov_error <- function(icov, cov, n, method = "loglik") {
  stopifnot(method %in% c("loglik", "bic"))

  # negative log likelihood = - log |Theta| + tr(Theta * Sigma)
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
