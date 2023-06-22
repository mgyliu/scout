#' gridge_select
#' @description
#' Uses scout::gridge2 to estimate an inverse *precision* matrix given a
#' data matrix.
#' @param X feature matrix, \eqn{n \times p}
#' @param nlambda number of lambdas to optimize over
#' @param lambda.min.ratio smallest value of lambda as a fraction of lambda_max
#' @param crit criteria to select the optimal lambda. one of "bic" or "loglik"
#' #' @return list of:
#' * `w`: matrix - inverse precision matrix estimate based on best lambda
#' * `best_lambda`: numeric - best lambda selected based on `crit`
#' * `lambda`: numeric vector - sequence of lambdas that was used for selection
#' * `errors`: numeric vector - `crit` values corresponding to each value in
#' `lambda`
#' @export
gridge_select <- function(X, nlambda = 100, lambda.min.ratio = 0.1,
                          crit = "bic", standardize = TRUE) {
  if (standardize) {
    sdx <- apply(X, 2, sd)
    X <- scale(X, T, T)
  } else {
    sdx <- rep(1, ncol(X))
    X <- scale(X, T, F)
  }

  cov_X <- cov(X)
  # Use same lambda sequence as in glasso
  lambda_seq <- sort(huge_glasso_lambda_seq(cov_X, nlambda, lambda.min.ratio))

  gr.out <- gridge2(X, rho = lambda_seq)
  browser()
  errors <- unlist(lapply(gr.out$ws, function(icov_inv) {
    icov_error(solve(icov_inv), cov_X, nrow(X), method = crit)
  }))

  best_idx <- which.min(errors)
  best_lambda <- lambda_seq[best_idx]
  w <- gr.out$ws[[best_idx]]

  list(
    w = w, best_lambda = best_lambda,
    lambda = lambda_seq, errors = errors
  )
}
