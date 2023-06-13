#' @param X n x p data matrix
#' @param Y n x 1 response vector
#' @param p2 regularization type in second step. one of NULL or 1
#' @param K number of folds to use in CV, only applicable for p2 = 1
#' @param nlambda1 number of regularization params to use in glasso step
#' @param lambda1.min.ratio smallest value of lambda1 as a fraction of lambda1_max
#' @param nlambda2 number of regularization params to use in lasso step
#' @param lambda2.min.ratio smallest value of lambda2 as a fraction of lambda2_max
#' @return list of
#' * `g.res`: result of graphical lasso step
#' * `cv.res`: result of cross-validated lasso step, or NA if p2 = NULL
#' * `mod`: final scout model using cross-validated parameters
scout_1something_stepwise <- function(X, Y, p2, K = 5,
                                      nlambda1 = 100, lambda1.min.ratio = 0.1,
                                      nlambda2 = 100, lambda2.min.ratio = 0.1,
                                      standardize = TRUE) {
  # First use glasso to find best lambda1
  g.res <- glasso_select(
    X,
    nlambda = nlambda1,
    lambda.min.ratio = lambda1.min.ratio,
    standardize = standardize,
    crit = "bic"
  )

  if (is.null(p2)) {
    cv.res <- NA
    mod <- scout(
      X, Y,
      p1 = 1, p2 = p2,
      lam1s = g.res$best_lambda, standardize = standardize
    )
  } else if (p2 == 1) {
    cv.res <- cv.scout(
      X, Y,
      K = K,
      p1 = 1, p2 = p2,
      lam1s = c(g.res$best_lambda),
      nlambda2 = nlambda2,
      lambda2_min_ratio = lambda2.min.ratio,
      standardize = standardize
    )
    mod <- scout(
      X, Y,
      p1 = 1, p2 = p2,
      lam1s = g.res$best_lambda,
      lam2s = cv.res$bestlam2,
      standardize = standardize
    )
  } else {
    stop(glue::glue("scout_lasso_stepwise not implemented for p2 = {deparse(p2)}"))
  }

  list(g.res = g.res, cv.res = cv.res, mod = mod)
}
