#' Puts together glasso_select and scout1something
#' @param X n \times p data matrix
#' @param Y n \times 1 response vector
#' @param p2 regularization type in second step. one of NULL or 1
#' @param K number of folds to use in CV. only applicable for p2 = 1
#' @param nlambda1 number of regularization params to use in glasso step
#' @param nlambda2 number of regularization params to use in lasso step
#' @param alternateCov one of "cellwise" or "default"
#' @return list of
#' * `g.res`: result of graphical lasso step
#' * `cv.res`: result of cross-validated lasso step, or NA if p2 = NULL
#' * `betahat`: coefficient estimate#' * `yhat`: predictions on test data
#' * `mod`: final scout model using cross-validated parameters
scout_glasso_stepwise <- function(X, Y, p2 = 1, K = 10,
                                  nlambda1 = 100, lambda1.min.ratio = 0.1,
                                  nlambda2 = 100, lambda2.min.ratio = 0.01,
                                  alternateCov = "default") {
  # First use glasso to find best lambda1
  g.res <- glasso_select(X, nlambda = nlambda1, lambda.min.ratio = lambda1.min.ratio, standardize = TRUE, crit = "bic", alternateCov = alternateCov)

  # Use best lambda from glasso in scout as lambda1
  if (is.null(p2)) {
    cv.res <- NA
    mod <- scout(X, Y, p1 = 1, p2 = p2, lam1s = g.res$best_lambda, alternateCov = alternateCov, trace = F)
  } else if (p2 == 1) {
    cv.res <- cv2.scout(
      X, Y,
      p1 = 1, p2 = p2,
      lam1s = c(g.res$best_lambda),
      nlambda2 = nlambda2,
      lambda2_min_ratio = lambda2.min.ratio,
      K = K, alternateCov = alternateCov
    )
    mod <- scout(X, Y, p1 = 1, p2 = p2, lam1s = g.res$best_lambda, lam2s = cv.res$bestlam2, alternateCov = alternateCov, trace = F)
  } else {
    stop(glue::glue("scout_lasso_stepwise not implemented for p2 = {deparse(p2)}"))
  }

  list(g.res = g.res, cv.res = cv.res, mod = mod)
}

scout_gridge_stepwise <- function(X, Y, p2 = 1, K = 10, nlambda1 = 100, lambda1.min.ratio = 0.1,
                                  nlambda2 = 100, lambda2.min.ratio = 0.01,
                                  alternateCov = "default") {
}
