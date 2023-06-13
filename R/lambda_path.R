#' Get the lambda1 path dynamically, from the data. The data is assumed to
#' be standardized already.
#'
#' Since scout only uses lambda1 to regularize the XX precision matrix,
#' this function only needs the X data
#'
#' @param X data matrix (n x p)
#' @param p1 type of regularization in first scout step. Can be NULL or 1.
#' Ignored if greater than 1
#' @param nlambda how many lambda1's do you want to regularize over?
#' @param lambda_min_ratio smallest value for lambda1 as a fraction of
#' lambda1_max.
#' @param desc_frac by how much to multiply lambda_max when tuning it in
#' bisection search for the "ideal" lambda_max
#' @param max_iter maximum number of iterations to use to find the lowest
#' lambda_max in bisection search
get_lambda1_path <- function(X, p1,
                             nlambda = 100, lambda_min_ratio = 0.001,
                             desc_frac = 0.8, max_iter = 100) {
  if (is.null(p1)) {
    return(0)
  } else if (p1 == 1) {
    lambda_max <- get_lambda1_max_glasso(X, desc_frac, max_iter)
    return(get_lambda_path(lambda_max, nlambda, lambda_min_ratio))
  } else {
    warning(paste0("get_lambda1_path not implemented for p1 == ", p1, ". Returning default sequence which may not work the best."))
    return(seq(0.01, 0.2, len = nlambda))
  }
}

#' Get the lambda2 path dynamically, from the data. The data is assumed to
#' be standardized already
#'
#' This function does not depend on lambda1. Empirically, the lowest possible
#' lambda2_max is the same regardless of lambda1. lambda1 controls the cov_x
#' value (coming form GLASSO). But we are just using a regular cov_x estimator
#' and the lambda2_max is the same.
#'
#' @param X data matrix (n x p)
#' @param Y data matrix (n x 1)
#' @param p2 type of regularization in second scout step. Can be NULL or 1.
#' Ignored if greater than 1
#' @param cov_x estimate for covariance matrix of X. This should be the $w
#' portion of the glasso::glasso output. This argument is needed since
#' lambda2_max depends on the step 1 result using lambda1
#' @param nlambda how many lambda2's do you want to regularize over?
#' @param lambda_min_ratio smallest value for lambda2 as a fraction of
#' lambda2_max.
#' @param desc_frac by how much to multiply lambda_max when tuning it in
#' bisection search for the "ideal" lambda_max
#' @param max_iter maximum number of iterations to use to find the lowest
#' lambda_max in bisection search
get_lambda2_path <- function(X, Y, p2,
                             nlambda = 100, lambda_min_ratio = 0.001,
                             desc_frac = 0.8, max_iter = 100) {
  if (is.null(p2)) {
    return(0)
  } else if (p2 == 1) {
    lambda_max <- get_lambda2_max_lasso(X, Y)
    return(get_lambda_path(lambda_max, nlambda, lambda_min_ratio))
  } else {
    warning(paste0("get_lambda2_path not implemented for p2 == ", p2, ". Returning default sequence which may not work the best."))
    return(seq(0.001, 0.2, len = nlambda))
  }
}

is_off_diag_zero <- function(mat) {
  diag(mat) <- 0
  all(mat == 0)
}

#' Smallest value of lambda1 which makes cov(X,X) have all zero off-diagonal
#' entries. The diagonal matrix will have diagonal entries
#' 1/(lambda1_max + x_i^T * x_i)
get_lambda1_max_glasso <- function(X,
                                   desc_frac = 0.8, max_iter = 100,
                                   use_huge = TRUE) {
  XTX <- abs(t(X) %*% X)
  diag(XTX) <- 0
  lambda1_max_tmp <- max(XTX)
  lambda1_max <- lambda1_max_tmp
  iter <- 0
  cond <- TRUE
  # lambda1_max is the last lambda1_max that makes all off-diagonal elements equal to zero
  while (cond & iter < max_iter) {
    covx <- cov(X)

    cond <- if (use_huge) {
      hg.out <- huge::huge.glasso(covx, lambda = lambda1_max_tmp, verbose = FALSE)
      is_off_diag_zero(hg.out$icov[[1]])
    } else {
      gg.out <- glasso::glasso(covx, rho = lambda1_max_tmp)
      is_off_diag_zero(gg.out$wi)
    }

    if (cond) {
      lambda1_max <- lambda1_max_tmp
      lambda1_max_tmp <- lambda1_max_tmp * desc_frac
    }
    iter <- iter + 1
  }

  lambda1_max
}

get_lambda2_max_lasso <- function(X, Y, desc_frac = 0.8, max_iter = 100) {
  lambda2_max_tmp <- max(abs(t(X) %*% Y) * 2)
  lambda2_max <- lambda2_max_tmp
  iter <- 0
  cond <- TRUE
  while (cond & iter < max_iter) {
    l_one_res <- lasso_one(
      cov(X),
      cov(X, Y),
      rho = lambda2_max_tmp
    )
    cond <- all(l_one_res$beta == 0)
    if (cond) {
      lambda2_max <- lambda2_max_tmp
      lambda2_max_tmp <- lambda2_max_tmp * desc_frac
    }
    iter <- iter + 1
  }
  # rlog::log_info(paste0("In get_lambda2_max_glasso: Number of iterations: ", iter, ". max_iter was ", max_iter))
  lambda2_max
}

get_lambda_path <- function(lambda_max, nlambda, lambda_min_ratio) {
  lambdapath <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio),
    length.out = nlambda
  ))
  digits <- 10
  lambdapath
}
