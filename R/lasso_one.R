#' crossProdLasso
#'
#' @title
#' Performs the lasso on the cross product matrices X'X and X'y
#'
#' @description
#' Perform L1-regularized regression of y onto X using only the cross-product
#' matrices X'X and X'y. In the case of covariance-regularized regression, this
#' is useful if you would like to try out something other than L1 or L2
#' regularization of the inverse covariance matrix.
#'
#' Suppose you use your own method to regularize X'X.
#' Then let Sigma denote your estimate of the population covariance matrix.
#' Now say you want to minimize
#'    beta' Sigma  beta - 2 beta' X'y + lambda ||beta||_1
#' in order to get the regression estimate beta, which maximizes the
#' second scout criterion when an L_1 penalty is used. You can do this by
#' calling crossProdLasso(Sigma, X'y, rho).
#'
#' If you run crossProdLasso(X'X, X'y, rho) then it should give the same
#' result as lars(X, y)
#'
#' Notice that the X'X that you pass into this function must be POSITIVE
#' SEMI DEFINITE (or positive definite) or the problem is not convex and
#' the algorithm will not converge.
#'
#' @details
#' If your xtx is simply X'X for some X, and your xty is simple X'y
#' with some y, then the results will be the same as running lars on data
#' (X,y) for a single shrinkage parameter value.
#'
#' Note that when you use the scout function with p2=1, the crossProdLasso
#' function is called internally to give the regression coefficients, after the
#' regularized inverse covariance matrix is estimated. It is provided
#' here in case it is useful to the user in other settings.
#'
#' @param xtx A pxp matrix, which should be an estimate of a covariance
#' matrix. This matrix must be POSITIVE SEMI DEFINITE (or positive definite)
#' or the problem is not convex and the algorithm will not converge
#' @param xty A px1 vector, which is generally obtained via X'y.
#' @param rho Must be non-negative; the regularization parameter you are using
#' @param thr Convergence threshold.
#' @param maxit How many iterations to perform?
#' @param beta.init If you're running this over a range of rho values,
#' then set beta.init equal to the solution you got for a previous rho
#' value. It will speed things up.
#'
#' @return list(beta) A px1 vector with the regression coefficients
crossProdLasso <- function(xtx, xty, rho, thr = 1e-4, maxit = 100, beta.init = NULL) {
  return(list(beta = lasso_one(xtx, xty, rho, thr, maxit, beta.init)$beta))
}

lasso_one <- function(w, ss, rho, thr = 1.0e-4, maxit = 100, trace = F, beta.init = NULL) {
  # does lasso fit of a single response variable  on predictors,
  # via coordinate descent. Data is inner products w and ss

  n <- length(ss)

  if (length(rho) == 1) {
    rho <- rep(rho, n)
  }

  itrace <- 1 * trace

  if (is.null(beta.init)) {
    beta.init <- rep(0, n)
  }
  mode(rho) <- "single"
  mode(ss) <- "single"
  mode(w) <- "single"
  mode(n) <- "integer"
  mode(maxit) <- "integer"
  mode(itrace) <- "integer"
  mode(thr) <- "single"
  mode(beta.init) <- "single"

  junk <- .Fortran("lasso7", rho, n, as.matrix(w), ss, thr, xx = beta.init, PACKAGE = "scout")

  return(list(beta = junk$xx))
}
