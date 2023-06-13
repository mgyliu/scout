# Outer loop: glasso
# Inner loop: lasso_one
# p1 = 1, p2 = 1
#' @param x "training" x data
#' @param y "training" y data
#' @param xtest "testing" x data
#' @param ytest "testing" y data
#' Pseudocode:
# for (lam1 in lam1s) {
#   cv_err_min <- 1e5
#   cv_err_curr <- 0
#   j <- 1
#   while (cv_err_curr <= cv_err_min) {
#     lam2 <- lam2s[j]
#     beta <- find_beta()
#     cv_err_curr <- compute_cv_err(beta, y, ytest)
#     if (cv_err_curr < cv_err_min) {
#       cv_err_min <- cv_err_curr
#     }
#   }
# }
scout1something_alt <- function(x, y, xtest, ytest, p2, lam1s, lam2s, rescale, trace, alternateCov = "default") {

}
