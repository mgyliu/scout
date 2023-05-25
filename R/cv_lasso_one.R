# Assuming p2 = 1
cv_lasso_one <- function(x, y, lam2s) {
  for (j in 1:length(lam2s)) {
    if (j == 1) {
      l_one_res <- lasso_one(covx, rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j])
      beta <- l_one_res$beta
    }
    if (j != 1) {
      if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
        l_one_res <- lasso_one(covx, rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j], beta.init = beta)
        beta <- l_one_res$beta
        # if got zero for smaller value of lambda 2,
        # then no need to keep computing!!!
      }
    }
  }
}
