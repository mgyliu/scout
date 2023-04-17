
scout1something <- function(x, y, p2, lam1s, lam2s, rescale, trace, intercept = FALSE, alternateCov = "default") {
  if (ncol(x) > 500) {
    rlog::log_warn("You are running scout with p1=1 and ncol(x) > 500. This will be slow. You may want to re-start and use p1=2, which is much faster.")
  }

  if (min(lam1s) == 0 && min(lam2s) == 0 && ncol(x) >= nrow(x)) {
    stop("don't run w/lam1=0 and lam2=0 when p>=n")
  }

  if (sum(order(lam2s) == (1:length(lam2s))) != length(lam2s)) {
    stop("Error!!!! lam2s must be ordered!!!")
  }

  # Init a 3D matrix of NAs
  betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))

  if (sum(lam1s > 0 & lam1s < 1e-4) > 0 && ncol(x) >= nrow(x)) {
    rlog::log_warn("Non-zero lam1s that were smaller than 1e-4 were increased to 1e-4 to avoid problems with graphical lasso that can occur when p>n.")
    lam1s[lam1s < 1e-4 & lam1s != 0] <- 1e-4
  }

  # Iterate through lambda 1's
  for (i in 1:length(lam1s)) {
    if (trace) cat(i, fill = F)
    g.out <- NULL

    # On first lambda_1
    # rlog::log_info("Starting GLASSO algorithm")
    if (i == 1 || is.null(g.out$w) || is.null(g.out$wi)) {
      if (lam1s[i] != 0) {
        # Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
        g.out <- glasso::glasso(rob_cov(x, alternateCov = alternateCov), rho = lam1s[i])
      }
      if (lam1s[i] == 0) {
        g.out <- list(w = rob_cov(x, alternateCov = alternateCov), wi = NULL)
      }
    } else if (i != 1 && !is.null(g.out$w) && !is.null(g.out$wi)) {
      g.out <- glasso::glasso(rob_cov(x, alternateCov = alternateCov), rho = lam1s[i], start = "warm", w.init = g.out$w, wi.init = g.out$wi)
    }
    # rlog::log_info("Finished GLASSO algorithm")
    # TODO: maybe estimate the nearest PD matrix to g.out$wi

    # g.out values
    # * w: estimated covariance matrix
    # * wi: estimated inverse covariance matrix

    # Iterate through lambda 2's
    # Step 2 in the algorithm
    for (j in 1:length(lam2s)) {
      if (p2 == 0 || lam2s[j] == 0) {
        # beta = sigma_{xx}^{-1} * sigma_{xy}
        beta <- g.out$wi %*% rob_cov(x, y, alternateCov = alternateCov)
      } else if (p2 == 1 && lam2s[j] != 0) {
        if (j == 1) {
          l_one_res <- lasso_one(g.out$w, rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j])
          beta <- l_one_res$beta
        }
        if (j != 1) {
          if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
            l_one_res <- lasso_one(g.out$w, rob_cov(x, y, alternateCov = alternateCov), rho = lam2s[j], beta.init = beta)
            beta <- l_one_res$beta
            # if got zero for smaller value of lambda 2,
            # then no need to keep computing!!!
          }
        }
      }
      if (rescale && sum(abs(beta)) != 0) {
        # Step 4: rescale beta_hat* = c * beta_hat
        beta <- beta * lsfit(x %*% beta, y, intercept = intercept)$coef
      }
      betamat[i, j, ] <- beta
    }
  }

  list(betamat = betamat, g.out = g.out)
}

scout1something_huge <- function(x, y, p2, lam1s, lam2s, rescale, trace, intercept = FALSE, alternateCov = "default") {
  if (ncol(x) > 500) {
    rlog::log_warn("You are running scout with p1=1 and ncol(x) > 500. This will be slow. You may want to re-start and use p1=2, which is much faster.")
  }

  if (min(lam1s) == 0 && min(lam2s) == 0 && ncol(x) >= nrow(x)) {
    stop("don't run w/lam1=0 and lam2=0 when p>=n")
  }

  if (sum(order(lam2s) == (1:length(lam2s))) != length(lam2s)) {
    stop("Error!!!! lam2s must be ordered!!!")
  }

  # Init a 3D matrix of NAs
  betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))

  if (sum(lam1s > 0 & lam1s < 1e-4) > 0 && ncol(x) >= nrow(x)) {
    rlog::log_warn("Non-zero lam1s that were smaller than 1e-4 were increased to 1e-4 to avoid problems with graphical lasso that can occur when p>n.")
    lam1s[lam1s < 1e-4 & lam1s != 0] <- 1e-4
  }

  # Iterate through lambda 1's
  for (i in 1:length(lam1s)) {
    if (trace) cat(i, fill = F)
    g.out <- NULL

    covx <- rob_cov(x, alternateCov = alternateCov)
    icovx <- NULL
    # On first lambda_1
    # rlog::log_info("Starting GLASSO algorithm")
    if (lam1s[i] != 0) {
      # Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
      hg.out <- huge::huge.glasso(rob_cov(x, alternateCov = alternateCov), lambda = lam1s[i], scr = FALSE, verbose = FALSE)
      # huge.glasso returns a list in $icov of length nlambda. Default nlambda is 10
      # Since we're just passing in one lambda value, get the first item in the list
      icovx <- hg.out$icov[[1]] # Will be (p x p) matrix
    }
    # rlog::log_info("Finished GLASSO algorithm")
    # TODO: maybe estimate the nearest PD matrix to g.out$wi

    # Iterate through lambda 2's
    # Step 2 in the algorithm
    for (j in 1:length(lam2s)) {
      if (p2 == 0 || lam2s[j] == 0) {
        # beta = sigma_{xx}^{-1} * sigma_{xy}
        beta <- icovx %*% rob_cov(x, y, alternateCov = alternateCov)
      } else if (p2 == 1 && lam2s[j] != 0) {
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
      if (rescale && sum(abs(beta)) != 0) {
        # Step 4: rescale beta_hat* = c * beta_hat
        beta <- beta * lsfit(x %*% beta, y, intercept = intercept)$coef
      }
      betamat[i, j, ] <- beta
    }
  }

  list(betamat = betamat, g.out = list(w = covx, wi = icovx))
}
