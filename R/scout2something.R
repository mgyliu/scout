
scout2something <- function(x, y, p2, lam1s, lam2s, rescale, trace) {
  if (sum(order(lam2s) == (1:length(lam2s))) != length(lam2s)) {
    stop("Error!!!! lam2s must be ordered!!!")
  }
  if (min(lam1s) == 0 && min(lam2s) == 0 && ncol(x) >= nrow(x)) stop("don't run w/lam1=0 and lam2=0 when p>=n")
  g.out <- NULL
  betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))
  for (i in 1:length(lam1s)) {
    if (trace) cat(i, fill = F)
    if (lam1s[i] != 0) {
      if (i == 1 || is.null(g.out)) g.out <- gridge(x, rho = lam1s[i])
      if (i != 1 && !is.null(g.out)) g.out <- gridge(x, rho = lam1s[i], v = g.out$svdstuff$v, thetas = g.out$svdstuff$thetas, u = g.out$svdstuff$u)
      for (j in 1:length(lam2s)) {
        if (p2 == 0) {
          beta <- diag(rep(g.out$wistuff$firstdiag, ncol(x))) %*% cov(x, y) + g.out$wistuff$v %*% (diag(g.out$wistuff$diagsandwich) %*% ((t(g.out$wistuff$v)) %*% cov(x, y)))
        } else if (p2 != 0 && p2 == 1) {
          if (j == 1) {
            beta <- lasso_one(
              diag(rep(g.out$wstuff$firstdiag, ncol(x))) + g.out$wstuff$v %*% diag(g.out$wstuff$diagsandwich) %*% t(g.out$wstuff$v),
              cov(x, y),
              rho = lam2s[j]
            )$beta
          }
          if (j != 1) {
            if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
              beta <- lasso_one(
                diag(rep(g.out$wstuff$firstdiag, ncol(x))) + g.out$wstuff$v %*% diag(g.out$wstuff$diagsandwich) %*% t(g.out$wstuff$v),
                cov(x, y),
                rho = lam2s[j], beta.init = beta
              )$beta
              # If you got zero for a smaller value of
              # lambda2, then no need to keep computing!!!!!!!
            }
          }
        }
        if (rescale && sum(abs(beta)) != 0) beta <- beta * lsfit(x %*% beta, y, intercept = FALSE)$coef
        betamat[i, j, ] <- beta
      }
    } else if (lam1s[i] == 0) {
      if (p2 == 0) betamat[i, 1, ] <- lsfit(x, y, intercept = FALSE)$coef
      if (p2 == 1) {
        for (j in 1:length(lam2s)) {
          if (lam2s[j] == 0) beta <- lsfit(x, y, intercept = FALSE)$coef
          if (lam2s[j] != 0) beta <- lasso_one(cov(x), cov(x, y), rho = lam2s[j])$beta
          if (sum(abs(beta)) != 0 && rescale) {
            betamat[i, j, ] <- beta * lsfit(x %*% beta, y, intercept = FALSE)$coef
          } else {
            betamat[i, j, ] <- beta
          }
        }
      }
    }
  }
  return(betamat)
}


scout2something_new <- function(x, y, p2, lam1s, lam2s, rescale, trace) {
  betamat <- array(NA, dim = c(length(lam1s), length(lam2s), ncol(x)))

  gr.out <- gridge2(x, rho = lam1s)
  ws <- gr.out$ws

  betamat <- sapply(lam2s, function(lam2) {
    sapply(ws, function(w) {
      if (p2 == 1) {
        lasso_one(w, cov(x, y), rho = lam2)$beta
      }
    }, simplify = "array")
  }, simplify = "array") |>
    aperm(c(2, 3, 1))

  return(betamat)
}
