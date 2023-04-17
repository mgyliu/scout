test_that("scout1something runs for multiple lam1s and lam2s", {
  nlam <- 5
  lam1s <- rev(exp(seq(log(1.5), log(1.5 * 0.001), length.out = nlam)))
  lam2s <- rev(exp(seq(log(1.5), log(1.5 * 0.001), length.out = nlam)))

  n <- 10
  p <- 7
  x <- matrix(rnorm(n * p), nrow = n)
  y <- x %*% 1:p + rnorm(n)
  scout_res <- scout1something(x, y,
    p2 = 1,
    lam1s = lam1s, lam2s = lam2s,
    rescale = TRUE, trace = FALSE,
    intercept = FALSE, alternateCov = "default"
  )
  expect_equal(dim(scout_res$betamat), c(nlam, nlam, p))
  expect_equal(dim(scout_res$g.out$w), c(p, p))
  expect_false(is.null(scout_res$g.out$wi))
  expect_equal(dim(scout_res$g.out$wi), c(p, p))
})

test_that("scout1something_huge runs for multiple lam1s and lam2s", {
  nlam <- 5
  lam1s <- rev(exp(seq(log(1.5), log(1.5 * 0.001), length.out = nlam)))
  lam2s <- rev(exp(seq(log(1.5), log(1.5 * 0.001), length.out = nlam)))

  n <- 10
  p <- 7
  x <- matrix(rnorm(n * p), nrow = n)
  y <- x %*% 1:p + rnorm(n)
  scout_res <- scout1something_huge(x, y,
    p2 = 1,
    lam1s = lam1s, lam2s = lam2s,
    rescale = TRUE, trace = FALSE,
    intercept = FALSE, alternateCov = "default"
  )
  expect_equal(dim(scout_res$betamat), c(nlam, nlam, p))
  expect_equal(dim(scout_res$g.out$w), c(p, p))
  expect_false(is.null(scout_res$g.out$wi))
  expect_equal(dim(scout_res$g.out$wi), c(p, p))
})
