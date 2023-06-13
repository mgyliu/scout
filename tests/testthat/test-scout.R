test_that("scout runs for all p1 in 0,1,2 and p2 in 0,1", {
  n <- 10
  p <- 7
  x <- matrix(rnorm(n * p), nrow = n)
  y <- x %*% 1:p + rnorm(n)

  scout_res <- scout(x, y, p1 = NULL, p2 = NULL, trace = FALSE)
  expect_true(is.null(scout_res$gout))
  scout_res <- scout(x, y, p1 = NULL, p2 = 1, trace = FALSE)
  expect_true(is.null(scout_res$gout))

  scout_res <- scout(x, y, p1 = 1, p2 = NULL, trace = FALSE)
  expect_false(is.null(scout_res$gout))
  scout_res <- scout(x, y, p1 = 1, p2 = 1, trace = FALSE)
  expect_false(is.null(scout_res$gout))

  scout_res <- scout(x, y, p1 = 2, p2 = NULL, trace = FALSE)
  expect_false(is.null(scout_res$gout))
  scout_res <- scout(x, y, p1 = 2, p2 = 1, trace = FALSE)
  expect_false(is.null(scout_res$gout))
})

test_that("rescale_betas works as expected", {
  n <- 10
  p <- 7
  x <- matrix(rnorm(n * p), nrow = n)
  y <- x %*% 1:p + rnorm(n)

  std_result <- rob_standardize(x, y)
  sdx <- std_result$sdx
  sdy <- std_result$sdy

  res.rescale <- scout(x, y, p1 = 1, p2 = NULL, lam1s = 0.1, lam2s = 0.1, trace = FALSE, rescale_betas = TRUE)
  res.rescale.coef <- as.numeric(res.rescale$coefficients)
  res.norescale <- scout(x, y, p1 = 1, p2 = NULL, lam1s = 0.1, lam2s = 0.1, trace = FALSE, rescale_betas = FALSE)
  res.norescale.coef <- as.numeric(res.norescale$coefficients)

  expect_equal(res.rescale.coef, res.norescale.coef * sdy / sdx)
})
