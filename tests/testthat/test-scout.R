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
