make_test_xy_data <- function(n = 10, p = 10) {
  x <- MASS::mvrnorm(n, rep(0, p), diag(p) * 0.5)
  beta <- rep(1, p)
  y <- x %*% beta + rnorm(n, 0, 0.1)
  list(x = x, y = y)
}

test_that("scout_glasso_stepwise works for p2 = NULL", {
  data <- make_test_xy_data()

  res.default <- scout_glasso_stepwise(data$x, data$y, p2 = NULL, nlambda1 = 5, nlambda2 = 5)
  expect_named(res.default, c("g.res", "cv.res", "mod", "X_wrap_res", "Y_wrap_res"))

  res.cw <- scout_glasso_stepwise(data$x, data$y, p2 = NULL, nlambda1 = 5, nlambda2 = 5, alternateCov = "cellwise")
  expect_named(res.cw, c("g.res", "cv.res", "mod", "X_wrap_res", "Y_wrap_res"))
})

test_that("scout_glasso_stepwise works for p2 = 1", {
  data <- make_test_xy_data()

  res.default <- scout_glasso_stepwise(data$x, data$y, p2 = 1, K = 3, nlambda1 = 5, nlambda2 = 5)
  expect_named(res.default, c("g.res", "cv.res", "mod", "X_wrap_res", "Y_wrap_res"))
})

test_that("scout_glasso_stepwise errors if p2 is not 1 or NULL", {
  data <- make_test_xy_data()
  expect_error(
    scout_glasso_stepwise(data$x, data$y, p2 = 2),
    "not implemented"
  )
})
