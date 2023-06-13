make_test_xy_data <- function(n = 10, p = 10) {
  x <- MASS::mvrnorm(n, rep(0, p), diag(p) * 0.5)
  x_test <- MASS::mvrnorm(n, rep(0, p), diag(p) * 0.5)
  beta <- rep(1, p)
  y <- x %*% beta + rnorm(n, 0, 0.1)
  list(x = x, y = y, x_test = x_test)
}

test_that("use_intercept argument of predict works as expected", {
  data <- make_test_xy_data()
  model <- scout(data$x, data$y, p1 = 1, p2 = 1, lam1s = 0.1, lam2s = 0.1, trace = F)
  # use_intercept is TRUE by default
  preds_with_icp <- as.vector(predict(model, data$x_test))
  preds_no_icp <- as.vector(predict(model, data$x_test, use_intercept = FALSE))
  expect_equal(preds_with_icp, preds_no_icp + as.vector(model$intercepts))
})

test_that("predict works if newx is a vector instead of matrix", {
  data <- make_test_xy_data(p = 1)
  model <- scout(data$x, data$y, p1 = 1, p2 = 1, lam1s = 0.1, lam2s = 0.1, trace = F)
  expect_error(predict(model, as.vector(data$x_test)), NA)
})

test_that("predict errors out if newx does not have the right dimensions", {
  data <- make_test_xy_data(p = 10)
  test_data_12 <- make_test_xy_data(p = 12)
  model <- scout(data$x, data$y, p1 = 1, p2 = 1, lam1s = 0.1, lam2s = 0.1, trace = F)
  # Too many features in test X matrix
  expect_error(
    predict(model, as.vector(test_data_12$x_test)),
    "newx does not have the correct number of features"
  )
  # Test X matrix is a vector but model wasn't fitted with 1 feature
  expect_error(
    predict(model, 1:nrow(data$x)),
    "newx does not have the correct number of features"
  )
})
