mean_cw <- function(x) {
  locScale.x <- cellWise::estLocScale(x)
  locScale.x$loc
}

scale_cw <- function(x) {
  locScale.x <- cellWise::estLocScale(x)
  locScale.x$scale
}

#' Returns list of x, y, meanx, meany, sdx, sdy
#' @param x A matrix, n x p
#' @param y A vector or n x 1 matrix
#' @param standardize TRUE or FALSE, indicating whether standardization should occur.
#' @param alternateCov A string indicating which alternate covariance estimator to use. Currently implemented options are "default" and "cellwise"
#' @returns list with x, y, meanx, meany, sdx, and sdy
rob_standardize <- function(x, y, standardize = TRUE, alternateCov = "default") {
  c(
    rob_standardize_matrix(x, standardize, alternateCov),
    rob_standardize_vector(y, standardize, alternateCov)
  )
}

rob_standardize_matrix <- function(x, standardize = TRUE, alternateCov = "default") {
  if (alternateCov == "default") {
    # Default
    meanx <- apply(x, 2, mean)
    sdx <- apply(x, 2, sd)
  } else if (alternateCov == "cellwise") {
    # Use robust centering and scaling
    meanx <- apply(x, 2, mean_cw)
    sdx <- apply(x, 2, scale_cw)
  } else {
    stop(paste("alternateCov = ", alternateCov, " not yet implemented in rob_standardize_matrix function", sep = ""))
  }

  if (standardize) {
    x <- scale(x, center = meanx, scale = sdx)
  } else {
    x <- scale(x, center = meanx, F)
    sdx <- rep(1, ncol(x))
  }

  list(x = x, meanx = meanx, sdx = sdx)
}

rob_standardize_vector <- function(y, standardize = TRUE, alternateCov = "default") {
  if (alternateCov == "default") {
    # Default
    meany <- mean(y)
    sdy <- sd(y)
  } else if (alternateCov == "cellwise") {
    # Use robust centering and scaling
    meany <- mean_cw(y)
    sdy <- scale_cw(y)
  } else {
    stop(paste("alternateCov = ", alternateCov, " not yet implemented in rob_standardize function", sep = ""))
  }

  if (standardize == FALSE) {
    sdy <- 1
  }
  # Scale and center the Y regardless
  # Original standardize code is here: https://github.com/cran/scout/blob/master/R/scout.R#L136-L147
  y <- (y - meany) / sdy

  list(y = y, meany = meany, sdy = sdy)
}
