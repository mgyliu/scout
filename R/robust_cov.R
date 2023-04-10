cov_cellwise <- function(x, y = NULL) {
  locScale.x <- cellWise::estLocScale(x)
  # the wrapped data is stored in $Xw. covariances get computed on this.
  Xw.x <- cellWise::wrap(x, locScale.x$loc, locScale.x$scale)$Xw
  if (is.null(y)) {
    cov_cellwise <- cov(Xw.x)
  } else {
    locScale.y <- cellWise::estLocScale(y)
    Xw.y <- cellWise::wrap(y, locScale.y$loc, locScale.y$scale)$Xw
    cov_cellwise <- cov(Xw.x, Xw.y)
  }
  return(cov_cellwise)
}

cov_mcd <- function(x, y = NULL) {
  if (is.null(y)) {
    cov_mcd <- robustbase::covMcd(x)$cov # Cov(X,X)
  } else {
    cov_full <- robustbase::covMcd(cbind(x, y))$cov
    # cov_full is (p+q) x (p+q)
    # Cov(X,Y) is the top-right matrix of cov_full here
    cov_mcd <- as.matrix(cov_full[1:p, (p + 1):(p + q)])
  }
  return(cov_mcd)
}

cov_mve <- function(x, y = NULL) {
  if (is.null(y)) {
    cov_mve <- rrcov::CovMve(x)$cov
  } else {
    cov_full <- rrcov::CovMve(cbind(x, y))$cov
    # cov_full is (p+q) x (p+q)
    # Cov(X,Y) is the top-right matrix of cov_full here
    cov_mve <- as.matrix(cov_full[1:p, (p + 1):(p + q)])
  }
  return(cov_mve)
}

rob_cov <- function(x, y = NULL, alternateCov = "default") {
  stopifnot("alternateCov must be a character type" = class(alternateCov) == "character")

  n <- NROW(x)
  p <- NCOL(x)
  q <- NCOL(y)

  if (alternateCov == "default") {
    return(stats::cov(x, y))
  }

  if (alternateCov == "cellwise") {
    return(cov_cellwise(x, y))
  }

  if (alternateCov == "mcd") {
    return(cov_mcd(x, y))
  }

  if (alternateCov == "mve") {
    return(cov_mve(x, y))
  }

  stop(paste("alternateCov = ", alternateCov, " not yet implemented!", sep = ""))
}
