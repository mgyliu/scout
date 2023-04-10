
# x: data
# rho: lambda_1
# v: svdstuff$v 
# thetas: svdstuff$thetas [eigenvalues]
# u: svdstuff$u [eigenvectors]
gridge <- function(x, rho=0, v=NULL, thetas=NULL, u=NULL){
  x <- x/sqrt(nrow(x)-1) # stdize x
  p <- 2*rho # 2*lambda_1
  if(is.null(v)||is.null(thetas)||is.null(u)){
    covarswap <- x%*%t(x) # sample cov <- TODO replace here
    eigenswap <- eigen(covarswap)
    keep <- ((eigenswap$values) > 1e-11)
    v <- t(
      diag(1/sqrt(eigenswap$values[keep])) %*% 
      t(eigenswap$vectors[,keep]) %*% 
      x) # TODO: what is this
    thetas <- eigenswap$values[keep]
    u <- eigenswap$vectors[,keep]
  }
  lambda <- sqrt(4*p)/2
  diagmat <- (-lambda+ (-thetas + sqrt(thetas^2 + 4*p))/2)
  dbar <- .5*(thetas+sqrt(thetas^2+4*p))-sqrt(p)
  
  # diag prob getting penalized version of var-covar matrix

  return(list(
    svdstuff = list(u=u, v=v, thetas=thetas), 
    wistuff = list( # precision matrix
      v=v, 
      firstdiag=(1/sqrt(p)), 
      diagsandwich = (-((1/p)*(1/(1/sqrt(p) + 1/dbar))))
    ), 
    wstuff = list( # covariance matrix
      v=v, 
      firstdiag=lambda,
      diagsandwich=(thetas+diagmat))
    )
  )
}
