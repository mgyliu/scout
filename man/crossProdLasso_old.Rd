\name{crossProdLasso}
\alias{crossProdLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs the lasso on the cross product matrices X'X and X'y}
\description{Perform L1-regularized regression of y onto X using only
  the cross-product matrices X'X and X'y. In the case of
  covariance-regularized regression, this is useful if you would like to
  try out something other than L1 or L2 regularization of the inverse
  covariance matrix.
 
  Suppose you
  use your own method to regularize X'X. Then let Sigma denote your estimate of the population covariance matrix.  Now say you want to minimize
  beta' Sigma  beta - 2 beta' X'y + lambda ||beta||_1
in order to get the regression estimate beta, which maximizes the
second scout criterion when an L_1 penalty is used. You can do this by
calling crossProdLasso(Sigma, X'y,rho).

If you run crossProdLasso(X'X,X'y,rho) then it should give the same
result as lars(X,y)

Notice that the xtx that you pass into this function must be POSITIVE
SEMI DEFINITE (or positive definite) or the problem is not convex and
the algorithm will not converge.
}
\usage{
crossProdLasso(xtx,xty,rho,thr=1e-4,maxit=100,beta.init=NULL)
}
\details{
If your xtx is simply X'X for some X, and your xty is simple X'y
with some y, then the results will be the same as running lars on data
(X,y) for a single shrinkage parameter value.

Note that when you use the scout function with p2=1, the crossProdLasso
function is
  called internally to give the regression coefficients, after the
  regularized inverse covariance matrix is estimated. It is provided
  here in case it is useful to the user in other settings.
}
\arguments{
  \item{xtx}{A pxp matrix, which should be an estimate of a covariance
    matrix. This matrix must be POSITIVE
SEMI DEFINITE (or positive definite) or the problem is not convex and
the algorithm will not converge.}
  \item{xty}{A px1 vector, which is generally obtained via X'y.}
  \item{rho}{Must be non-negative; the regularization parameter you are
    using. }
  \item{thr}{Convergence threshold.}
  \item{maxit}{How many iterations to perform?}
  \item{beta.init}{If you're running this over a range of rho values,
    then set beta.init equal to the solution you got for a previous rho
    value. It will speed things up.}
}
\value{
  \item{beta}{A px1 vector with the regression coefficients.}
}
\references{Witten, DM and Tibshirani, R (2008) Covariance-regularized
  regression and classification for high-dimensional problems. Journal
  of the Royal Statistical Society, Series B 71(3): 615-636. <http://www-stat.stanford.edu/~dwitten>}
\author{FORTRAN code by Jerry Friedman.   R interface by Daniela M. Witten and Robert Tibshirani}
\note{The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman.
}
\examples{
set.seed(1)
#data(diabetes)
#attach(diabetes)
x2 <- matrix(rnorm(10*20),ncol=20)
y <- rnorm(10)
# First, let's do scout(2,1) the usual way).
scout.out <- scout(x2,y,p1=2,p2=1)
print(scout.out)



# Now, suppose I want to do develop a covariance-regularized regression
# method as in Section 3.2 of Witten and Tibshirani (2008). It will work
# like this:
# 1. Develop some positive definite estimate of Sigma
# 2. Find \beta by minimize \beta^T \Sigma \beta - 2 \beta^T X^T y +
# \lamda ||\beta||_1
# 3. Re-scale \beta.

# Step 1:
regcovx <- cov(x2)*(abs(cov(x2))>.005) + diag(ncol(x2))*.01

# Step 2:
betahat <-  crossProdLasso(regcovx, cov(x2,y), rho=.02)$beta
# Step 3:
betahat.sc <- betahat*lsfit(x2\%*\%betahat, y, intercept=FALSE)$coef
print(betahat.sc)

# Try a different value of rho:
betahat2 <- crossProdLasso(regcovx,cov(x2,y),rho=.04,beta.init=betahat)$beta
plot(betahat,betahat2, xlab="rho=.02",ylab="rho=.04")
#detach(diabetes)
}
