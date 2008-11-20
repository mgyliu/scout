\name{cv.scout}
\alias{cv.scout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform cross-validation for covariance-regularized regression, aka the Scout.}
\description{This function returns cross-validation error rates for a range
  of lambda1 and lambda2 values, and also makes beautiful CV plots if
  plot=TRUE.
}
\usage{
cv.scout(x, y, K= 10, lam1s=seq(0.001,.2,len=10),lam2s=seq(0.001,.2,len=10),p1=2,p2=1,trace = TRUE, plot=TRUE,plotSE=FALSE,rescale=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\details{Pass in a data matrix x and a vector of outcomes y; it will perform
(10-fold) cross-validation over a range of lambda1 and lambda2
values. By default, Scout(2,1) is performed.
}
\arguments{
  \item{x}{A matrix of predictors, where the rows are the samples and
    the columns are the predictors}
  \item{y}{A matrix of observations, where length(y) should equal
nrow(x)}
  \item{K}{Number of cross-validation folds to be performed; default is
10}
  \item{lam1s}{The (vector of) tuning parameters for regularization of the
    covariance matrix. Can be NULL if p1=NULL, since then no covariance
    regularization is taking place. If p1=1 and nrow(x)<ncol(x), then the no value in lam1s
    should be smaller than 1e-3, because this will cause graphical lasso
    to take too long. Also, if ncol(x)>500 then we really do not
    recommend using p1=1, as graphical lasso can be uncomfortably slow.}
  \item{lam2s}{The (vector of) tuning parameters for the $L_1$ regularization of
    the regression coefficients, using the regularized covariacne
    matrix. Can be NULL if p2=NULL. (If p2=NULL, then non-zero lam2s
    have no effect). A value of 0 will result in no
    regularization. }
  \item{p1}{The $L_p$ penalty for the covariance regularization. Must be
  one of 1, 2, or NULL. NULL corresponds to no covariance
  regularization. }
  \item{p2}{The $L_p$ penalty for the estimation of the regression
    coefficients based on the regularized covariance matrix. Must be one
  of 1 (for $L_1$ regularization) or NULL (for no regularization).}
  \item{trace}{Print out progress as we go? Default is TRUE.}
  \item{plot}{If TRUE (by default), makes beautiful CV plots.}
  \item{plotSE}{Should those beautiful CV plots also display std error
    bars for the CV? Default is FALSE}
  \item{rescale}{Scout rescales coefficients, by default, in order to
    avoid over-shrinkage}
  \item{...}{Additional parameters}
}
\value{
  \item{folds}{The indices of the members of the K test sets are
    returned.}
  \item{cv}{A matrix of average cross-validation errors is returned.}
  \item{cv.error}{A matrix containing the standard errors of the
    elements in "cv", the matrix of average cross-validation errors.}
  \item{bestlam1}{Best value of lam1 found via cross-validation.}
  \item{bestlam2}{Best value fo lam2 found via cross-validation.}
  \item{lam1s}{Values of lam1 considered.}
  \item{lam2s}{Values of lam2 considered.}
}
\references{Witten, DM and Tibshirani, R (2008) Covariance-regularized
  regression and classification for high-dimensional problems. Journal
  of the Royal Statistical Society, Series B 71(3): 615-636. <http://www-stat.stanford.edu/~dwitten>}
\author{Daniela M. Witten and Robert Tibshirani}
\note{ 
}
\seealso{\link{scout}, \link{predict.scoutobject} }
\examples{ 
data(diabetes)
attach(diabetes)
par(mfrow=c(2,1))
par(mar=c(2,2,2,2))
cv.sc <- cv.scout(x2,y,p1=2,p2=1)
print(cv.sc)
cv.la <- cv.lars(x2,y)
print(c("Lars minimum CV is ", min(cv.la$cv)))
print(c("Scout(2,1) minimum CV is ", min(cv.sc$cv)))
detach(diabetes)
}
\keyword{ }
\keyword{ }
