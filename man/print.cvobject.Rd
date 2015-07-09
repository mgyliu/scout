\name{print.cvobject}
\alias{print.cvobject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print function for scout}
\description{
A function to print CV output for scout
}
\usage{
\method{print}{cvobject}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The results of a call to the "cv.scout" function.}
  \item{...}{Additional arguments; ignored.}
}
\references{Witten, DM and Tibshirani, R (2008) Covariance-regularized
  regression and classification for high-dimensional problems.  Journal
  of the Royal Statistical Society, Series B 71(3): 615-636.}
\author{Daniela M. Witten and Robert Tibshirani}
\seealso{\link{scout}, \link{cv.scout} }
\examples{
library(lars)
data(diabetes)
attach(diabetes)
# Split data into training and test set
training <- sample(nrow(x2),floor(nrow(x2)/2))
xtrain <- x2[training,]
ytrain <- y[training]
# Done splitting data into training and test set
# Do cross-validation to determine best tuning parameter values for Scout(1,1)
\dontrun{cv.out <- cv.scout(xtrain,ytrain,p1=1,p2=1, lam1s=seq(0.001,0.1), K=4)}
\dontrun{print(cv.out)}
# Done cross-validation
detach(diabetes)
}
