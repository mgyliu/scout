\name{print.scoutobject}
\alias{print.scoutobject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Print function for scout}
\description{
A function to print scout output
}
\usage{
\method{print}{scoutobject}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The results of a call to the "scout" function.}
  \item{...}{additional arguments; these are ignored.}
}
\references{Witten, DM and Tibshirani, R (2008) Covariance-regularized
  regression and classification for high-dimensional problems.  Journal
  of the Royal Statistical Society, Series B 71(3): 615-636. <http://www-stat.stanford.edu/~dwitten>}
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
xtest <- x2[-training,]
ytest <- y[-training]
# Done splitting data into training and test set
# Fit Model
scout.object <- scout(xtrain,ytrain,p1=1,p2=1,lam1s=c(0.001,0.1),lam2s=c(0.01,0.2))
print(scout.object)
# Done Fitting Model
detach(diabetes)
}
