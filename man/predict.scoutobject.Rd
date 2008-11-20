\name{predict.scoutobject}
\alias{predict.scoutobject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction function for covariance-regularized regression, aka the Scout.}
\description{
A function to perform prediction, using an x matrix and the output of
the "scout" function.
}
\usage{
predict.scoutobject(object, newx, ...)
}
%- maybe also 'usage' for other objects documented here.
\details{
}
\arguments{
  \item{object}{The results of a call to the "scout" function. The
    coefficients that are part of this object will be used for
    making predictions.}
  \item{newx}{The new x at which predictions should be made. Can be a
    vector of length ncol(x), where x is the data on which scout.obj was
    created, or a matrix with ncol(x) columns.}
  \item{...}{Additional arguments to predict}
}
\value{
  \item{yhat}{If newx was a vector, then a  matrix will be returned,
    with dimension length(lam1s)xlength(lam2s) (where lam1s and lam2s
    are attributes of scout.obj). The (i,j) element of this matrix will
    correspond to tuning parameter values (lam1s[i], lam2s[j]). If newx
    is a matrix, then an array of dimension
    nrow(newx)xlength(lam1s)xlength(lam2s) will be returned.}
}
\references{Witten, DM and Tibshirani, R (2008) Covariance-regularized
  regression and classification for high-dimensional problems.  Journal
  of the Royal Statistical Society, Series B 71(3): 615-636. <http://www-stat.stanford.edu/~dwitten>}
\author{Daniela M. Witten and Robert Tibshirani}
\note{ 
}
\seealso{\link{scout}, \link{cv.scout} }
\examples{
data(diabetes)
attach(diabetes)
# Split data into training and test set
training <- sample(nrow(x2),floor(nrow(x2)/2))
xtrain <- x2[training,]
ytrain <- y[training]
xtest <- x2[-training,]
ytest <- y[-training]
# Done splitting data into training and test set
# Do cross-validation to determine best tuning parameter values for Scout(1,1)
cv.out <- cv.scout(xtrain,ytrain,p1=1,p2=1, lam1s=seq(0.001,.15,len=10))
print(cv.out)
# Done cross-validation
# Fit Model
scout.object <- scout(xtrain,ytrain,p1=1,p2=1,lam1s=cv.out$bestlam1,lam2s=cv.out$bestlam2)
print(scout.object)
# Done Fitting Model
# Predict on test data, and report MSE
yhats <- predict(scout.object,xtest)
print(mean((yhats-ytest)^2))
detach(diabetes)
}
\keyword{ }
\keyword{ }
