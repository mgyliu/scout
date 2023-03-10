mean_cw <- function(x) {
  locScale.x <- cellWise::estLocScale(x)
  locScale.x$loc
}

scale_cw <- function(x) {
  locScale.x <- cellWise::estLocScale(x)
  locScale.x$scale
}

mycov <- function(x, y = NULL, alternateCov = NULL) {
  n = NROW(x)
  p = NCOL(x)
  q = NCOL(y)

  if (is.null(alternateCov)) {
    return(stats::cov(x, y))
  }

  # At this point we know alternateCov is not null
  stopifnot("alternateCov must be a character type" = class(alternateCov) == "character")

  options = c("mcd", "cellwise", "mve")
  stopifnot(alternateCov %in% options)

  if (alternateCov == "cellwise") {
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

  if (alternateCov == "mcd") {    
    if (is.null(y)) {
      cov_mcd <- robustbase::covMcd(x)$cov  # Cov(X,X)
    } else {
      cov_full <- robustbase::covMcd(cbind(x,y))$cov 
      # cov_full is (p+q) x (p+q)
      # Cov(X,Y) is the top-right matrix of cov_full here
      cov_mcd <- as.matrix(cov_full[1:p, (p+1):(p+q)])
    }
    return(cov_mcd)
  }

  if (alternateCov == "mve") {
    if (is.null(y)) {
      cov_mve <- rrcov::CovMve(x)$cov
    } else {
      cov_full <- rrcov::CovMve(cbind(x,y))$cov
      # cov_full is (p+q) x (p+q)
      # Cov(X,Y) is the top-right matrix of cov_full here
      cov_mve <- as.matrix(cov_full[1:p, (p+1):(p+q)])
    }
    return(cov_mve)
  }

  stop(paste("alternateCov = ", alternateCov, " not yet implemented!", sep = ""))
}

# Returns list of x, y, meanx, meany, sdx, sdy
mystandardize <- function(x, y, standardize = TRUE, alternateCov = NULL) {
  # Default
  meanx <- apply(x, 2, mean)
  meany <- mean(y)
  
  sdx <- apply(x, 2, sd)
  sdy <- sd(y)

  # Use robust centering and scaling  
  if (!is.null(alternateCov)) {
    meanx <- apply(x, 2, mean_cw)
    meany <- mean_cw(y)

    sdx <- apply(x, 2, scale_cw)
    sdy <- scale_cw(y)
  }

  if (standardize) {
    x <- scale(x, center = meanx, scale = sdx)
  } else {
    x <- scale(x, center = meanx, F)
    sdx <- rep(1, ncol(x))
    sdy <- 1
  }
  # Scale and center the Y regardless
  # Original standardize code is here: https://github.com/cran/scout/blob/master/R/scout.R#L136-L147
  y <- (y - meany)/sdy

  return(list(x=x, y=y, meanx=meanx, meany=meany, sdx=sdx, sdy=sdy))
}

# yhat: array of dim (n_lam1 x n_lam2 x n_obs)
# ytrue: vector of length n_obs
# returns: (nlam1 x nlam2) matrix of computed metrics
#          corresponds to the metric for each lam1, lam2 combination
compute_cv_metric <- function(yhat, ytrue, metric = "mse") {
  lengths_match <- all(apply(yhat, c(1,2), length) == length(ytrue))
  stopifnot("yhat and ytrue dimension mismatch" = lengths_match)

  # For each (lam1, lam2), compute yhat - ytrue
  # resids is (n_lam1 x n_lam2 x n_omit)
  resids <- sweep(yhat, MARGIN = 3, STATS = ytrue, FUN = "-")
  
  if (metric == "mse") {
    sq_resids <- resids^2
    # For each (lam1, lam2), compute the MSE
    return(apply(sq_resids, c(1,2), mean))
  } 

  else if (metric == "tau_size") {
    return(apply(resids, c(1,2), pense::tau_size))
  }

  # median absolute prediction error
  else if (metric == "median_ape") {
    return(apply(resids, c(1,2), function(r) median(abs(r))))
  }

  # mean absolute prediction error
  else if (metric == "mean_ape") {
    return(apply(resids, c(1,2), function(r) mean(abs(r))))
  }

  stop(paste("cv metric", metric, "not yet implemented"))
}

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


scout1something <- function(x, y, p2, lam1s, lam2s, rescale,trace, intercept = FALSE, alternateCov = NULL){
  if (ncol(x)>500) {
    print("You are running scout with p1=1 and ncol(x) > 500. This will be slow. You may want to re-start and use p1=2, which is much faster.")
  }

  if (min(lam1s)==0 && min(lam2s)==0 && ncol(x)>=nrow(x)) {
    stop("don't run w/lam1=0 and lam2=0 when p>=n")
  }

  if (sum(order(lam2s)==(1:length(lam2s)))!=length(lam2s)) {
    stop("Error!!!! lam2s must be ordered!!!")
  }

  # Init a 3D matrix of NAs
  betamat <- array(NA, dim=c(length(lam1s), length(lam2s), ncol(x)))

  if (sum(lam1s>0 & lam1s<1e-4)>0 && ncol(x)>=nrow(x)){
    warning("Non-zero lam1s that were smaller than 1e-4 were increased to 1e-4 to avoid problems with graphical lasso that can occur when p>n.")
    lam1s[lam1s < 1e-4 & lam1s!=0] <- 1e-4
  }
  
  # Iterate through lambda 1's
  for (i in 1:length(lam1s)) {
    if(trace) cat(i,fill=F)
    g.out <- NULL

    # On first lambda_1
    if (i==1 || is.null(g.out$w) || is.null(g.out$wi)) {
      if (lam1s[i]!=0) {
        # Estimates a sparse inverse covariance matrix using a lasso (L1) penalty
        g.out <- glasso::glasso(mycov(x, alternateCov = alternateCov), rho=lam1s[i])
      }
      if (lam1s[i]==0) {
        g.out <- list(w = mycov(x, alternateCov = alternateCov), wi = NULL)
      }
    } else if(i!=1 && !is.null(g.out$w) && !is.null(g.out$wi)) {
      g.out <- glasso::glasso(mycov(x, alternateCov = alternateCov), rho=lam1s[i], start="warm", w.init=g.out$w, wi.init=g.out$wi)
    }

    # g.out values
    # * w: estimated covariance matrix
    # * wi: estimated inverse covariance matrix

    sigma_xx_hat <- g.out$w 
    theta_xx_hat <- g.out$wi

    # Iterate through lambda 2's
    # Step 2 in the algorithm
    for(j in 1:length(lam2s)){
      if (p2 == 0 || lam2s[j] == 0) {
        beta <- g.out$wi %*% mycov(x,y, alternateCov = alternateCov)
      } else if (p2 == 1 && lam2s[j] != 0) {
        if (j==1) {
          l_one_res <- lasso_one(g.out$w, mycov(x,y, alternateCov = alternateCov), rho=lam2s[j])
          beta <- l_one_res$beta
        }
        if (j!=1) {
          if(sum(abs(beta))!= 0 || lam2s[j] < lam2s[j-1]) { 
            l_one_res <- lasso_one(g.out$w, mycov(x,y, alternateCov = alternateCov), rho=lam2s[j], beta.init=beta)
            beta <- l_one_res$beta
            # if got zero for smaller value of lambda 2,
            # then no need to keep computing!!!
          }
        }  
      }
      if(rescale && sum(abs(beta))!=0) {
        # Step 4: rescale beta_hat* = c * beta_hat
        beta <- beta*lsfit(x%*%beta,y,intercept=intercept)$coef
      }
      betamat[i,j,] <- beta
    }
  }
  return(betamat)
}


scout2something <- function(x, y, p2, lam1s, lam2s,rescale, trace, intercept = FALSE, alternateCov = NULL){
  if(sum(order(lam2s)==(1:length(lam2s)))!=length(lam2s)){
    stop("Error!!!! lam2s must be ordered!!!")
  }
  if(min(lam1s)==0 && min(lam2s)==0 && ncol(x)>=nrow(x)) stop("don't run w/lam1=0 and lam2=0 when p>=n")
  g.out <- NULL
  betamat <- array(NA, dim=c(length(lam1s), length(lam2s), ncol(x)))
  for(i in 1:length(lam1s)){
    if (trace) cat(i,fill=F)
    if (lam1s[i]!=0) {
      if (i==1 || is.null(g.out)) {
        g.out <- gridge(x, rho=lam1s[i])
      }
      if (i!=1 && !is.null(g.out)) {
        g.out <- gridge(x, rho=lam1s[i], v=g.out$svdstuff$v, thetas=g.out$svdstuff$thetas, u=g.out$svdstuff$u)
      }
      for(j in 1:length(lam2s)){
        if (p2==0){
          beta <- diag(rep(g.out$wistuff$firstdiag, ncol(x))) %*% mycov(x,y, alternateCov = alternateCov) + 
            g.out$wistuff$v %*% (diag(g.out$wistuff$diagsandwich) %*% ((t(g.out$wistuff$v)) %*% mycov(x,y, alternateCov = alternateCov)))
        } else if (p2!=0 && p2==1) {
          if (j == 1) {
            l_one_res <- lasso_one(
              diag(rep(g.out$wstuff$firstdiag, ncol(x))) + g.out$wstuff$v %*% diag(g.out$wstuff$diagsandwich) %*% t(g.out$wstuff$v), 
              mycov(x,y, alternateCov = alternateCov), 
              rho=lam2s[j]
            )
            beta <- l_one_res$beta 
          }
          if (j != 1) {
            if (sum(abs(beta))!=0 || lam2s[j]<lam2s[j-1]) {
              l_one_res <- lasso_one(
                diag(rep(g.out$wstuff$firstdiag, ncol(x))) + g.out$wstuff$v %*% diag(g.out$wstuff$diagsandwich) %*% t(g.out$wstuff$v), 
                mycov(x,y, alternateCov = alternateCov), 
                rho=lam2s[j], 
                beta.init=beta
              )
              beta <- l_one_res$beta
              # If you got zero for a smaller value of
              # lambda2, then no need to keep computing!!!!!!!
            }
          }  
        }
        if(rescale && sum(abs(beta))!=0) {
          beta <- beta*lsfit(x%*%beta,y,intercept=intercept)$coef
        }
        betamat[i,j,] <- beta
      }
    } else if(lam1s[i]==0){
      if(p2==0) betamat[i,1,] <- lsfit(x,y,intercept=intercept)$coef
      if(p2==1){
        for(j in 1:length(lam2s)){
          if (lam2s[j]==0) { 
            beta <- lsfit(x,y,intercept=intercept)$coef
          }
          if (lam2s[j]!=0) { 
            l_one_res <- lasso_one(mycov(x, alternateCov = alternateCov),mycov(x,y, alternateCov = alternateCov), rho=lam2s[j])
            beta <- l_one_res$beta
          }
          if(sum(abs(beta))!=0 && rescale){
            betamat[i,j,] <- beta*lsfit(x%*%beta,y,intercept=intercept)$coef
          } else {
            betamat[i,j,] <- beta
          }
        }
      }
    }
  }
  return(betamat)
}


predict.scoutobject <- function(object, newx, ...){
  if(class(object)!="scoutobject") stop("Object must be of class 'scoutobject', created by a call to 'scout' function.")
   scout.obj <- object
   if(!is.matrix(newx) && length(newx)!=length(scout.obj$coef[1,1,])) stop("newx should be a single observation in vector form, or a matrix of observations. If in matrix form, then there should be one observation per row, with the features on the columns, just like the x matrix.")
   lam1s <- scout.obj$lam1s
   lam2s <- scout.obj$lam2s
   interceptmat <- scout.obj$intercepts
   betamat <- scout.obj$coefficients
   if(is.matrix(newx)) yhat <- array(dim=c(length(lam1s),length(lam2s),nrow(newx)))
   if(!is.matrix(newx)) yhat <- matrix(0,nrow=length(lam1s),ncol=length(lam2s))
   for(i in 1:length(lam1s)){
     for(j in 1:length(lam2s)){
       if(is.matrix(newx)) yhat[i,j,] <- interceptmat[i,j] + newx%*%matrix(betamat[i,j,],ncol=1)
       if(!is.matrix(newx)) yhat[i,j] <- interceptmat[i,j] + matrix(newx,nrow=1)%*%matrix(betamat[i,j,],ncol=1)
     }
   }
   return(yhat)
}

scout <- function(
  x, y, newx = NULL, p1 = 2, p2 = 1, 
  lam1s = seq(.001, .2, len=10), lam2s = seq(.001, .2, len=10), 
  rescale = TRUE, trace = TRUE, standardize = TRUE, 
  intercept = FALSE, alternateCov = NULL) {
  call <- match.call()
  if(!is.null(p1) && p1!=1 && p1!=2) stop("p1 must be 1, 2, or NULL.")
  if(!is.null(p2) && p2!=1) stop("p1 must be 1 or NULL.")
  if((sum(is.na(x)) + sum(is.na(y)))>0) stop("Please fix the NAs in your data set first. Missing values can be imputed using library 'impute'.")
  x <- as.matrix(x)
  if(min(apply(x,2,sd))==0) stop("Please do not enter an x matrix with variables that are constant.")

  # do standardization if needed
  std_result <- mystandardize(x, y, standardize, alternateCov)
  meanx <- std_result$meanx
  meany <- std_result$meany
  sdx <- std_result$sdx
  sdy <- std_result$sdy
  x <- std_result$x 
  y <- std_result$y

  # Want to re-order lam2s in increasing order
  if(length(lam2s)>0){
    lam2s.orig <- lam2s
    s2.out <- sort(lam2s.orig,index.return=TRUE)
    lam2s <- lam2s.orig[s2.out$ix]
  }
  # Done re-ordering lam2s
  # Want to re-order lam1s in increasing order
  if(length(lam1s)>0){
    lam1s.orig <- lam1s
    s1.out <- sort(lam1s.orig,index.return=TRUE)
    lam1s <- lam1s.orig[s1.out$ix]
  }
  # Done re-ordering lam1s
  if(is.null(lam1s) || is.null(p1) || sum(lam1s)==0){ p1 <- 0; lam1s=c(0); lam1s.orig=c(0) }
  if(is.null(lam2s) || is.null(p2) || sum(lam2s)==0){ p2 <- 0; lam2s=c(0); lam2s.orig=c(0) }
  if(p1!=0 && p1!=1 && p1!=2) stop("p1 must be 1, 2, or 0")
  if(p2!=0 && p2!=1) stop("p2 must be 1 or 0")
  if(!is.null(lam1s) && min(lam1s)<0) stop("lam1s cannot be negative")
  if(!is.null(lam2s) && min(lam2s)<0) stop("lam2s cannot be negative")
  if(ncol(x) >= nrow(x) && p1==0 && p2==0){
      stop("p1 and p2 cannot both be zero when p>n.")
  }

  # Conditions on p1 start here
  if(p1==0){
    betamat <- array(NA,dim=c(1,length(lam2s),ncol(x)))
    for(j in 1:length(lam2s)){
      if(trace) cat(j,fill=F)
      if(lam2s[j]==0){
        if(ncol(x)>=nrow(x)) stop("Cannot do Least Squares when ncol(x)>=nrow(x)")
        beta <- lsfit(x,y,intercept=intercept)$coef
        betamat[1,j,] <- beta
      } else {
        if(j==1) {
          l_one_res <- lasso_one(mycov(x, alternateCov = alternateCov), mycov(x,y, alternateCov = alternateCov), rho=lam2s[j])
          beta <- l_one_res$beta
          the_cov <- l_one_res$the_cov
        }
        if(j!=1) {
          if(sum(abs(beta))!=0 || lam2s[j]<lam2s[j-1]){ 
            l_one_res <- lasso_one(mycov(x, alternateCov = alternateCov), mycov(x,y, alternateCov = alternateCov), rho=lam2s[j], beta.init=beta)
            beta <- l_one_res$beta
            the_cov <- l_one_res$the_cov
            # if got zero for smaller value of lambda 2,
            # then no need to keep computing!!!
          }
        }
        if(rescale && sum(abs(beta))!=0) {
          beta <- beta*lsfit(x%*%beta,y,intercept=intercept)$coef
        }
        betamat[1,j,] <- beta
      }
    }  
  } else if(p1==1){
    betamat <- scout1something(x, y, p2, lam1s, lam2s, rescale, trace, alternateCov = alternateCov)
  } else if (p1==2){
    betamat <- scout2something(x, y, p2, lam1s, lam2s, rescale, trace, alternateCov = alternateCov)
  } 
  interceptmat <- matrix(meany,nrow=length(lam1s),ncol=length(lam2s))
  for(i in 1:length(lam1s)){
    for(j in 1:length(lam2s)){
     interceptmat[i,j] <- interceptmat[i,j] -  sum((sdy*meanx/sdx)*betamat[i,j,])
    }
  }
  betamat <- sweep(betamat, 3, sdy/sdx, "*")
  betamat <- array(
   betamat[rank(lam1s.orig),rank(lam2s.orig),],
   dim = c(length(lam1s.orig), length(lam2s.orig), ncol(x))
  )
  interceptmat <- matrix(interceptmat[rank(lam1s.orig),rank(lam2s.orig)],nrow=length(lam1s.orig),ncol=length(lam2s.orig))
  scout.obj <- (list(intercepts=interceptmat,coefficients=betamat,p1=p1,p2=p2,lam1s=lam1s.orig,lam2s=lam2s.orig, yhat=NULL,call=call))
  class(scout.obj) <- "scoutobject"
  if(!is.null(newx)){
    yhat <- predict.scoutobject(scout.obj,newx)
    scout.obj$yhat <- yhat
  }
  return(scout.obj)
}

print.scoutobject <- function(x,...){
  if(class(x)!="scoutobject") stop("Class of x must be 'scoutobject', created by call to function 'scout'.")
  scout.obj <- x
  cat("Call:\n",fill=F)
  dput(scout.obj$call)
  p1 <- scout.obj$p1
  p2 <- scout.obj$p2
  if(p1==0) p1 <- "NULL"
  if(p2==0) p2 <- "NULL"
  cat("\n Scout(",p1,",",p2,") was performed with lambda1 = (",scout.obj$lam1s,") and lambda2 = (", scout.obj$lam2s," ).\n")
  cat("\n Number of non-zero coefficients for each (lambda1, lambda2) pair:\n")
  mat <- apply(x$coef!=0,c(1,2),sum)
  dimnames(mat) <- list(round(x$lam1s,3),round(x$lam2s,3))
  print(mat,quote=F)
  invisible()
}   


cv.folds <- function(n, folds = 10){
  split(sample(1:n), rep(1:folds, length = n))
}

cv.scout <- function(
  x, y, 
  K = 10, 
  lam1s = seq(0.001, 0.2, len=10), 
  lam2s = seq(0.001, 0.2, len=10), 
  p1 = 2, p2 = 1,
  trace = TRUE, plot=TRUE, plotSE = FALSE, 
  rescale = TRUE, intercept = FALSE, alternateCov = NULL, ...) {
  call <- match.call()
  if(K==1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
  if(K > length(y)/2) stop("Please choose a value of K between 2 and length(y)/2.")
  if(p1==0 && p2==0) stop("Why would you want to cross-validate least squares?")
  if(is.null(p1)) lam1s <- 0
  if(is.null(p2)) lam2s <- 0
  lam1s <- c(lam1s)
  lam2s <- c(lam2s)
  if(length(lam1s)<2 && length(lam2s)<2) stop("Not a reasonable range of lambdas over which to be cross-validating")
  all.folds <- cv.folds(length(y), K)

  if (length(lam1s)>1 && length(lam2s)>1) {
    residmat <- array(0, dim=c(length(lam1s), length(lam2s), K))

    for (i in seq(K)) {
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <- scout(x[ - omit,  ], y[ - omit], newx=x[omit,], p1=p1, p2=p2, lam1s=lam1s,  lam2s=lam2s,rescale=rescale, trace=trace, intercept = intercept, alternateCov = alternateCov)
      # For each (lam1, lam2), compute yhat - ytrue
      # resids is (n_lam1 x n_lam2 x n_omit)
      resids <- sweep(fit$yhat, MARGIN = 3, STATS = y[omit], FUN = "-")
      sq_resids <- resids^2
      # For each (lam1, lam2), compute the MSE
      residmat[,,i] <- apply(sq_resids, c(1,2), mean)
    }

    # Now residmat is a (n_lam1 x n_lam2 x K) matrix 
    # For each lam1, lam2, and fold, residmat stores the MSPE from evaluating
    # the model on that fold
    # 
    # Next line averages out the MSPE over the folds
    cv <- apply(residmat, c(1,2), mean)
    # TODO: cv.error using `var` which is not robust.
    # use robust measures inside CV here
    cv.error <- sqrt(apply(residmat, c(1,2), var)/K)
    object<-list(p1=p1,p2=p2,lam1s=lam1s,lam2s=lam2s, cv = cv, cv.error = cv.error, call=call, bestlam1=lam1s[which.min(apply(cv,1,min))], bestlam2=lam2s[which.min(apply(cv,2,min))])
    if(plot){
      myrainbow <- rainbow(length(lam1s)*1.15)
      if(!plotSE) plot(0,0,xlim=range(lam2s),ylim=range(cv),col="white",main="CV Error", xlab="lam2", ylab="CV Error")
      if(plotSE) plot(0,0,xlim=range(lam2s),ylim=range(c(cv, cv+cv.error, cv-cv.error)),col="white",main="CV Error", xlab="lam2", ylab="CV Error")
      for(i in 1:length(lam1s)){
        lines(lam2s, cv[i,],col=myrainbow[i])
        points(lam2s, cv[i,], col=myrainbow[i])
        if(plotSE){
          points(lam2s,cv[i,]+cv.error[i,],col=myrainbow[i],type="l",lty="dashed")
          points(lam2s,cv[i,]-cv.error[i,],col=myrainbow[i],type="l",lty="dashed")
        }
      }
      legend("topright",pch=15,col=myrainbow[1:length(lam1s)],paste("lam1=",sep="",round(lam1s,4)))
    }
  } else if(length(lam1s)==1){
    lam1 <- lam1s[1]
    residmat <- matrix(0,nrow=length(lam2s),ncol=K)
    for(i in seq(K)){
      if(trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <-  scout(x[ - omit,  ], y[ - omit], newx=x[omit,], p1=p1, p2=p2, lam1s=lam1, lam2s=lam2s, rescale=rescale, trace = trace, intercept = intercept, alternateCov = alternateCov)
      residmat[,i] <- apply(sweep(fit$yhat[1,,],2,y[omit],"-")^2,1,mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat,1,var)/K) 
    object<-list(p1=p1,p2=p2,lam1s=lam1,lam2s=lam2s, cv = cv, cv.error = cv.error,call=call, bestlam1=lam1, bestlam2=lam2s[which.min(cv)])
    if(plot){
      if(!plotSE) plot(0,0,xlim=range(lam2s),ylim=range(cv),col="white",main="CV Error",xlab="lam2",ylab="CV Error")
      if(plotSE) plot(0,0,xlim=range(lam2s),ylim=range(c(cv, cv+cv.error, cv-cv.error)),col="white",main="CV Error", xlab="lam2", ylab="CV Error")
      lines(lam2s, cv, col="red")
      points(lam2s,cv,col="red")
      legend("topright",pch=15,col="red", paste("lam1=",sep="",lam1))
      if(plotSE){
        points(lam2s,cv+cv.error,col="red",type="l",lty="dashed")
        points(lam2s,cv-cv.error,col="red",type="l",lty="dashed")
      }
    }
  } else if (length(lam2s)==1){
    lam2 <- lam2s[1]
    residmat <- matrix(0,nrow=length(lam1s),ncol=K)
    for(i in seq(K)){
      if (trace) cat("\n CV Fold", i, "\t")
      omit <- all.folds[[i]]
      fit <-  scout(x[ - omit,  ], y[ - omit], newx=x[omit,], p1=p1,p2=p2,lam1s=lam1s,lam2s=lam2, trace = trace, rescle = rescale, intercept = intercept, alternateCov = alternateCov)
      residmat[,i] <- apply(sweep(fit$yhat[,1,],2,y[omit],"-")^2,1,mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat,1,var)/K) 
    object<-list(p1=p1,p2=p2,lam1s=lam1s,lam2s=lam2 , cv = cv, cv.error = cv.error, call=call, bestlam1=lam1s[which.min(cv)], bestlam2=lam2)
    if(plot){
      if(!plotSE) plot(0,0,xlim=range(lam1s),ylim=range(cv),col="white",main="CV Error",xlab="lam1",ylab="CV Error")
      if(plotSE) plot(0,0,xlim=range(lam1s),ylim=range(c(cv, cv+cv.error, cv-cv.error)),col="white",main="CV Error", xlab="lam1", ylab="CV Error")
      lines(lam1s, cv, col="red")
      points(lam1s,cv,col="red")
      if(plotSE){
        points(lam1s,cv+cv.error,col="red",type="l",lty="dashed")
        points(lam1s,cv-cv.error,col="red",type="l",lty="dashed")
      }
      legend("topright",pch=15,col="red",paste("lam2=",sep="",lam2))
    }
  }
  if(trace) cat("\n")
  object$call <- call
  object$folds <- all.folds
  class(object) <- "cvobject"
  invisible(object)
}

# Same as cv.scout but without plotting functions
# Shorter because I got rid of the cases depending on number of lambda1 and lambda2
# The result is exactly the same
cv2.scout <- function(
  x, y, 
  K = 10, 
  lam1s = seq(0.001, 0.2, len=10), 
  lam2s = seq(0.001, 0.2, len=10), 
  p1 = 2, p2 = 1,
  trace = TRUE,
  rescale = TRUE, 
  intercept = FALSE,
  alternateCov = NULL, 
  cvmetric = "mse",
  ...) {
  call <- match.call()
  if(K==1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
  if(K > length(y)/2) stop("Please choose a value of K between 2 and length(y)/2.")
  if(p1==0 && p2==0) stop("Why would you want to cross-validate least squares?")
  if(is.null(p1)) lam1s <- 0
  if(is.null(p2)) lam2s <- 0
  lam1s <- c(lam1s)
  lam2s <- c(lam2s)
  if(length(lam1s)<2 && length(lam2s)<2) stop("Not a reasonable range of lambdas over which to be cross-validating")
  all.folds <- cv.folds(length(y), K)

  residmat <- array(0, dim=c(length(lam1s), length(lam2s), K))

  for (i in seq(K)) {
    if (trace) cat("\n CV Fold", i, "\t")
    omit <- all.folds[[i]]
    fit <- scout(x[ - omit,  ], y[ - omit], newx=x[omit,], p1=p1, p2=p2, lam1s=lam1s,  lam2s=lam2s,rescale=rescale, trace=trace, intercept = intercept, alternateCov = alternateCov)

    residmat[,,i] <- compute_cv_metric(fit$yhat, y[omit], cvmetric)
  }

  # Now residmat is a (n_lam1 x n_lam2 x K) matrix 
  # For each lam1, lam2, and fold, residmat stores the MSPE from evaluating
  # the model on that fold
  # 
  # Next line averages out the MSPE over the folds
  cv <- apply(residmat, c(1,2), mean)

  cv.error <- sqrt(apply(residmat, c(1,2), var)/K)
  object<-list(
    p1 = p1,
    p2 = p2,
    lam1s = lam1s,
    lam2s = lam2s, 
    cv = cv, 
    cv.error = cv.error, 
    call = call,
    folds = all.folds,
    bestlam1 = lam1s[which.min(apply(cv,1,min))], 
    bestlam2 = lam2s[which.min(apply(cv,2,min))])

  if(trace) cat("\n")
  class(object) <- "cvobject"
  invisible(object)
}

print.cvobject <- function(x, ...){
  if(class(x)!="cvobject") stop("Class of x must be 'cvobject', created by call to cv.scout.")
  cat("Call:\t", fill=F)
  dput(x$call)
  cat("\n Cross-validation MSE for each lambda1/lambda2 pair: \n")
  mat <- matrix(round(x$cv,2),nrow=length(x$lam1s), ncol=length(x$lam2s))
  dimnames(mat) <- list(round(x$lam1s,3),round(x$lam2s,3))
  print(mat,quote=F)
  invisible()
}





