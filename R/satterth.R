## calculates satterthwate's approximation
## L is a matrix or a vector
calcSatt <- function(rho, L, isF = TRUE, calcSS = TRUE){
  if(is.vector(L))
    result <- calcSatterth1DF(rho, L, isF)
  else
    result <- calcSatterthMultDF(rho, L)
  if(!calcSS) {
    result$ss <- NULL
    result$ms <- NULL
  }
  result
}

## calculates t-test with satterthwate's approximation
calcSatterth1DF <- function(rho, L, isF = TRUE){
  result <- matrix(0, nrow = 1, ncol = 4)
  colnames(result) <- c("df", "t value", "p-value", "sqrt.varcor")
  vss <- vcovLThetaL(rho$model)
  g <- mygrad(function(x)  vss(t(L), x), c(rho$thopt, rho$sigma))
  denom <- t(g) %*% rho$A %*% g
  varcor <- vss(t(L), c(rho$thopt, rho$sigma)) ## for the theta and sigma parameters
  ## df
  result[,1] <- 2*(varcor)^2/denom
  ## statistics
  result[,2] <- (L %*% rho$fixEffs)/sqrt(varcor) 
  result[,3] <- 2*(1 - pt(abs(result[,2]), df = result[,1]))
  result[,4] <- sqrt(varcor) 
  if(isF){
    F.stat <- result[,2]^2
    ss <- ms <- F.stat * rho$sigma^2
    return(list(ss = ss, ms = ms, denom = result[,1], Fstat = F.stat, 
                 pvalue = result[,3], ndf=1))
  }
  else
    return(result)
}

## calculates t test with KR approximation uses vcovAdj function
## from pbkrtest package
calcKR1DF <- function(rho, L){
  
  if (!requireNamespace("pbkrtest", quietly = TRUE)) 
    stop("pbkrtest package required for Kenward-Roger's approximations")
  result <- matrix(0, nrow = 1, ncol = 4)
  colnames(result) <- c("df", "t value", "p-value", "sqrt.varcor")
  Va <- pbkrtest::vcovAdj(rho$model)
  .ddf <- pbkrtest::get_ddf_Lb(rho$model, L)      
  b.hat <- rho$fixEffs
  Lb.hat <- sum(L * b.hat)
  Va.Lb.hat <- t(L) %*% Va %*% L
  t.stat <- as.numeric(Lb.hat / sqrt(Va.Lb.hat))
  p.value <- 2 * pt(abs(t.stat), df = .ddf, lower.tail = FALSE)
  result[,1] <- .ddf
  result[,2] <- t.stat
  result[,3] <- p.value
  result[,4] <- as.numeric(sqrt(Va.Lb.hat))
  result
}

## calculates F test with KR approximation uses KRmodcomp function 
## from pbkrtest package
calcKRMultDF <- function(rho, Lc){
  if (!requireNamespace("pbkrtest", quietly = TRUE)) 
    stop("pbkrtest package required for Kenward-Roger's approximations")
  
  if(is.vector(Lc))
    res.KR <- pbkrtest::KRmodcomp( rho$model, t(as.matrix(Lc)) )
  else
    res.KR <- pbkrtest::KRmodcomp( rho$model, Lc )
  
  ## calculate ms and ss
  ms <- res.KR$test[1,"stat"] * rho$sigma^2
  ss <- ms * res.KR$test[1,"ndf"]

  return(list(denom = res.KR$test[1,"ddf"], Fstat = res.KR$test[1,"stat"], 
              pvalue =  res.KR$test[1,"p.value"], ndf = res.KR$test[1,"ndf"], 
              ss = ss , ms = ms))
}


## calculates F statistics with Satterthwaite's approximation
## L is a matrix
calcSatterthMultDF <- function(rho, Lc) {
  # F statistics for tested term
  vcov.final <- as.matrix(vcov(rho$model))
  if(is.vector(Lc))
    C.theta.optim <- as.matrix(t(Lc) %*% vcov.final %*% Lc)    
  else
    C.theta.optim <- as.matrix(Lc %*% vcov.final %*% t(Lc))    
  
  invC.theta <- tryCatch({solve(C.theta.optim)}, error = function(e) { NULL })
  if(is.null(invC.theta))
    return(list(denom = 0, Fstat = NA, pvalue = NA, ndf=NA, ss = NA, ms = NA))
  
  q <- qr(C.theta.optim)$rank
  F.stat <- (t(Lc %*% rho$fixEffs) %*% invC.theta %*% (Lc %*% rho$fixEffs)) / q
  
  #df for F statistics for tested term
  svdec <- eigen(C.theta.optim) 
  
  PL <- t(svdec$vectors) %*% Lc

  vss2 <- vcovTheta(rho$model)
  theopt <- c(rho$thopt, rho$sigma)
  g <- mygrad(function(x)  vss2(x), theopt)
  
  if(class(g) == "numeric")
    mat.grad <- llply(1:length(theopt), function(x) matrix(g[x], 
                                                           ncol = ncol(vcov.final), 
                                                           nrow = nrow(vcov.final)))
  else
    mat.grad <- llply(1:length(theopt), function(x) matrix(g[, x], 
                                                           ncol = ncol(vcov.final), 
                                                           nrow = nrow(vcov.final)))
  
  nu.m.fun <- function(m){    
    den.nu <- unlist(llply(1:length(mat.grad), function(x) 
      as.matrix(t(PL[m,]) %*% mat.grad[[x]] %*% PL[m,])))   
    2*(svdec$values[m])^2/(t(den.nu) %*% rho$A %*% den.nu)
  }
  
  
  nu.m <- unlist(llply(1:length(svdec$values), .fun = nu.m.fun))
  
  nu.m[which(abs(2 - nu.m) < 1e-5)] <- 2.00001
  
  E <- sum( (nu.m/(nu.m-2)) * as.numeric(nu.m>2))
  nu.F <- 2 * E * as.numeric(E > q) / (E - q)
  
  pvalueF <- 1 - pf(F.stat,qr(Lc)$rank, nu.F)
  
  # calculate ss and ms
  if(is.na(F.stat))
    ms <- ss <- NA
  else{
    ms <- F.stat * rho$sigma^2
    ss <- ms * q
  }
  
  ## calculate ss from camp method proc glm
  ## ss <- getSS(Lc, rho$fixEffs ,ginv(rho$XtX)) 
  return( list(ss = ss, ms = ms, denom = nu.F, Fstat = F.stat, 
               pvalue = pvalueF, ndf=q))  
}


## calculates asymptotic variance covariance matrix of variance parameters based on theta
calcApvar <- function(rho){
  ## based on theta parameters and sigma
  dd <- devfunTheta(rho$model)
  h <- myhess(dd, c(rho$thopt, sigma = rho$sigma))  
  
  ch <- try(chol(h), silent=TRUE)
  if(inherits(ch, "try-error")) {
    message("Error in calculation of the Satterthwaite's approximation. The output of lme4 package is returned")
  }
  A <- 2*chol2inv(ch)
  
  eigval <- eigen(h, symmetric=TRUE, only.values=TRUE)$values
  isposA <- TRUE
  if(min(eigval) < sqrt(.Machine$double.eps)) ## tol ~ sqrt(.Machine$double.eps)
    
  isposA <- FALSE
  if(!isposA)
    print("Asymptotic covariance matrix A is not positive!")
  A
}


## devfun function as a function of optimal parameters
devfunTheta <- function(fm) 
{
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  nf <- length(fixef(fm)) 
  if (!isGLMM(fm)) 
    np <- np + 1L
  n <- nrow(fm@pp$V)
  
  ff <- updateModel(fm, devFunOnly = TRUE)
  reml <- getME(fm, "is_REML")
  
  envff <- environment(ff)
  
  if (isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")
      ff(thpars[-np])
      
      sigsq <- thpars[np]^2
      dev <- envff$pp$ldL2() + (envff$resp$wrss() + envff$pp$sqrL(1))/sigsq + n * 
        log(2 * pi * sigsq)      
      if(reml){
        p <- ncol(envff$pp$RX())
        dev <- dev + 2*determinant(envff$pp$RX())$modulus - p * log(2 * pi * sigsq)              
      }
      return(dev)     
    }
  }
  
  attr(ans, "thopt") <- fm@pp$theta
  class(ans) <- "devfunTheta"
  ans
}




## returns Lc %*% vcov as a function of theta parameters %*% t(Lc)
vcovLThetaL <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  nf <- length(fixef(fm))
  if (!isGLMM(fm)) 
    np <- np + 1L
  

  ff2 <- updateModel(fm, devFunOnly = TRUE) 

  envff2 <- environment(ff2)
  
  if (isLMM(fm)) {
    ans <- function(Lc, thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(Lc %*% as.matrix(vcov_out) %*% t(Lc)))        
    }
  } 
  class(ans) <- "vcovLThetaL"
  ans
}


## returns variance covariance of fixed effects as a function of theta parameters
vcovTheta <- function(fm)
{
  stopifnot(is(fm, "merMod"))
  
  np <- length(fm@pp$theta)
  if (!isGLMM(fm)) 
    np <- np + 1L
  
  
  ff2 <- updateModel(fm, devFunOnly = TRUE) 
  envff2 <- environment(ff2)
  
  if (isLMM(fm)) {
    ans <- function(thpars) {
      stopifnot(is.numeric(thpars), length(thpars) == np)
      
      sigma2 <- thpars[np]^2
      ff2(thpars[-np])
      
      #.Call("lmer_Deviance", pp$ptr(), resp$ptr(), thpars[-np], PACKAGE = "lme4")      
      vcov_out <- sigma2 * tcrossprod(envff2$pp$RXi()) 
      
      return(as.matrix(vcov_out))
      
    }
  } 
  class(ans) <- "vcovTheta"
  ans
}

