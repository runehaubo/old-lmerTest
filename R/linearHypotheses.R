################################################################################       
## type 3 hypothesis SAS
################################################################################
makeContrastType3SAS <- function(model, term, L) {
  eps <- 1e-8
  #apply rule 1 (Goodnight 1976)
  
  #find all effects that contain term effect
  model.term <- terms(model)
  fac <- attr(model.term,"factors")
  names <- attr(model.term,"term.labels")
  classes.term <- attr(terms(model, FALSE), "dataClasses")
  
  cols.eff <- which(colnames(L)==term)
  num.relate <- relatives(classes.term, term, names, fac)
  if( length(num.relate)==0 )
    colnums <- setdiff(1:ncol(L), cols.eff)
  if( length(num.relate)>0 )
  {
    cols.contain <- NULL
    for( i in 1:length(num.relate) )
      cols.contain <- c(cols.contain, which(colnames(L)==names[num.relate[i]]))
    colnums <- setdiff(1:ncol(L), c(cols.eff, cols.contain))   
  }
  
  for(colnum in colnums)
  {
    pivots <- which(abs(L[, colnum]) > eps)
    if(length(pivots) > 0)
    {
      L[pivots[1], ] <- L[pivots[1], ] / L[pivots[1], colnum]
      nonzeros <- setdiff(pivots, pivots[1])
      if(length(nonzeros) != 0)
      {
        for( nonzero in nonzeros )
        {
          L[nonzero, ] <- L[nonzero, ]-L[nonzero, colnum]*L[pivots[1], ]
        }
      }
      L[pivots[1], ] <- rep(0, ncol(L))
    }
  }
  
  nums <- which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L <- L[nums,]
  
  if(is.vector(L))
    return(L)
  
  #orthogonalization
  if( length(cols.eff)>1 )
    zero.rows <- which(apply(L[,cols.eff],1,function(y) sum(abs(y)))==0)
  else
    zero.rows <- which(L[,cols.eff]==0)
  
  for(zero.row in zero.rows) 
  {
    w <- L[zero.row,]
    for(i in setdiff(1:nrow(L),zero.row))
    {
      if(sum(abs(L[i,]))!=0)
        L[i,] <- L[i,] - ((w %*% L[i,])/(w %*% w)) %*% w
    }
    L[zero.row,] <- rep(0, ncol(L))
  }
  
  L[abs(L)<1e-6] <- 0
  
  nums <- which(apply(L,1,function(y) sum(abs(y)))!=0) 
  L <- L[nums,]
  return(L)
}

################################################################################
###  type 1 SS
################################################################################
doolittle <- function(x, eps = 1e-6) {
  if(!is.matrix(x)) stop("argument 'x' is not a matrix")
  if(ncol(x) != nrow(x))
    stop( "argument x is not a square matrix" )
  if (!is.numeric(x) )
    stop( "argument x is not numeric" )
  n <- nrow(x)
  L <- U <- matrix(0, nrow=n, ncol=n)
  diag(L) <- rep(1, n)
  for(i in 1:n) {
    ip1 <- i + 1
    im1 <- i - 1
    for(j in 1:n) {
      U[i,j] <- x[i,j]
      if (im1 > 0) {
        for(k in 1:im1) {
          U[i,j] <- U[i,j] - L[i,k] * U[k,j]
        }
      }
    }
    if ( ip1 <= n ) {
      for ( j in ip1:n ) {
        L[j,i] <- x[j,i]
        if ( im1 > 0 ) {
          for ( k in 1:im1 ) {
            L[j,i] <- L[j,i] - L[j,k] * U[k,i]
          }
        }
        L[j, i] <- if(abs(U[i, i]) < eps) 0 else L[j,i] / U[i,i]
      }
    }
  }
  L[abs(L) < eps] <- 0
  U[abs(U) < eps] <- 0
  list( L=L, U=U )
}




## construct design matrix for F test 
createDesignMat <- function(rho) {
  model.term <- terms(rho$model)
  fixed.term <- attr(model.term,"term.labels") 
  X.design <- names.design.mat <-  names.design <- NULL
  X.design.red <- model.matrix(rho$model)
  attr(X.design.red, "dataClasses") <- 
    attr(terms(rho$model, FALSE), "dataClasses")
  dd <- model.frame(rho$model) 
  
  for(i in 1:length(fixed.term)) {
    formula.term <- as.formula(paste("~", fixed.term[i], "- 1"))
    X.design <- cbind(X.design, model.matrix(formula.term, dd))
    names.design.mat <- c(names.design.mat, 
                          rep(fixed.term[i],
                          ncol(model.matrix(formula.term, dd))))
  }
  
  if(attr(model.term, "intercept") != 0){
    names.design <- c("(Intercept)", colnames(X.design))
    X.design <- cbind(rep(1, dim(X.design)[1]), X.design)
    names.design.mat <- c("(Intercept)", names.design.mat)
  }
  else
    names.design <- colnames(X.design)
  colnames(X.design) <- names.design.mat
  
  if(length(which(colSums(X.design)==0)) != 0){
     warning(paste("missing cells for some factors (combinations of factors) \n", 
                   "care must be taken with type",
                   as.roman(3) ,
                   " hypothesis "))
  }
  
  fullCoefs <- rep(0, ncol(X.design))
  fullCoefs <- setNames(fullCoefs, names.design) 
  if("(Intercept)" %in% names.design)
    names(fullCoefs)[1] <- "(Intercept)"
  fullCoefs[names(rho$fixEffs)] <- rho$fixEffs
  nums.Coefs <- which(names(fullCoefs) %in% names(rho$fixEffs))
  nums.Coefs <- setNames(nums.Coefs, names(fullCoefs[nums.Coefs])) 
  Xlist <- list(X.design.red = X.design.red, 
                trms = model.term,
                X.design = X.design,
                names.design = names.design,
                fullCoefs = fullCoefs,
                nums.Coefs = nums.Coefs)
  return(Xlist)
}


## caclulate the General contrast matrix for the hypothesis (as in SAS)
calcGeneralSet<- function(rho, type) {
  if(type == 3)
    L <- calcGeneralSetType3(rho) #basis.set3(rho)#calcGeneralSetType3(rho)  
  else
  {
    L <- calcGeneralSet12(rho$Xlist$X.design.red)
    if(type == 1)
      L <- makeContrastType1(L, rho$Xlist$X.design.red, 
                             rho$Xlist$trms, rho$test.terms)
    if(type == 2)
      L <- makeContrastType2(L, rho$Xlist$X.design.red, 
                             rho$Xlist$trms, rho$test.terms)
  }
  L
}

      
## caclulate the General contrast matrix for the hypothesis (as in SAS)
calcGeneralSetType3 <- function(rho) {
  xtx <- t(rho$Xlist$X.design) %*% rho$Xlist$X.design
  g2 <- matrix(0, ncol=ncol(xtx), nrow=nrow(xtx))
  
  inds <- rho$Xlist$nums.Coefs
  g2[inds, inds] <- solve(xtx[inds, inds])
  g2[abs(g2) < 1e-10] <- 0
  
  #general set of estimable function
  L <- g2 %*% xtx
  L[abs(L) < 1e-6] <- 0
  L
}


## find which effect contains effect term
relatives <- function(classes.term, term, names, factors) {
  ## checks if the terms have the same number of covariates (if any)
  checkCovContain <- function(term1, term2) {        
    num.numeric <- which(classes.term=="numeric")
    num.numeric.term1 <- which((num.numeric %in% which(factors[,term1]!=0))==TRUE)
    num.numeric.term2 <- which((num.numeric %in% which(factors[,term2]!=0))==TRUE)
    if((length(num.numeric.term1)>0 && length(num.numeric.term2)>0)||
       (length(num.numeric.term1)==0 && length(num.numeric.term2)==0))
      return(all(num.numeric.term2 == num.numeric.term1))
    else
      return(FALSE)
  }
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2]))) && checkCovContain(term1,term2)
  }
  if(length(names) == 1) return(NULL)
  which.term <- which(term==names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
                                        function(term2) is.relative(term, term2))]
}


makeContrastType2 <- function(L, X, trms, trms.lab) {
  fac <- attr(trms,"factors")
  #trms.lab <- attr(trms,"term.labels")
  trms.assig <- attr(X, "assign")
  if(attr(trms,"intercept") == 0)
    attr(X, "assign") <- setNames(trms.assig, rep(trms.lab, table(trms.assig)))
  else
    attr(X, "assign") <- setNames(trms.assig, rep(c("(Intercept)", trms.lab), 
                                                  table(trms.assig)))
  Lt.list <- lapply(trms.lab, function(trm)
    makeContrastType2.term(trm, L, X, trms.lab, fac))
  names(Lt.list) <- trms.lab
  Lt.list
}


makeContrastType2.term <- function(trm, L, X, trms.lab, fac) {
  #find all effects that contain term effect
  num.relate <- relatives(attr(X, "dataClasses"), trm, trms.lab, fac)
  contain <- trms.lab[num.relate]
  if(length(contain) == 0 && (which(trms.lab == trm) == length(trms.lab))){
    Lc <- L[which(names(attr(X, "assign")) == trm), , drop = FALSE] 
    colnames(Lc) <- colnames(X)
    Lc
  } else {
    ## columns of the X are rearranged in a way that columns corresponding to the 
    ## effects that do not contain effet term are put before the columns corresponding 
    ## to the term
    ind.indep <- which(names(attr(X, "assign")) != trm & 
                         !(names(attr(X, "assign")) %in% contain))
    new.X <- cbind(X[,ind.indep, drop = FALSE], X[,-ind.indep, drop = FALSE])
    attr(new.X, "assign") <- c(attr(X, "assign")[ind.indep], attr(X, "assign")[-ind.indep])
    ## doolittle transform t(new.X) %*% new.X
    L <- calcGeneralSet12(new.X)
    colnames(L) <- colnames(new.X)
    ## columns of L are rearranged to reflect the original order
    Lc <- L[which(names(attr(new.X, "assign")) == trm), , drop = FALSE]
    Lc <- Lc[, colnames(X), drop = FALSE]
  }
  Lc
}


makeContrastType1 <- function(L, X, trms, trms.lab) {
  asgn <- attr(X, "assign")
  #trms.lab <- attr(trms, "term.labels")
  p <- ncol(X)
  ind.list <- split(1L:p, asgn)
  df <- unlist(lapply(ind.list, length))
  Lt.master <- L
  Lt.list <- lapply(ind.list, function(i) Lt.master[i, , drop=FALSE])
  if(attr(trms,"intercept") == 0)
    names(Lt.list) <- trms.lab
  else
    names(Lt.list) <- c("(Intercept)", trms.lab)
  Lt.list
}


calcGeneralSet12 <- function(X) {
  p <- ncol(X)
  XtX <- crossprod(X)
  U <- doolittle(XtX)$U
  d <- diag(U)
  for(i in 1:nrow(U))
    if(d[i] > 0) U[i, ] <- U[i, ] / d[i]
  L <- U
  L
}

