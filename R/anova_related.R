calcANOVA <- function(model, ddf = "Satterthwaite", type = 3, 
                      change.contr = TRUE) {
  if(!(type %in% c(1, 2, 3)))  
    stop('Parameter type is wrongly specified') ## check type of hypothesis
  
  rho <- list() ## environment containing info about model
  if(type == 3)
    change.contr <- TRUE
  else
    change.contr <- FALSE
  rho <- rhoInit(rho, model, change.contr) ## save lmer outcome in rho envir variable

  if(ddf == "Satterthwaite")
    rho$A <- calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma  

  rho$Xlist <- createDesignMat(rho) ## X design matrix for fixed effects
  
  ## define the terms that are to be tested
  rho$test.terms <- attr(terms(rho$model),
                     "term.labels")[unique(attr(rho$Xlist$X.design.red, 
                                                "assign"))] 
  anova.table <- initAnovaTable(rho$model, rho$test.terms, FALSE) ## initialize anova table

  ## calculate general set of hypothesis matrix 
  L <- calcGeneralSet(rho, type) 
  
  resultFpvalueSS <- llply(rho$test.terms, calcFpvalueMAIN, L = L, 
                           rho = rho, 
                           ddf = ddf, type = type)
  ## fill anova table
  anova.table <- fillAnovaTable(resultFpvalueSS,  anova.table)
  anova.table
}

## calculate t test for the summary function for lmerMod object
calcSummary <- function(model, ddf = "Satterthwaite") {
  # Make rho-list:
  rho <- list() ## vector containing info about model
  rho <- rhoInit(rho, model, FALSE) ## save lmer outcome in rho envir variable
  # Make contrast matrix:
  Lmat <- diag(rep(1, length(rho$fixEffs)))
  if(ddf == "Satterthwaite") {
    rho$A <- calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma  
    result <- do.call(rbind, lapply(1:nrow(Lmat), function(i) # always a matrix
      calcSatterth1DF(rho=rho, L=Lmat[i, , drop=TRUE], isF=FALSE)))
  } else { # KR-method:
    result <- do.call(rbind, lapply(1:nrow(Lmat), function(i) 
      calcKR1DF(rho=rho, L=Lmat[i, , drop=TRUE])))
  }
  return(list(df = result[, "df"], tvalue = result[, "t value"], 
              tpvalue = result[, "p-value"]))
}


################################################################################
## function to calculate F stat and pvalues for a given term
################################################################################
calcFpvalueSS <- function(Lc, rho, ddf, type) {
  if(is.null(Lc))
    return(NULL) 
  ## BUG: check vases example from Per
  # calculate ss
  # ss = 1##getSS(Lc, fullCoefs, ginv(crossprod(X.design)))
  
  # for running rune's vcov function
  ## TODO: UNCOMMENT (commented because of use ...type3.red)
  if(is.vector(Lc)) {
    Lc <- Lc[rho$Xlist$nums.Coefs]
  } else {
    if(type == 3) Lc <- Lc[, rho$Xlist$nums.Coefs]
  }

  if(ddf=="Kenward-Roger")
    return(calcKRMultDF(rho, Lc)) ## apply KR approximation
  else
    return(calcSatt(rho, Lc))  ## apply satterthwaite's approximation of ddf
}


################################################################################
## function to calculate F stat and pvalues for a given term. MAIN
################################################################################
calcFpvalueMAIN <- function(term, L, rho, ddf, type)
{
  if(type == 3)
  {
    Lc <- makeContrastType3SAS(rho$model, term, L)    
    #non identifiable because of rank deficiency
    if(!length(Lc))
      return(list(denom=0, Fstat=NA, pvalue=NA, ndf=NA, ss = NA, ms = NA))
    else
      return(c(calcFpvalueSS(Lc, rho, ddf, type), list(name=term)))       
  }
  if(type == 1)
  {
    Lc <- L[[term]]
    return(c(calcFpvalueSS(Lc, rho, ddf, type), list(name=term))) 
  } 
  if(type == 2){
    Lc <- L[[term]]
    return(c(calcFpvalueSS(Lc, rho, ddf, type), list(name=term))) 
  }
}


################################################################################
## initialize anova table for F test 
################################################################################
initAnovaTable <- function(model, test.terms, isFixReduce)
{
  anova.table <- matrix(NA, nrow=length(test.terms), ncol=6)
  rownames(anova.table) <- test.terms
  colnames(anova.table) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF","F.value", 
                             "Pr(>F)")
  anm <- anova(model, ddf="lme4")
  colnames(anm) <- c("NumDF", "Sum Sq", "Mean Sq", "F.value")
  
  
  anova.table[rownames(anm), c("NumDF", "Sum Sq", "Mean Sq")] <- 
    as.matrix(anm[, c("NumDF", "Sum Sq", "Mean Sq")])
  
  if(isFixReduce)
  {
    if(nrow(anova.table)==1)
    {
      anova.table <- c(anova.table[,1:5], 0, anova.table[,6])
      anova.table <- matrix(anova.table, nrow=1, ncol=length(anova.table))
      colnames(anova.table) <- c("Sum Sq", "Mean Sq", "NumDF", "DenDF", 
                                 "F.value", "elim.num","Pr(>F)")
      rownames(anova.table) <- rownames(anm)
      return(anova.table)
    }
    elim.num <- rep(0, nrow(anova.table))
    anova.table <- cbind(anova.table[,1:5], elim.num, anova.table[,6])
    colnames(anova.table)[7] <- "Pr(>F)"
  }
  return(anova.table)    
}


## sum of squares based on the proc glm SAS
getSS <- function(L, coef, XtX.) {
  L.beta <- L %*% coef
  if(is.vector(L))
    var.L.beta <- t(L) %*% XtX. %*% L
  if(is.matrix(L))
    var.L.beta <- L %*% XtX. %*% t(L)
  ss <- c(t(L.beta) %*% ginv(var.L.beta) %*% L.beta)
  ss
}


################################################################################
## fill anova table
################################################################################
fillAnovaTable <- function(result, anova.table)
{
  for (i in 1:length(result))
  {
    if(!is.null(result[[i]]$name) && 
       !( result[[i]]$name %in% rownames(anova.table)))
      next
    anova.table[result[[i]]$name, 4] <- result[[i]]$denom
    anova.table[result[[i]]$name, 5] <- result[[i]]$Fstat
    anova.table[result[[i]]$name, which(colnames(anova.table)=="Pr(>F)")] <- 
      result[[i]]$pvalue
    if(!is.na(result[[i]]$ss)){
      anova.table[result[[i]]$name, "Sum Sq"] <- result[[i]]$ss
      anova.table[result[[i]]$name, "Mean Sq"] <- result[[i]]$ms
    }
    
  }
  anova.table <- as.data.frame(anova.table)
  anova.table$NumDF <- as.integer(anova.table$NumDF)
  anova.table
}
