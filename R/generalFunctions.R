step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, 
                 alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, 
                 fixed.calc=TRUE ,lsmeans.calc=TRUE, difflsmeans.calc=TRUE, 
                 test.effs=NULL, keep.effs = NULL,...)
{  
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  
  ddf <- checkNameDDF(ddf)
  
  result <- stepFun(model = model, ddf = ddf , type = type,  
                                  alpha.random = alpha.random, 
                                  alpha.fixed = alpha.fixed,
                                  reduce.fixed = reduce.fixed, 
                                  reduce.random = reduce.random,
                                  fixed.calc = fixed.calc, 
                                  lsmeans.calc = lsmeans.calc,
                                  difflsmeans.calc = difflsmeans.calc, 
                                  test.effs = test.effs,
                                  keep.effs = keep.effs, change.contr = TRUE)
  class(result) <- "step"
  result
}



print.step <- function(x, ...)
{
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n") 
    x$rand.table[,"p.value"] <- format.pval(x$rand.table[,"p.value"], digits=4, 
                                            eps=1e-7)
    x$rand.table[,"Chi.sq"] <- round(x$rand.table[,"Chi.sq"], 2)
    print(x$rand.table)     
  } 
  if(is.null(x$anova.table)){
    
  }else{
    if(nrow(x$anova.table) != 0)
    {
      if(class(x$model) == "lm" | class(x$model) == "gls")
      {
        cat("\nFixed effects:\n")
        print(x$anova.table)
        cat("\nLeast squares means:\n")
        print(x$lsmeans.table) 
        cat("\nFinal model:\n")
        print(x$model)
        return()
      }
      else
      {
        cat("\nFixed effects:\n")
        x$anova.table[,"Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"], 
                                                digits=4, eps=1e-7)
        x$anova.table[,c("Sum Sq","Mean Sq", "F.value")] <- 
          round(x$anova.table[,c("Sum Sq","Mean Sq", "F.value")],4)
        x$anova.table[,"DenDF"] <- round(x$anova.table[,"DenDF"],2)
        print(x$anova.table)          
        if(!is.null(x$lsmeans.table))
        {
          cat("\nLeast squares means:\n")
          printCoefmat(x$lsmeans.table, dig.tst=3 ,
                       tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),
                                 which(colnames(x$lsmeans.table)=="DF")), 
                       digits=3 , P.values = TRUE, has.Pvalue=TRUE)
        }
        if(!is.null(x$diffs.lsmeans.table))
        {
          cat("\n Differences of LSMEANS:\n")
          printCoefmat(x$diffs.lsmeans.table, dig.tst=1  ,
                       tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)==
                                            "Estimate")-1),
                                 which(colnames(x$diffs.lsmeans.table)=="DF")),
                       digits=3 , P.values = TRUE, has.Pvalue=TRUE)
        }
        
      }    
    }
    else
      print(x$anova.table)
  }
  
  cat("\nFinal model:\n")
  print(x$model@call) 
}


plot.step <- function(x, main = NULL, cex = 1.4, 
                      which.plot = c("LSMEANS", "DIFF of LSMEANS"),
                      effs = NULL, mult = TRUE, ...)
{
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0 && ("LSMEANS" %in% which.plot)){
    if(length(which.plot) == 1 && which.plot == "LSMEANS")
      return(plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS", main = main, cex = cex,
                         effs = effs, mult = mult))
  }         
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0 
     && ("DIFF of LSMEANS" %in% which.plot))
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS", 
                main = main, cex = cex, effs = effs, mult = mult)
}


lmer <- function(formula, data = NULL, REML = TRUE,
                 control = lmerControl(), start = NULL, verbose = 0L,
                 subset, weights, na.action, offset, contrasts = NULL,
                 devFunOnly = FALSE, ...)
{
  mc <- match.call()
  mc[[1]] <- quote(lme4::lmer)
  model <- eval.parent(mc)
  if(inherits(model, "merMod"))
    model <- as(model,"merModLmerTest")    
  return(model)
}




setMethod("anova", signature(object="merModLmerTest"),
          function(object, ..., ddf="Satterthwaite", type=3)  
          {
            mCall <- match.call(expand.dots = TRUE)
            dots <- list(...)
            modp <- if (length(dots))
              sapply(dots, is, "merModLmerTest") | sapply(dots, is, "merMod") | 
              sapply(dots, is, "lm") else logical(0)
            if (any(modp)) {
              return(callNextMethod())
            }
            else
            {
              cnm <- callNextMethod()
              if(!is.null(ddf) &&  ddf=="lme4") 
                return(cnm)              
                {
                  table <- cnm 
                  
                  ## errors in specifying the parameters
                  ddf <- checkNameDDF(ddf)
                  an.table <- tryCatch({calcANOVA(model=object, ddf=ddf, type=type)}
                                       , error = function(e) { NULL })
                  if(!is.null(an.table))
                  {
                    table <- an.table
                    
                    attr(table, "heading") <- 
                      paste("Analysis of Variance Table of type", as.roman(type) ,
                            " with ", ddf, 
                            "\napproximation for degrees of freedom")
                  }
                  else
                    message("anova from lme4 is returned\nsome computational error has occurred in lmerTest")
                  
                  
                  
                  class(table) <- c("anova", "data.frame")
                  return(table)
                }  

            }

          })

setMethod("summary", signature(object = "merModLmerTest"),
          function(object, ddf="Satterthwaite", ...){
            if(!is.null(ddf) && ddf=="lme4"){
              if(class(object) == "merModLmerTest")
                return(summary(as(object, "lmerMod")))
              #return(cl)
            }else{
              ## commented callNextMethod
              ## since it produces warning, summary cannot have multiple arguments
              ##cl <- callNextMethod()
              if(class(object) == "merModLmerTest")
                cl <- summary(as(object, "lmerMod"))
              #errors in specifying the parameters
              ddf <- checkNameDDF(ddf)
              
              tsum <- tryCatch( {calcSummary(object, ddf)}, 
                                error = function(e) { NULL })
              if(is.null(tsum)){
                message("summary from lme4 is returned\nsome computational error has occurred in lmerTest")
                return(cl)
              }
              coefs.satt <- cbind(cl$coefficients[,1:2, drop = FALSE], tsum$df, 
                                  tsum$tvalue, tsum$tpvalue)               
              cl$coefficients <- coefs.satt
              colnames(cl$coefficients)[3:5] <- c("df","t value","Pr(>|t|)")              
            }   
            
            cl$methTitle <- paste(cl$methTitle,  "\nt-tests use ", ddf, 
                                  "approximations to degrees of freedom")
            return(cl)
          })


calcSatterth <- function(model, L){
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  rho <- list()  ## vector containing info about model
  rho <- rhoInit(rho, model) ## save lmer outcome in rho vactor
  rho$A <- calcApvar(rho)
  calcSatt(rho, L, calcSS = FALSE)
}


rand <- function(model, ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  result <- testrand(model, reduce.random = FALSE, keep.effs = NULL,
                                 alpha.random = 0.1)
  res <- list(rand.table=result)
  class(res) <- "rand"
  res
}

print.rand <- function(x, ...)
{
  
  cat("Analysis of Random effects Table:\n")
  if(!is.null(x))
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,
                 tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),
                           which(colnames(x$rand.table)=="elim.num")), 
                 P.values=TRUE, has.Pvalue=TRUE)        
}

lsmeans <- function(model, test.effs=NULL, ddf = "Satterthwaite", ...){
  .Deprecated("lsmeansLT", package = "lmerTest") 
  lsmeansLT(model, test.effs=test.effs, ddf = ddf)
}


lsmeansLT <- function(model, test.effs=NULL, ddf = "Satterthwaite", ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  
  ddf <- checkNameDDF(ddf) 
  result <- lsmeans.calc(model, 0.05, test.effs = test.effs, 
                           lsmeansORdiff = TRUE, ddf)
  res <- list(lsmeans.table=result$summ.data, response=result$response)
  class(res) <- "lsmeansLT"
  res 
}

print.lsmeansLT <- function(x, ...)
{
  
  cat("Least Squares Means table:\n")
  printCoefmat(data.matrix(x$lsmeans.table), dig.tst=1, 
               tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),
                         which(colnames(x$lsmeans.table)=="DF")), digits=3 , 
               P.values=TRUE, has.Pvalue=TRUE)       
}

# plot.lsmeans <- function(x, main = NULL, cex = 1.4, effs = NULL, mult = TRUE, ...){
#   .Deprecated("plot.lsmeansLT", package = "lmerTest") 
#   plot.lsmeansLT(x, main = main, cex = cex, effs = effs, mult = mult)
# }

plot.lsmeansLT <- function(x, main = NULL, cex = 1.4, effs = NULL, mult = TRUE, ...)
{
  
  #plots for LSMEANS
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS", main = main, cex = cex,
                effs = effs,  mult = mult)     
}

difflsmeans <- function(model, test.effs=NULL, ddf = "Satterthwaite", ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  
  ddf <- checkNameDDF(ddf)
  result <- lsmeans.calc(model, 0.05, test.effs = test.effs, 
                         lsmeansORdiff = FALSE, ddf)  
  res <- list(diffs.lsmeans.table=result$summ.data, 
              response=result$response)
  class(res) <- "difflsmeans"
  res 
}


print.difflsmeans <- function(x, ...)
{
  
  cat("Differences of LSMEANS:\n")
  printCoefmat(data.matrix(x$diffs.lsmeans.table), dig.tst=1, 
               tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),
                         which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3,
               P.values=TRUE, has.Pvalue=TRUE)
  
}



plot.difflsmeans <- function(x, main = NULL, cex = 1.4, effs = NULL, 
                             mult = TRUE, ...)
{
  
  #plots for DIFF of LSMEANS
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS", 
                main = main, cex = cex, effs = effs, mult = mult)   
}