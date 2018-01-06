
##########################################################################
## Create rho vector containing info about  mixed model ##################
##########################################################################
rhoInit <- function(rho, model, change.contr = FALSE, mf.final = NULL)
{
  if(change.contr)
    rho$model <- updateModel(model, mf.final = mf.final, change.contr = change.contr)
  else
    rho$model <- model
  rho$fixEffs <- fixef(rho$model)
  rho$sigma <- sigma(rho$model)
  rho$thopt <- getME(rho$model, "theta")
  return(rho)  
}

       
############################################################################
#get formula for model 
############################################################################
getFormula <- function(model, withRand=TRUE)
{
  fmodel <- formula(model)  
 
  if(withRand)
    return(fmodel)
  
  fm <- paste(fmodel)
  fmodel.red <- paste(fm[2],fm[1], 
                        paste(fm[3], 
                              paste(unlist(lapply(names(.fixedrand(model)$randeffs),
                                                  function(x) paste("(",x, ")"))), 
                                    collapse = " - "), 
                              sep = "-"))
  return(update(fmodel, fmodel.red))
}



emptyAnovaLsmeansTAB <- function()
{
  result <- NULL
  anova.table <-  matrix(ncol=5,nrow=0)
  colnames(anova.table) <- c("Estimate", "Standard Error", "DF", "F-value", 
                             "p-value")
  anova.table <- as.data.frame(anova.table)
  result$TAB.fixed <- anova.table
  lsmeans.summ <-  matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ) <- c("Estimate", "Standard Error", "DF", "t-value", 
                              "Lower CI", "Upper CI", "p-value")
  lsmeans.summ <- as.data.frame(lsmeans.summ)
  result$TAB.lsmeans <- lsmeans.summ
  return(result)
}


## save results for fixed effects for model with only fixed effects
saveResultsFixModel <- function(result, model, type = 3)
{
  if(type==3)
    result$anova.table <- drop1(model, test="F")
  else{
    result$anova.table <- anova(model)
  }
  result$model <- model
  lsmeans.summ <-  matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                              "Lower CI", "Upper CI", "p-value")
  result$lsmeans.table <- lsmeans.summ
  result$diffs.lsmeans.table <- lsmeans.summ
  return(result)
}

getREML <- function(model)
{
  if(inherits(model,"merMod"))
    return(getME(model, "is_REML"))

}


updateModel <- function(object, mf.final = NULL, ..., change.contr = FALSE) {
  if (is.null(call <- getCall(object)))
    stop("object should contain a 'call' component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(mf.final))
    call$formula <- update.formula(formula(object), mf.final)
  if(any(grepl("sample", call))){
    call <- as.list(call)[-which(names(as.list(call)) %in% c("data", "subset"))]
    call[["data"]] <- quote(model.frame(object))
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
    }
  }
  if(change.contr){
    mm <- model.matrix(object)
    contr <- attr(mm,"contrasts")
    ## change contrasts for F tests calculations
    ## list of contrasts for factors
    if(change.contr && length(which(unlist(contr)!="contr.SAS")) > 0)
    {
      names.facs <- names(contr)
      l.lmerTest.private.contrast <- as.list(rep("contr.SAS",length(names.facs)))
      names(l.lmerTest.private.contrast) <- names(contr)
      call[["contrasts"]] <- l.lmerTest.private.contrast
    }
    else if(!is.null(contr)){
      call[["contrasts"]] <- contr[names(contr) %in% attr(terms(call$formula),"term.labels")]
    }
  }
  call <- as.call(call)
  ff <- environment(formula(object))
  pf <- parent.frame()  ## save parent frame in case we need it
  sf <- sys.frames()[[1]]
  ff2 <- environment(object)
  tryCatch(eval(call, envir=ff),
           error=function(e) {
             tryCatch(eval(call, envir=sf),
                      error=function(e) {
                        tryCatch(eval(call, envir=pf),
                         error=function(e) {
                           eval(call, envir=ff2)
                         })})})
}


checkNameDDF <- function(ddf){
  ddfs <- c("Satterthwaite", "Kenward-Roger")
  ind.ddf <- pmatch(tolower(ddf), tolower(ddfs))
  if(is.na(ind.ddf))  
    stop('Parameter ddf is wrongly specified')  
  else
    ddf <- ddfs[ind.ddf]
  ddf
}



checkForEstim <- function(mat, rho){
  if (!requireNamespace("estimability", quietly = TRUE)) 
    stop("estimability package required to check for estimability")
  contr <- attr(rho$model@pp$X, "contrasts")
  mm <- model.matrix(terms(rho$model), rho$model@frame, contrasts.arg=contr)
  nbasis <- estimability::nonest.basis(mm)
  return(which(apply(mat, 1, function(x) estimability::is.estble(x, nbasis))))
}


