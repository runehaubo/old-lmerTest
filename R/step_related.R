stepFun <- function(model, ddf = "Satterthwaite", type = 3, 
                    alpha.random = 0.1, alpha.fixed = 0.05, 
                    reduce.fixed = TRUE, reduce.random = TRUE, 
                    fixed.calc = TRUE, lsmeans.calc = TRUE, 
                    difflsmeans.calc = TRUE, 
                    test.effs = NULL, keep.effs = NULL, 
                    change.contr = TRUE)
{

  ## check type of hypothesis
  if(!(type %in% c(1,2,3)))  
    stop('Parameter type is wrongly specified') 
  
  ## check keep.effs 
  if(!is.null(keep.effs)){
    model.effs <- .fixedrand(model)
    keep.effs1 <- .getKeepEffs(keep.effs, model.effs) 
    if(length(unlist(keep.effs1)) == 0)
      message(paste("No ", keep.effs, "exist among effects in the model"))
    keep.effs <- keep.effs1
  }    
  
  result <- NULL
  result$response <- rownames(attr(terms(model),"factors"))[1]
  result$call <- model@call ## save the call of the model 
  result$corr.intsl <- checkCorr(model)   
  anova.table <- NULL
  
  ## save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
  {
    result <- saveResultsFixModel(result, model)
    result$rand.table <- NULL
    return(result)
  }
  
  ## analysis of the random part  
  result.rand <- elimRandEffs(model, alpha.random, reduce.random, 
                              keep.effs$randeffs)  
  model <- result.rand$model
  ## convert rand table to data frame
  rt <- as.data.frame(result.rand$TAB.rand)
  rt$Chi.DF <- as.integer(rt$Chi.DF)
  if(!is.null(rt$elim.num))
    rt$elim.num <- as.integer(rt$elim.num)
  
  result$rand.table <- rt
  if(!fixed.calc){
    result$model <- model
    return(result)
  }
  
  ## save results for fixed effects for model with only fixed effects
  if(class(model) == "lm" | class(model) == "gls")
    return(saveResultsFixModel(result, model, type))
 
  ## perform reduction of fixed effects for model with mixed effects
  stop = FALSE
  is.first.anova <- TRUE
  is.first.sign <- TRUE  
  rho <- list() ## environment containing info about model
  rho <- rhoInit(rho, model, change.contr)
  
  while(!stop) {      
    ## if there are no fixed terms
    if(nrow(anova(model, ddf="lme4"))==0) {
      if(is.null(anova.table)) return(emptyResult(result, model))
      break
    }   
    
    rho$A <- calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma  
    rho$Xlist <- createDesignMat(rho) ## X design matrix for fixed effects
  
    test.terms <- attr(terms(rho$model),"term.labels") ## define the terms that are to be tested
    
    ## initialize anova table
    if(is.first.anova) {
      anova.table <- initAnovaTable(rho$model, test.terms, reduce.fixed) 
      is.first.anova <- FALSE
      elim.num <- 1
    }
    
    ## calculate general set of hypothesis matrix 
    L <- calcGeneralSet(rho, type) 

    ## compute (list of) F-values etc.:
    resultFpvalueSS <- lapply(test.terms, function(tt)
      calcFpvalueMAIN(term=tt, L=L, rho=rho, ddf=ddf, type=type))
    
    ## fill anova table
    anova.table <- fillAnovaTable(resultFpvalueSS,  anova.table)
    
    if(!reduce.fixed) break     
    else {
      resNSelim <- elimNSFixedTerm(rho, anova.table, alpha.fixed, 
                                   elim.num, keep.effs$fixedeffs, change.contr)
      if(is.null(resNSelim))
        break
      else {
        rho <- resNSelim$rho
        anova.table <- updateAnovaTable(resNSelim)
        elim.num <- elim.num + 1
      }        
    }
  }
  ## convert anova table to data frame
  anova.table <- as.data.frame(anova.table)
  anova.table$NumDF <- as.integer(anova.table$NumDF)
  if(!is.null(anova.table$elim.num))
    anova.table$elim.num <- as.integer(anova.table$elim.num)
  result$anova.table <- anova.table

  ## if in step function least squares means of diffs of LSMEANS are required
  if(lsmeans.calc) {
    lsmeans.tab <- calcLSMEANS(rho, alpha.fixed, test.effs = test.effs,
                               lsmeansORdiff = TRUE)
    result$lsmeans.table <- lsmeans.tab$summ.data
  } else {
    result$lsmeans.table <- NULL
  }
  if(difflsmeans.calc) {
    lsmeans.tab <- calcLSMEANS(rho, alpha.fixed, test.effs = test.effs, 
                               lsmeansORdiff=FALSE)
    result$diffs.lsmeans.table <- lsmeans.tab$summ.data
  } else {
    result$diffs.lsmeans.table <- NULL
  }
  
  ## format anova.table and random.table according to elim.num column
  result$anova.table <- formatElimNumTable(result$anova.table) 
  result$rand.table <- formatElimNumTable(result$rand.table) 
  
  ## save model
  if(inherits(rho$model, "merMod"))
    model <- as(rho$model,"merModLmerTest")
  result$model <- rho$model
  return(result)
}


################################################################################
## find NS effect from the model (starting from highest order interactions)
################################################################################
getNSFixedTerm <- function(anova.table, alpha, keep.effs = NULL) {
  
  pv.max <- 0
  
  if(length(which(anova.table[,"elim.num"]==0))==1)
    terms.compare <- rownames(anova.table)[anova.table[,"elim.num"]==0]
  else
    terms.compare <- getTermsToCompare(anova.table[anova.table[,"elim.num"]==0,], 
                                       keep.effs)
  
  for(tcmp in terms.compare)
  {
    if((!tcmp %in% rownames(anova.table)) || is.na(anova.table[tcmp,"Pr(>F)"]))
      next
    ind <- which(rownames(anova.table)==tcmp)
    if(anova.table[ind, which(colnames(anova.table)=="Pr(>F)")]>=pv.max)
    {
      ns.term <- tcmp
      pv.max <- anova.table[ind,which(colnames(anova.table)=="Pr(>F)")]
    }
  }  
  if(pv.max >= alpha)
    return(ns.term)  
  else
    return(NULL)  
}


################################################################################
## eliminate NS fixed effect from the model
################################################################################
elimNSFixedTerm <- function(rho, anova.table, alpha, elim.num, 
                            keep.effs = NULL, change.contr)
{
  ns.term <- getNSFixedTerm(anova.table, alpha, keep.effs = keep.effs)
  if( is.null(ns.term) )
    return(NULL)
  anova.table[ns.term, "elim.num"] <- elim.num
  fm <- formula(rho$model)
  fm[3] <- paste(fm[3], "-", ns.term)
  
  mf.final <- as.formula(paste(fm[2], fm[1], fm[3], sep=""))

  rho <- rhoInit(rho, rho$model, change.contr = change.contr, mf.final = mf.final)

  
  return(list(rho = rho, anova.table = anova.table))
}

###############################################################################
## get terms to compare in anova.table
###############################################################################
getTermsToCompare <- function(anova.table, keep.effs = NULL) {
  #order.terms <- attr(terms(model),"order")
  #allterms <- attr(terms(model),"term.labels")
  anova.table.upd <- anova.table[complete.cases(anova.table), , drop=FALSE]
  order.terms <- orderterms(anova.table.upd)
  allterms <- rownames(anova.table.upd )
  ind.hoi <- which(order.terms == max(order.terms))
  ind.terms.contain <- getIndTermsContained(allterms, ind.hoi)
  
  #get the rest of the terms to compare
  allterms.rest <- allterms[-c(ind.terms.contain, ind.hoi)]
  if( length(allterms.rest)==0 )
    terms.compare <- allterms[ind.hoi]
  else
  {
    #get highest order terms in the remaining ones
    order.rest <- unlist(lapply(allterms.rest, function(x) 
      length(unlist(strsplit(x,":")))))
    ind.hoi.rest <- which(order.rest == max(order.rest))
    gtc <- getIndTermsContained(allterms.rest, ind.hoi.rest)
    if( !is.null(gtc) )
      terms.compare <- c(allterms[ind.hoi], 
                         allterms.rest[-getIndTermsContained(allterms.rest, 
                                                             ind.hoi.rest)])
    else
      terms.compare <- c(allterms[ind.hoi], allterms.rest)
  }
  
  return(setdiff(terms.compare, keep.effs))
}

###############################################################################
# get terms contained 
###############################################################################
getIndTermsContained <- function(allterms, ind.hoi)
{
  
  terms.hoi.split <- strsplit(allterms[ind.hoi],":")
  ind.terms.contain <- NULL
  #check which of the terms are contained in the highest order terms
  for(i in (1:length(allterms))[-ind.hoi]) 
  {
    isContained<-FALSE
    for(j in 1:length(terms.hoi.split))
    {
      #if the term is contained in some of the highest order interactions then 
      #we cannot test it for significance
      if(length(which(unlist(strsplit(allterms[i],":")) %in% 
                      terms.hoi.split[[j]] == FALSE))==0)
      {
        isContained <- TRUE
        break
      }                
    }
    if(isContained)
      ind.terms.contain <- c(ind.terms.contain,i)
    
  }
  # if there are no terms that are contained in the maximum order effects
  # then compare all the terms between each other for the maximum p value
  if( is.null(ind.terms.contain) )
    return(NULL)
  return(ind.terms.contain)
}



orderterms <- function(anova.table) {
  return(unlist(lapply(rownames(anova.table), function(x) 
    length(unlist(strsplit(x,":"))))))
}


.findKeepRandEff <- function(keep.eff, model.effs) {
  randTerms <- getRandTermsTable(names(model.effs))
  is.scalar <- unlist(lapply(randTerms, function(x) x$sl.part == "1"))
  randTermsScal <- randTerms[is.scalar]
  randTermsSl <- randTerms[!is.scalar]
  keep.eff.scal <- .findKeepEff(keep.eff, names(randTermsScal))
  keep.eff.slope <- .findKeepEffSlope(keep.eff, randTermsSl)
  c(keep.eff.scal, keep.eff.slope)
}


.getKeepEffs <- function(keep.effs, model.effs){
  #rand.terms.table <- getRandTermsTable(rand.terms)
  randeffs <- unlist(lapply(keep.effs, .findKeepRandEff, 
                            unlist(model.effs$randeffs)))
  fixedeffs <- unlist(lapply(keep.effs, .findKeepEff, 
                             unlist(model.effs$fixedeffs)))
  return(list(randeffs = randeffs, fixedeffs = fixedeffs))
}


.findKeepEffSlope <- function(keep.eff, randTermsSl) {
  ind.Slope <- unlist(lapply(randTermsSl, function(x) x$sl.part == keep.eff)) 
  if(sum(ind.Slope) == 0){
    ind.Gr <- unlist(lapply(randTermsSl, function(x) x$gr.part == keep.eff)) 
    if(sum(ind.Gr) > 0)
      return(c(names(randTermsSl[ind.Slope]), keep.eff))
  } 
  c(names(randTermsSl[ind.Slope]))
}


.getKeepInter <- function(model.eff, split.keep.eff) {
  if(setequal(unlist(strsplit(model.eff, ":")), split.keep.eff))
    return(model.eff)
}

.findKeepEff <- function(keep.eff, model.effs) {
  ## for main terms
  if(length(grep(":", keep.eff)) == 0){
    if(keep.eff %in% model.effs)
      return(keep.eff)
  }
  else{
    ## for interaction terms
    return(unlist(lapply(model.effs, .getKeepInter, 
                         unlist(strsplit(keep.eff, ":")))))
  }   
}


## format table according to elim.num column
formatElimNumTable <- function(table) {
  if("elim.num" %in% colnames(table)) {
    table[which(table[,"elim.num"]==0),"elim.num"] <- 1000
    table <- table[with(table, order(elim.num, decreasing=FALSE)),]
    table[,"elim.num"] <- as.character(table[,"elim.num"])
    table[which(table[,"elim.num"]=="1000"),"elim.num"] <- "kept"
  }
  return(table)
}


updateAnovaTable <- function(resNSelim){
  anm <- anova(resNSelim$rho$model)
  anova.table <- resNSelim$anova.table
  anova.table[rownames(anm), c("Sum Sq", "Mean Sq", "NumDF")] <-
    as.matrix(anm[, c("Sum Sq", "Mean Sq", "Df")])
  anova.table
}


emptyResult <- function(result, model){
  lsmeans.summ <-  matrix(ncol=7,nrow=0)
  colnames(lsmeans.summ) <- c("Estimate", "Standard Error", "DF", 
                              "t-value", "Lower CI", "Upper CI", "p-value")
  lsmeans.summ <- as.data.frame(lsmeans.summ)
  result$lsmeans.table <- lsmeans.summ
  result$diffs.lsmeans.table <- lsmeans.summ
  result$model <- model
  result$anova.table <- anova(model, ddf="lme4")
  result
}


