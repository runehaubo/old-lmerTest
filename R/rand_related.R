################################################################################
## functions related to testing random effects
################################################################################

## the main function to test / reduce random effects
testrand <- function(model, reduce.random = FALSE, keep.effs = NULL,
                     alpha.random = 0.1){
 result.rand <- elimRandEffs(model, alpha.random, reduce.random) 
  model <- result.rand$model
  ## convert rand table to data frame
  rt <- as.data.frame(result.rand$TAB.rand)
  rt$Chi.DF <- as.integer(rt$Chi.DF)
  if(!is.null(rt$elim.num))
    rt$elim.num <- as.integer(rt$elim.num)
  rt
}


## function performs LRT test on random effect
.doLRT <- function(rnm, rand.terms, fmodel, model){
  nt <- rnm$gr.part
  rnm <- list(rnm)
  names(rnm) <- nt
  
  fm <- paste(fmodel)
  mf.final <- fmElimRandTerm(rnm, rand.terms, fm)
  is.present.rand <- checkPresRandTerms(mf.final)
  
  # no more random terms in the model
  if(!is.present.rand)
  {
    return(compareMixVSFix(model, mf.final, data, rnm))
    
  } 
  model.red <- update(model, mf.final, getME(model, "is_REML"))
  anova.red <- suppressMessages(anova(model, model.red, refit = FALSE))
  return(saveInfoForTerm(rnm, anova.red$Chisq[2], anova.red$"Chi Df"[2], 
                         anova.red$'Pr(>Chisq)'[2], model.red))
}



## eliminate NS random terms 
elimRandEffs <- function(model, alpha, reduce.random, 
                         keep.effs = NULL)
{
  isInitRand <- TRUE
  elim.num <- 1
  stop <- FALSE
  while(!stop)
  {
    fmodel <- formula(model)    
    rand.terms <- sapply(getRandTerms(fmodel), 
                         function(x) substr(x,2,nchar(x)-1), USE.NAMES = FALSE) 
    rand.terms.table <- getRandTermsTable(rand.terms)
    
    if(isInitRand)
    {
      rand.table <- initRandTable(names(rand.terms.table), reduce.random)
      isInitRand <- FALSE
    }      
    
    fm <- paste(fmodel)    
    
    
    infoForTerms <- llply(rand.terms.table, .fun = .doLRT, rand.terms, fmodel, 
                          model)
    
    ## find the maximal p-value if the reduction is required
    if(reduce.random)     
      infoForTermElim <- .findElimTerm(infoForTerms, keep.effs)    
    else{
      rand.table <- updateRandTable(infoForTerms, rand.table, 
                                    reduce.random = reduce.random)
      model.last <- model
      break
    }
    
    if(length(infoForTermElim)!=0 && infoForTermElim[[1]]$pv > alpha)
    {
      rand.table <- updateRandTable(infoForTermElim, rand.table, elim.num, 
                                    reduce.random)
      elim.num <- elim.num + 1      
    }
    else
    {
      rand.table <- updateRandTable(infoForTerms, rand.table, 
                                    reduce.random=reduce.random)
      model.last <- model
      break
    }
    
    model <- infoForTermElim[[1]]$model.red # model.final 
    
    if(is(model, "lm")){
      model.last <- model
      break
    }
    
    
  }
  return(list(model = model.last, TAB.rand = rand.table))
}


## compare mixed model versus fixed
compareMixVSFix <- function(model, mf.final, data, name.term)
{
  
  model.red <- refitLM(model)
  
  l.fix <- -2*logLik(model, REML=TRUE)[1]
  l.red <- -2*logLik(model.red, REML=TRUE)[1]
  
  p.chisq <- 1 - pchisq (l.red -l.fix ,1)
  infoForTerm <- saveInfoForTerm(name.term, l.red -l.fix, 1, p.chisq, 
                                 model.red = model.red)
  return(infoForTerm)
}



## fill a row for the random matrix
fillRowRandTable <- function(term, rand.table, elim.num, reduce.random)
{
  #   nrow.term <- which(rownames(rand.table)==names(term))
  #   um"] <- elim.num  
  #   return(rand.table)
  
  #
  nrow.term <- which(gsub(" ","",rownames(rand.table))==names(term))
  rand.table[nrow.term, "Chi.sq"] <- term[[1]]$chisq
  rand.table[nrow.term, "Chi.DF"] <- term[[1]]$chisq.df
  rand.table[nrow.term, "p.value"] <- term[[1]]$pv
  if(reduce.random)
  {
    rand.table[nrow.term, "elim.num"] <- elim.num 
    if(term[[1]]$chisq.df==2 && elim.num!=0)
    {
      rand.table.upd <- matrix(0,ncol=4, nrow=1)
      colnames(rand.table.upd) <- c("Chi.sq", "Chi.DF", "elim.num", "p.value")
      
      rownames(rand.table.upd) <- paste(paste(rep(" ", 
                                                  substring.location( names(term), 
                                                                      term[[1]]$term[[1]]$gr.part)$first-2),
                                              collapse=""), term[[1]]$term[[1]]$gr.part, 
                                        collapse="")
      if(nrow.term==nrow(rand.table))
      {
        rand.table <- rbind(rand.table, rand.table.upd)
        #rownames(rand.table)[1] <- term$term
      }
      else
      {
        rnames <- c(rownames(rand.table)[1:nrow.term], rownames(rand.table.upd),
                    rownames(rand.table)[(nrow.term+1):nrow(rand.table)])
        rand.table <- rbind(rand.table[1:nrow.term,], rand.table.upd, 
                            rand.table[(nrow.term+1):nrow(rand.table),])  
        rownames(rand.table) <- rnames 
      } 
    }      
    
  }
  
  rand.table  
}

## update table for random terms
updateRandTable <- function(infoForTerm, rand.table, 
                            elim.num=0, reduce.random)
{  
  if(!is.null(infoForTerm$term))
  {   
    rand.table <- fillRowRandTable(infoForTerm, rand.table, elim.num, 
                                   reduce.random)  
  } 
  else
  {
    for(i in 1:length(infoForTerm))
    {
      iterm <- infoForTerm[i]
      rand.table <- fillRowRandTable(iterm, rand.table, 
                                     elim.num, reduce.random) 
    }
  }  
  rand.table
}


## function get the slope part of a random term
findSlopePart <- function(term)
{
  sub.loc.div <- substring.location(term," |")
  slopepart <- substring2(term,1,sub.loc.div$first)
  grouppart <- substring2(term,sub.loc.div$last, nchar(term))
  parts <- unlist(strsplit(slopepart, split=c("[[:punct:]]")))
  parts <- gsub(" ","", parts , fixed=TRUE)
  parts
}

getSlGrParts <- function(term)
{
  randTerm.split <- unlist(strsplit(term, "\\|")) 
  sl.part <- sapply( strsplit(randTerm.split[1], split=c("[[:punct:]]")), 
                     function(x) gsub(" ","", x , fixed=TRUE))
  gr.part <- gsub(" ","", randTerm.split[2] , fixed=TRUE)
  return(list(sl.part=sl.part, gr.part=gr.part))
}


getRowNamesForRandTable <- function(term)
{
  sl.part <- getSlGrParts(term)$sl.part
  gr.part <- getSlGrParts(term)$gr.part
  if(length(sl.part) > 1 || sl.part !="1")
  {
    sl.part <- sl.part[! sapply(sl.part, function(x) x=="1" || x=="0")]
    listterms <- function(x)
    {
      lst <- list(sl.part=x, gr.part=gr.part)
      lst
    }
    la <- lapply(sl.part, listterms)
    names(la) <- sapply(sl.part, 
                        function(x) paste(c(x, gr.part), collapse=":"), USE.NAMES=FALSE)     
  }
  else
  {
    lst <- list(sl.part=sl.part, gr.part=gr.part) 
    la <- list(lst)
    names(la) <- gr.part     
  } 
  
  la
}

## initialize table for random terms
initRandTable <- function(names, reduce.random)
{ 
  
  if(reduce.random)
  {
    rand.table <- matrix(0,ncol=4, nrow=length(names))
    colnames(rand.table) <- c("Chi.sq", "Chi.DF", "elim.num", "p.value")
    rownames(rand.table) <- names
  }
  else
  {
    rand.table <- matrix(0,ncol=3, nrow=length(names))
    colnames(rand.table) <- c("Chi.sq", "Chi.DF", "p.value")
    rownames(rand.table) <- names
  }
  return(rand.table)
}


## return formula with the simplified random structure
fmElimRandTerm <- function(rnm, rand.terms, fm)
{
  
  nSameGr <- length(rand.terms[sapply(rand.terms, 
                                      function(x) (getSlGrParts(x)$gr.part==rnm[[1]]$gr.part))])
  term.simplify <-  rand.terms[sapply(rand.terms, 
                                      function(x) (getSlGrParts(x)$gr.part==rnm[[1]]$gr.part) && 
                                        (rnm[[1]]$sl.part %in% findSlopePart(x)))] 
  sub.loc.div <- substring.location(term.simplify," |")
  slopepart <- gsub(" ","", substring2(term.simplify, 1, sub.loc.div$first) , 
                    fixed=TRUE)
  mf.term.simpl <- as.formula(paste(fm[2],fm[1],paste(slopepart, "-", 
                                                      rnm[[1]]$sl.part), sep=""))
  mf.term.simpl <- update.formula( mf.term.simpl, mf.term.simpl)
  if(rnm[[1]]$sl.part != "1" && as.character(mf.term.simpl[3])!="1 - 1" 
     && nSameGr==1)
  {
    
    fm[3] <- paste(fm[3], "-", paste("(",term.simplify,")", sep=""), "+" , 
                   paste("(",as.character(mf.term.simpl[3]),
                         "|",rnm[[1]]$gr.part, ")", sep=""))
  }
  else
  {
    fm[3] <- paste(fm[3], "-", paste("(",term.simplify,")", sep=""), sep="")
  }
  
  mf.final <-  as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  mf.final <- update.formula( mf.final, mf.final)
  mf.final
}


getRandTermsTable <- function(rand.terms)
{
  rand.terms.table <- NULL
  for(rand.term in rand.terms)
  {
    rand.terms.table <- c(rand.terms.table, getRowNamesForRandTable(rand.term))
  }
  rand.terms.table
}

  
## save info (pvalues, chisq val, std, var...) for term
saveInfoForTerm <- function(term, chisq, chisq.df, pv, model.red = NULL)
{
  term.info <- NULL 
  term.info$pv <- pv
  term.info$chisq <- chisq
  term.info$chisq.df <- chisq.df
  term.info$term <- term
  if(!is.null(model.red))
    term.info$model.red <- model.red
  return(term.info)
}


.findElimTerm <- function(infoForTerms, keep.effs){  
  infoForTermsUpd <- infoForTerms[setdiff(names(infoForTerms), keep.effs)]
  ind <- which.max(lapply(infoForTermsUpd, function(x) x$pv))
  return(infoForTermsUpd[ind])  
}



## the same function as in SensMixed package
## list of random and fixed terms of the model
.fixedrand <- function(model)
{  
  effs <- attr(terms(formula(model)), "term.labels")
  neffs <- length(effs)
  randeffs <- effs[grep(" | ", effs)]
  randeffs <- sapply(randeffs, function(x) substring(x, 5, nchar(x)))
  fixedeffs <- effs[!(effs %in% names(randeffs))]
  return(list(randeffs=randeffs, fixedeffs=fixedeffs))
}



## create reduce slopes model
createModelRedSlopes <- function(x, term, fm, model)
{
  fm[3] <- paste(fm[3], "-", term, "+" , paste(x, collapse="+"))
  mf.final <-  as.formula(paste(fm[2],fm[1],fm[3], sep=""))
  mf.final <- update.formula(mf.final,mf.final)
  model.red <- updateModel(model, mf.final, change.contr = FALSE)
  # model.red <- update(model, formula = mf.final, 
  #                     contrasts = l.lmerTest.private.contrast)
  return(model.red)
}

## get the random terms
getRandTerms <- function(fmodel)
{
  terms.fm <- attr(terms(fmodel),"term.labels")
  ind.rand.terms <- which(unlist(lapply(terms.fm,
                                        function(x) 
                                          substring.location(x, "|")$first))!=0)
  return(unlist(lapply(terms.fm[ind.rand.terms],
                       function(x) paste("(",x,")",sep=""))))
}

## check if there are no random terms in the model
checkPresRandTerms <- function(mf.final)
{
  sub.loc.rand <- substring.location(paste(mf.final)[3], "|")
  if(sub.loc.rand$first==0 && sub.loc.rand$last==0)
    return(FALSE)
  return(TRUE)
}

## check if there are correlations between intercepts and slopes
checkCorr <- function(model)
{
  corr.intsl <- FALSE
  modelST <- getME(model, "ST")
  lnST <- length(modelST)
  for(i in 1:lnST)
  {    
    if(nrow(modelST[[i]])>1)
      corr.intsl <- TRUE
  } 
  return(corr.intsl) 
}

