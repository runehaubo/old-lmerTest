testSat <- FALSE
testKen <- TRUE

if(testSat){
  ## checking on big data sets

  require(lmerTest)

  load(system.file("testdata", "ttDamon.RData", package="lmerTest"))

  testFirstR02.test=lmer(RT2LogR ~ condition*v_matri_plus_n_indef_freq_sq+
                           (1+condition|Subject),ttDamon, REML=TRUE)

  an.sat <- anova(testFirstR02.test)

  TOL <- 1e-2 # for the check

  ## with 4 decimals should agree with SAS output
  ## numbers before decimals should agree with SAS output
  stopifnot(
    all.equal(an.sat[,"DenDF"], c(76.206, 6074.6, 6146.1), tol = TOL),
    all.equal(an.sat[, "F.value"], c(9.2182, 2.8360, 6.6544), tol = TOL)
    , TRUE)

  ## check with calcSatterth function
  # L <- matrix(0, nrow=3, ncol=8)
  # L[1,2] <- L[2,3] <- L[3,4] <- 1
  # system.time(calcSatterth(testFirstR02.test, L))
  # if(require(pbkrtest))
  #   system.time(KRmodcomp(testFirstR02.test, L))

  if(testKen){
    if(require(pbkrtest))
      an.kenw <- anova(testFirstR02.test, ddf="kenw")
    stopifnot(ncol(an.kenw) == 6)
  }

}



