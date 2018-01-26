require(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

## rand step lsmeans difflsmeans should work only with inheritance of lmerMod class

## TODO uncomment the following whenever fixed deepcopy thing
## with the lme4
ifTest <- TRUE

if(ifTest){
  gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
               data = cbpp, family = binomial)
  
  ## should not work with glmer models
  assertError(rand(gm1))
  assertError(step(gm1))
  assertError(lsmeans(gm1))
  assertError(difflsmeans(gm1))
  
  ## should not work with nlmer models
  # startvec <- c(Asym = 200, xmid = 725, scal = 350)
  # nm1 <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
  #               Orange, start = startvec)
  # assertError(rand(nm1))
  # assertError(step(nm1))
  # assertError(lsmeans(nm1))
  # assertError(difflsmeans(nm1))
  
  
  ## should work with lmer from lme4 package (class lmerMod)
  m <- lme4::lmer(Coloursaturation ~ TVset*Picture+
                    (1|Assessor), data=TVbo)
  rand(m)
  step(m)
  lsmeansLT(m)
  difflsmeans(m)
}


