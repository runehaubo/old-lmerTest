require(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

m <- lme4::lmer(Informed.liking ~ Product + (1|Consumer) , data = ham)

initRho <- function(model){
  lctr <- as.list(rep("contr.SAS",1))
  names(lctr) <- "Product"
  rho <- list()
  model <- update(model, contrasts = lctr)
  return(model)
}

m2 <- initRho(m)

assertError(update(m2))
