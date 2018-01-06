require(lmerTest)

m <- lme4::lmer(Informed.liking ~ Product + (1|Consumer) , data = ham)

initRho <- function(model){
  lctr <- as.list(rep("contr.SAS",1))
  names(lctr) <- "Product"
  rho <- list()
  model <- update(model, contrasts = lctr)
  return(model)
}

m2 <- initRho(m)

tools::assertError(update(m2))
