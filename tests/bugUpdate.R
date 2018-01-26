require(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

#### ASSERTS ERRORS
m <- lmer(Coloursaturation ~ TVset*Picture +
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)
m2 <- update(m)
assertError(stopifnot(class(m)==class(m2), TRUE), TRUE)




