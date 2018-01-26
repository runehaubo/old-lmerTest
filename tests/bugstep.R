require(lmerTest)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...) 
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

m <- lmer(Coloursaturation ~ TVset*Picture+
            (1|Assessor)+(1|Assessor:TVset), data=TVbo)
step(m)

#does not with lmerTest attached
m.lm <- lm(Coloursaturation ~ TVset*Picture, data=TVbo)
assertError(step(m.lm))

#stats package needs to implicitely specified then
stats::step(m.lm)


assertError(lsmeans(m.lm))

