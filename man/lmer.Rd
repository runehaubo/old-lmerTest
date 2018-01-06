\name{lmer}
\alias{lmer}
\title{Fit Linear Mixed-Effects Models}
\description{
  Fit a linear mixed model
}

\details{

  This \code{lmer} function is an overloaded function of \code{lmer} (\code{merMod} class from \pkg{lme4} package).
    
}
\value{
  An object of class \code{"\linkS4class{merModLmerTest}"} 
}
\seealso{
  \code{\linkS4class{merModLmerTest}} class
}
\examples{

library(lmerTest)

## linear mixed models
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

# anova table the same as of class merMod but with additional F statistics and  
# p-values calculated based on Satterthwaite's approximations
anova(fm1)

# anova table the same as of class merMod but with additional F statistics and  
# p-values calculated based on Kenward-Roger's approximations
\dontrun{
if(requireNamespace("pbkrtest", quietly = TRUE))
anova(fm1, ddf="Kenward-Roger")

# anova table the same as of class merMod
anova(fm1, ddf="lme4")
}

# gives summary of merModLmerTest class. The same as of class merMod but with
# additional p-values calculated based on Satterthwate's approximations
summary(fm1)

## multiple comparisons statistics. The one from lme4 package
\dontrun{
anova(fm1, fm2)
}
}
\keyword{models}
\keyword{methods}
