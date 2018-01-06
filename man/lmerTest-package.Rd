\name{lmerTest-package}
\alias{lmerTest-package}
\alias{lmerTest}
\docType{package}
\title{
 The package performs different kinds of tests on lmer objects, such as \verb{F} 
 tests of  types I - III hypotheses for the fixed part, 
likelihood ratio tests for the random part, least squares means 
(population means) and differences of least squares means for the factors 
of the fixed part with corresponding plots. The package also provides with 
a function step, that preforms backward elimination of non-significant effects, 
starting from the random effects, and then fixed ones. 
}
\description{
 The package provides anova function, that gives data frame similar to what 
 gives \pkg{lme4} package but with \verb{p}-values calculated from 
 \verb{F} statistics of types I - III hypotheses. There are two options for 
 denominator degrees of freedom of \verb{F} statistics: "Satterthwaite" and 
 "Kenward-Roger". The calculation of anova with Kenward-Roger's approximation 
 is based on function from \pkg{pbkrtest} package, the calculation of 
 Satterthwaite's approximation is based on SAS proc mixed theory (see reference).
 The package also provides \code{summary} function, 
 which gives the same as \pkg{lme4} package summary function but with 
 \verb{p}-values and degrees of freedom added for the \verb{t}-test 
 (based on Satterthwaite approximation for denominator degrees of freedom).
The tests on random effects are performed using likelihood ratio tests.
}

\details{
\tabular{ll}{
Package: \tab lmerTest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-01-10\cr
License: \tab GPL\cr
}

The calculation of statistics for the fixed part was developed according to SAS Proc Mixed Theory (see reference).   


}
\author{
Alexandra Kuznetsova <alku@dtu.dk>, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen 
}
\references{

Alexandra Kuznetsova, Per B. Brockhoff and Rune H.B. Christensen (2017)
lmerTest Package: Tests in Linear Mixed Effects Models.
\emph{Journal of Statistical Software}, 82(13), 1--26. doi:10.18637/jss.v082.i13


SAS Technical Report R-101 1978 Tests of Hypotheses in Fixed-Effects Linear Models \emph{Copyright (C)}  (SAS Institute Inc., Cary, NC, USA)

Goodnight, J.H. 1976 General Linear Models Procedure (S.A.S. Institute, Inc.)

Schaalje G.B., McBride J.B., Fellingham G.W. 2002 Adequacy of approximations to distributions of test Statistics in complex mixed linear models 
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%%~~ directory ~~
%%\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
\examples{
## import lmerTest package
library(lmerTest)

## an object of class merModLmerTest
m <- lmer(Informed.liking ~ Gender+Information+Product +(1|Consumer), data=ham)

## gives summary of lmer object. The same as of class merMod but with 
## additional p-values calculated based on Satterthwate's approximations
summary(m)

## anova table the same as of class merMod but with additional F statistics and 
##  and denominator degrees of freedom and 
## p-values calculated based on Satterthwaite's approximations
anova(m)

## anova table the same as of class merMod but with additional F statistics and  
## denominator degrees of freedom and 
## p-values calculated based on Kenward-Roger's approximations
\dontrun{
if(requireNamespace("pbkrtest", quietly = TRUE))
anova(m, ddf = "Kenward-Roger")
}

## anova table of class merMod 
anova(m, ddf="lme4")

## backward elimination of non-significant effects of model m
st <- step(m)

plot(st)

}

