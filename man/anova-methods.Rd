\name{anova-methods}
\docType{methods}
\alias{anova,ANY-method}
\alias{anova,merModLmerTest-method}
\alias{anova.merModLmerTest}
\title{Methods for function \code{anova} in package \pkg{lmerTest}}
\description{
Methods for Function \code{anova} in Package \pkg{lmerTest}
}


\usage{\S4method{anova}{merModLmerTest}(object, \dots , ddf="Satterthwaite", 
type=3)}
\arguments{
\item{object}{object of class "merModLmerTest"}
\item{\dots}{object of class "merModLmerTest". Then the model comparison statistisc will be calculated}
\item{ddf}{By default the Satterthwaite's approximation to degrees of freedom is calculated. If ddf="Kenward-Roger", then the Kenward-Roger's approximation is calculated using \code{KRmodcomp} function from \pkg{pbkrtest} package. If ddf="lme4" then the anova table that comes from \pkg{lme4} package is returned.}
\item{type}{type of hypothesis to be tested. Could be type=3 or type=2 or type = 1 
(The definition comes from SAS theory) }
}


%}}

\details{
 According to (Goodnight, J.H. 1976) the behaviour of the type 3 hypothesis is not fully studied for the situations with missing cells (where observations are missing at some factor-level combination). A warning is returned in such cases.
 }

\references{
SAS Technical Report R-101 1978 Tests of Hypotheses in Fixed-Effects Linear Models \emph{Copyright (C)}  (SAS Institute Inc., Cary, NC, USA)

Goodnight, J.H. 1976 General Linear Models Procedure (S.A.S. Institute, Inc.)

Schaalje G.B., McBride J.B., Fellingham G.W. 2002 Adequacy of approximations to distributions of test Statistics in complex mixed linear models 
}

\keyword{methods}
\examples{
#import lmerTest package
library(lmerTest)

m.ham <- lmer(Informed.liking ~ Product*Information*Gender 
+ (1|Consumer), data = ham)

# type 3 anova table with  denominator degrees of freedom
# calculated based on Satterthwaite's approximation
anova(m.ham)

# type 1 anova table with denominator degrees of freedom
# calculated based on Satterthwaite's approximation
\dontrun{
anova(m.ham, type = 1)
}

# type3 anova table with additional F statistics and denominator degrees of freedom
# calculated based on Kenward-Roger's approximation
if(require(pbkrtest))
anova(m.ham, ddf = "Kenward-Roger")

\dontrun{
# anova table, that is returned by lme4 package
anova(m.ham, ddf = "lme4")
}

\dontshow{
\dontrun{
fm2 <- lmer(Preference ~ sens2 + I(sens1^2)  +
                           (1+sens2|Consumer), data=carrots)
anova(fm2)
}
}
}

