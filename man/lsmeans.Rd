\name{lsmeans}
\alias{lsmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object. 
}
\description{
 Produces a data frame which resembles to what SAS software gives in proc mixed statement. The approximation of degrees of freedom is Satterthwate's. This is a deprecated function, use lsmeansLT function instead.
}
\usage{
lsmeans(model, test.effs = NULL, ddf="Satterthwaite",  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
   linear mixed effects model (lmer object).
}


\item{test.effs}{
    
     charachter vector specyfying the names of terms to be tested. If NULL all the terms are tested.
}

\item{ddf}{By default the Satterthwaite's approximation to degrees of freedom is calculated. If ddf="Kenward-Roger", then the Kenward-Roger's approximation is calculated using \code{KRmodcomp} function from \pkg{pbkrtest} package. If ddf="lme4" then the anova table that comes from \pkg{lme4} package is returned}

  \item{\dots}{
     other potential arguments.
}
}

\value{
 Produces Least Squares Means (population means) table with p-values and Confidence intervals. 
 }
\references{
 \pkg{doBy} package, \pkg{gplots} package
 }
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\note{
  For construction of the contrast matrix \code{popMatrix} function from \pkg{doBy} package was used. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{step}}, \code{\link{rand}}, \code{\link{lsmeansLT}}, \code{\link{difflsmeans}}
}
\examples{


## import lme4 package and lmerTest package
library(lmerTest)

## specify lmer model
m1 <- lmer(Informed.liking ~ Gender*Information +(1|Consumer), data=ham)

\dontrun{
## calculate least squares means for interaction Gender:Information
lsmeans(m1, test.effs="Gender:Information")
m <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor), data=TVbo)
plot(lsmeans(m))
lsmeans(m, test.effs="TVset")
}
}
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
