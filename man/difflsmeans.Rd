\name{difflsmeans}
\alias{difflsmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculates Differences of Least Squares Means and Confidence Intervals for the factors of a fixed part of mixed effects model of lmer object. 
}
\description{
 Produces a data frame which resembles to what SAS software gives in proc mixed statement. The approximation for degrees of freedom is Satterthwaite's.
}

\usage{
difflsmeans(model, test.effs=NULL, ddf="Satterthwaite",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
   linear mixed effects model (lmer object).
}

\item{test.effs}{
     charachter vector specifying  names of terms to be tested. If NULL all the terms are tested.
}

\item{ddf}{By default the Satterthwaite's approximation to degrees of freedom is calculated. If ddf="Kenward-Roger", then the Kenward-Roger's approximation is calculated using \code{KRmodcomp} function from \pkg{pbkrtest} package. If ddf="lme4" then the anova table that comes from \pkg{lme4} package is returned.}

  \item{\dots}{
     other potential arguments.
}
}

\value{
 Produces Differences of Least Squares Means (population means) table with p-values and Confidence intervals. 
 }
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lsmeansLT}}, \code{\link{step}}, \code{\link{rand}}
}
\examples{

## import lme4 package and lmerTest package
library(lmerTest)

## specify lmer model
m1 <- lmer(Informed.liking ~ Gender*Information +(1|Consumer), data=ham)

## calculate least squares means for interaction Gender:Information
difflsmeans(m1, test.effs="Gender:Information")
m <- lmer(Coloursaturation ~ TVset*Picture + (1|Assessor), data=TVbo)
plot(difflsmeans(m, test.effs="TVset"))

}
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
