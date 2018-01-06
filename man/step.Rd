\name{step}
\alias{step}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs backward elimination of non-significant effects of linear mixed effects model:
}
\description{
performs automatic backward elimination of all effects of linear mixed effect model. First backward elimination of the random part is performed following by backward elimination of the fixed part. Finally LSMEANS (population means) and differences of LSMEANS for the fixed part of the model are calculated and the final model is provided. The \verb{p}-values for the fixed effects are calculated from \verb{F} test based on Sattethwaite's or Kenward-Roger approximation), \verb{p}-values for the random effects are based on likelihood ratio test.  All analysis may be performed on \code{lmer} object of \pkg{lme4} package. 
}
\usage{
step(model, ddf = "Satterthwaite", type = 3, alpha.random = 0.1, alpha.fixed = 0.05, 
reduce.fixed = TRUE, reduce.random = TRUE, fixed.calc = TRUE, lsmeans.calc = TRUE,
difflsmeans.calc = TRUE, test.effs = NULL,  keep.effs = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{model}{
    
     linear mixed effects model (lmer object).
}

\item{ddf}{
    
     approximation for denominator degrees of freedom. By default Satterthwaite's 
     approximation. ddf="Kenward-Roger"" calculates Kenward-Roger approximation
}

\item{type}{
   type of hypothesis to be tested (SAS notation). Either type=1 or type=3. 
}

\item{alpha.random}{
    
     significance level for elimination of the random part (for LRT test)
}

\item{alpha.fixed}{
     significance level for elimination of the fixed part (for \code{F} test and \code{t}-test for least squares means)
}

\item{reduce.fixed}{
    
    logical for whether the reduction of the fixed part is required
}

\item{reduce.random}{
    
    logical for whether the reduction of the random part is required
}

\item{fixed.calc}{
    
    logical for whether the calculation of the table for fixed effects is needed. 
    If FALSE then only the analysis of random effects is done
}

\item{lsmeans.calc}{    
    logical for whether the calculation of LSMEANS(population means) is required
}

\item{difflsmeans.calc}{    
    logical for whether the calculation of differences of LSMEANS is required
}

\item{test.effs}{
    
    charachter vector specifying the names of terms to be tested in LSMEANS. 
    If NULL all the terms are tested. If \code{lsmeans.calc==FALSE} then 
    LSMEANS are not calculated.
}

\item{keep.effs}{
    
    charachter vector specifying the names of terms to be kept in the model
    even if being non-significant
}


  \item{\dots}{
     other potential arguments.
}
}
\details{
Elimination of all effects is done one at a time.
Elimination of the fixed part is done by the principle of marginality that is: the highest order interactions are tested first: if they are significant, the lower order effects are not tested for significance. 
The step function of lmerTest overrides the one from stats package for lm objects. So if the
lmerTest is attached and one wants to call step fof lm object, then needs to use stats::step
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{rand.table}{data frame with value of Chi square statistics, p-values  for the likelihood ratio test for random effects}
  \item{anova.table }{data frame with tests for whether the model fixed terms are significant (Analysis of Variance)}
  \item{lsmeans.table}{Least Squares Means data frame with p-values and Confidence intervals}
  \item{diffs.lsmeans.table}{Differences of Least Squares Means data frame with p-values and Confidence intervals}
  \item{model}{Final model - object of merLmerTest(contains mer class) or gls (after all the required reduction has been performed)}
  }
%% ...


\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\note{
  For the random coefficient models: in the random part if correlation is present between slope and intercept, then the simplified model will contain just an intercept. That is if the random part of the initial model is \code{(1+c|f)}, then this model is compared to \code{(1|f)} by using LRT. If there are multiple slopes, then the the slope with the highest \verb{p}-value (and higher then alpha level) is eliminated. That is if the random part of the initial model has the following form \code{(1+c1+c2|f)}, then two simplified models are constracted and compared to the initial one: the first one has \code{(1+c1|f)} in the random part and the second one has: \code{(1+c2|f)}. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rand}}, \code{\link{lsmeans}}, \code{\link{difflsmeans}}
}
\examples{
#import lme4 package and lmerTest package
library(lmerTest)

\dontrun{
m <- lmer(Informed.liking ~ Product*Information*Gender+ 
(1|Consumer) + (1|Product:Consumer), data=ham)


#elimination of non-significant effects
s <- step(m)

#plot of post-hoc analysis of the final model
plot(s)

m <- lmer(Coloursaturation ~ TVset*Picture+
(1|Assessor)+(1|Assessor:TVset), data=TVbo)

step(m, keep.effs = "Assessor")

}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}