\name{merModLmerTest-class}
\title{Mixed Model Representations}
\docType{class}
% Classes
\alias{merModLmerTest-class}
%\alias{anova.merModLmerTest}
%\alias{coef.mer-class}
%\alias{summary.merModLmerTest}
%\alias{ranef.mer-class}
% Methods
%\alias{anova,merLmerTest-method}

%\alias{summary,merLmerTest-method}
%\alias{summary,summary.merLmerTest-method}

\description{The \code{merModLmerTest} \emph{contains} \code{merMod} class of \pkg{lme4} package and overloads \code{anova} and \code{summary} functions.  
}

\section{Objects from the Class}{Objects can be created via the
  \code{\link{lmer}}  functions.
} 
 
%\section{Methods}{
%anova \code{signature(object = "merLmerTest", \dots)}: returns the data frame with F test for each term in formula and corresponding p-values. For multiple arguments, model comparison statistics. type=c(1,3) (type=3 is the default one) constructs corresponding hypothesis F test and calculates  p-values, sum of squares and means squares. The definition of the types of hypotheses comes from SAS theory proc mixed. The approximations for degrees of freedom by default are Satterthwaite's (ddf="Satterthwaite"). If ddf="Kenward-Roger", then the Kenward-Roger's approximations are calculated using \code{KRmodcomp} function from \code{pbkrtest} package. If ddf="lme4" then the anova table that comes from \code{lme4} package is returned. 
%method.grad  is an argument for the approximation method for the grad function, which is used in calculation of denominator degrees of freedom. method.grad could be "simple" or "Richardson". method.grad="simple" is the default one and the fastest one, method.grad="Richardson" gives more accurate results.
%      \code{\link{anova}}
%}
 
\seealso{
  \code{\link{lmer}()}
}
\examples{

(m <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))

## type 3 anova table with  denominator degrees of freedom
## calculated based on Satterthwaite's approximation
anova(m)

## type 1 anova table with denominator degrees of freedom
## calculated based on Satterthwaite's approximation
\dontrun{
anova(m, type=1)
}

## type3 anova table with additional F statistics and denominator degrees of freedom
## calculated based on Kenward-Roger's approximation
\dontrun{
if(requireNamespace("pbkrtest", quietly = TRUE))
anova(m, ddf="Kenward-Roger")
}

## anova table, that is returned by lme4 package
anova(m, ddf="lme4")

## summary of merModLmerTest object. Returns the same as merMod object but with an 
## additional column of p values for the t test. 
summary(m)

}
\keyword{classes}
