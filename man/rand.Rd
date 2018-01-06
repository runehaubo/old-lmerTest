\name{rand}
\alias{rand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs likelihood ratio test on random effects of linear mixed effects model. 
}
\description{
Returns a data frame with  values of Chi square statistics and corresponding 
p-values of likelihood ratio tests. 
}
\usage{
rand(model, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{model}{
    
     linear mixed effects model (lmer object).
}


  \item{\dots}{
     other potential arguments.
}
}

\value{
Produces a data frame with LR tests for the random terms.
}

\details{
 The columns of the data are:
 
 
 Chisq: The value of the chi square statistics
 
 Chi Df: The degrees of freedom for the test
 
 p.value: The p-value of the likelihood ratio test for the effect

}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{step}}, \code{\link{lsmeans}}, \code{\link{difflsmeans}}
}
\examples{

#import lme4 package and lmerTest package
library(lmerTest)

#lmer model with correlation between intercept and slopes
#in the random part
m <- lmer(Preference ~ sens2+Homesize+(1+sens2|Consumer), data=carrots)

# table with  p-values for the random effects
rand(m)
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}