\name{calcSatterth}
\alias{calcSatterth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 F-test based on the Satterthwaite's approximation for denominator degrees of freedom.
}
\description{
 Produces a list with the values for an approximate F-test based on the Satterthwaite's approximation.
}

\usage{
calcSatterth(model, L)
}

\arguments{
    \item{model}{
     linear mixed effects model (lmer object).
     }
     
     \item{L}{
     hypothesis contrast matrix or a vector
     }
    
    \item{\dots}{
       other potential arguments.}
}

\value{
A list with the results from the F test
 \item{denom}{numeric. Denominator degrees of freedom, calculated with the Satterthwaite's approximation}
 \item{Fstat}{numeric. F statistic}
  \item{pvalue}{numeric. p-value of the corresponding F test}
  \item{ndf}{numeric. Numerator degrees of freedom}
}

\details{
 F test for the null hypothesis \eqn{H_0: L \beta
  = 0}, where  \eqn{\beta} is a vector of the same length as
\code{fixef(model)}
 }

\references{
Schaalje G.B., McBride J.B., Fellingham G.W. 2002 Adequacy of approximations to distributions of test Statistics in complex mixed linear models 
}

\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}


\seealso{
\code{\link{anova}}
}
\examples{
## import lme4 package and lmerTest package
library(lmerTest)

## specify lmer model for the sleepstudy data from the lme4 package
m <- lmer(Reaction ~ Days + (1 + Days|Subject), sleepstudy)


L <- cbind(0,1) ## specify contrast vector
calcSatterth(m, L) ## calculate F test

## specify model for the ham data
m.ham <- lmer(Informed.liking ~ Product + (1|Consumer), data = ham)

## specify contrast vector for testing product effect
L <- matrix(0, ncol = 4, nrow = 3)
L[1,2] <- L[2,3] <- L[3,4] <- 1
calcSatterth(m.ham, L)

## by using anova function we get the same result
anova(m.ham)
}
\keyword{function}

