\name{ham}
\alias{ham}
\title{
 Conjoint study of dry cured ham
}
\description{
 One of the purposes of the study was to investigate the effect of information given to the consumers measured in hedonic liking for the hams. Two of the hams were Spanish and two were Norwegian, each origin representing different salt levels and different aging time. The information about origin was given in such way that both true and false information was given. essentially a 4*2 design with 4 samples and 2 information levels. A total of 81 Consumers participated in the study.    
}
\usage{
ham
}
\format{
\describe{
  \item{\code{Consumer}}{
    factor with 81 levels: numbering identifying consumers
  }
  \item{\code{Product}}{
   factor with four levels
  }
  \item{\code{Informed.liking}}{
    numeric: hedonic liking for the products
  }
  \item{\code{Information}}{
    factor with two levels
  }
  \item{\code{Gender}}{
    factor with two levels (gender)
  }
  \item{\code{Age}}{
    numeric: age of Consumer
  } 
}}

\references{
 "Alternative methods for combining design variables and consumer preference with information about attitudes and demographics in conjoint analysis" . T. Naes, V.Lengard, S. Bolling Johansen, M. Hersleth
}
\examples{

#import lmerTest package
library(lmerTest)

m <- lmer(Informed.liking ~ Product*Information*Gender 
+ (1|Product:Consumer) , data=ham)

#anova table with p-values with Satterthwaite's approximation for denominator 
#degrees of freedom
anova(m)

#analysis of random and fixed parts and post hoc 
#analysis of Product and Information effects
step(m, reduce.random=FALSE, reduce.fixed=FALSE,
test.effs=c("Product", "Information"))

}

\keyword{datasets}
