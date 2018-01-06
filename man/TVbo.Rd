\name{TVbo}
\alias{TVbo}
\title{
  TV dataset
}
\description{
  The \code{TVbo} dataset comes from Bang and Olufsen company. The main purpose was to test products, specified by two attributes Picture and TVset. 15 different response variables (characteristics of the product) were assessed by trained panel list.
}
\usage{
TVbo
}
\format{
\describe{
  \item{\code{Assessor}}{
    factor: numbering identifying assessors
  }
  \item{\code{TVset}}{
    factor: attribute of the product
  }
  \item{\code{Picture}}{
    factor: attribute of the product
  }
  \item{\code{15 Characteristics of the product}}{
    numeric variables: Coloursaturation, Colourbalance, Noise, Depth, Sharpness, Lightlevel, Contrast, Sharpnessofmovement, Flickeringstationary, Flickeringmovement,
 Distortion, Dimglasseffect, Cutting, Flossyedges, Elasticeffect   
  }
  
}}
\source{
  Bang and Olufsen company
}

\examples{


## import lme4 package and lmerTest package
library(lmerTest)

\dontshow{
m <- lmer(log(Coloursaturation)~TVset*Picture+
(1|Assessor)+(1|Assessor:TVset), data=TVbo)
step(m, test.effs="TVset", reduce.fixed=TRUE, reduce.random=TRUE)
}

\dontrun{
m <- lmer(Coloursaturation ~ TVset*Picture+
(1|Assessor)+(1|Assessor:TVset), data=TVbo)

step(m, test.effs="TVset", reduce.fixed=FALSE, reduce.random=TRUE)
}
}

\keyword{datasets}
