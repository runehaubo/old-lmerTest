require(lmerTest)

load(system.file("testdata","MAMexampledata_red.RData", package="lmerTest"))

lm2 <- lmer(att1 ~ Product + x + Assessor:x + (1|Assessor) + (1|Assessor:Product), 
            data = MAMexampledata_red)
an1 <- anova(lm2, type=1)
stopifnot(dim(an1)==c(2,6))

## check that results agree with lme4
an.lme4 <- anova(lm2, ddf="lme4")
stopifnot(all.equal(an.lme4[,4], an1[,5]))
