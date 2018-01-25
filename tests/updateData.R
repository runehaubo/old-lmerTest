library(lmerTest)

if(require(nlme))
  data(Oats, package = "nlme")
Oats <- transform(Oats, nitfac = factor(nitro))
Oats$sample1 <- Oats$yield
Oats.lmer <- lmer(yield ~ Variety + nitfac + (1|Block/Variety), data = Oats)
Oats.lmer2 <- lmer(sample1 ~ Variety + nitfac + (1|Block/Variety), data = Oats)

anova(Oats.lmer2)

### Fit the same model to a random subset:
set.seed(12345)
ss <- sample(1:72, 51)
randSub.lmer = update(Oats.lmer, subset = ss)
randSub.lmer <- as(randSub.lmer, "merModLmerTest")
randSub.lmer2 = update(Oats.lmer, subset = ss)
randSub.lmer2 <- as(randSub.lmer2, "merModLmerTest")

lsm1 <- lmerTest::lsmeansLT(randSub.lmer, "Variety")
lsm2 <- lmerTest::lsmeansLT(randSub.lmer, "Variety")

stopifnot(all.equal(lsm1$lsmeans.table, lsm2$lsmeans.table))

an1 <- anova(randSub.lmer)
an2 <- anova(randSub.lmer)

stopifnot(all.equal(an1, an2))

## test for Variety
##        (Intercept) Variety Variety nitfac nitfac nitfac
##[1,]           0       1       0      0      0      0
##[2,]           0       0       1      0      0      0


L <- matrix(0, nrow=2, ncol=6)
L[1,2] <- L[2,3] <-  1
stopifnot(all.equal(unlist(calcSatterth(Oats.lmer, L)), 
          unlist(calcSatterth(Oats.lmer, L))))

stopifnot(all.equal(unlist(calcSatterth(randSub.lmer, L)), 
                    unlist(calcSatterth(randSub.lmer, L))))


## however using sample together with transformation of the variables in formula
## results in error 
Oats.lmer3 <- lmer(log(yield) ~ Variety + nitfac + (1|Block/Variety), data = Oats)
randSub.lmer3 <- as(update(Oats.lmer3, subset = ss), "merModLmerTest")
an <- anova(randSub.lmer3)
stopifnot(ncol(an) == 6L)
set.seed(12345)
randSub.lmer3 <- as(update(Oats.lmer3, subset = sample(1:72, 51)), "merModLmerTest")
an <- anova(randSub.lmer3)
stopifnot(ncol(an) == 4)

