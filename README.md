# lmerTest - Tests in Linear Mixed Effects Models

## Main features

- Implements Satterthwaite's degrees-of-freedom for _t_ and _F_ tests of fixed-effects in linear mixed models
- Provides _p_-values for `summary` _t_-tests and `anova` _F_-tests when for `lmer`-objects (class `merMod`) from the **lme4** package to provide _p_-values
- Provides ANOVA Type III _F_-tests (type I and II are also available)
- Provides an anova-like table for random-effects based on likelihood-ratio tests
- A `step` function implements automatic model selection

## Citation

To cite lmerTest in publications use:

Kuznetsova A., Brockhoff P.B. and Christensen R.H.B. (2017). "lmerTest Package: Tests in Linear Mixed Effects Models." Journal of Statistical Software, 82(13), pp. 1–26. doi: 10.18637/jss.v082.i13.

Corresponding BibTeX entry:

  @Article{,
    title = {{lmerTest} Package: Tests in Linear Mixed Effects Models},
    author = {Alexandra Kuznetsova and Per B. Brockhoff and Rune H. B.
      Christensen},
    journal = {Journal of Statistical Software},
    year = {2017},
    volume = {82},
    number = {13},
    pages = {1--26},
    doi = {10.18637/jss.v082.i13},
  }


