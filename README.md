# lmerTest - Tests in Linear Mixed Effects Models

## Main features

- Provides _p_-values for `summary` _t_-tests and `anova` _F_-tests for `lmer`-objects (class `merMod`) from the **lme4** package
- Degrees-of-freedom for _t_ and _F_-tests are based on Satterthwaite's method
- Provides Type III ANOVA tables (type I and II are also available)
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
