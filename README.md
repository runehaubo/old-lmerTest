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

## Discovered a bug?

Please raise a new issue! Preferably add code that illustrates the problem using one of the datasets from the **lmerTest**.

## Installation

Basically there are two options for installing **lmerTest**:

1. Released (stable version) from CRAN: in **R** run `install.packages("lmerTest")`.
2. Development version from GitHub: First load the **devtools** package (and install it if you do not have it) and install the default (master) branch:
```
library("devtools")
install_github("runehaubo/lmerTest")
```
If you haven't already installed a previous version of **lmerTest** you need to also install dependencies (other packages that **lmerTest** depends on and requires you to install to function properly). We recommend that you install **lmerTest** from CRAN (using `install.packages("lmerTest")`) before installing from GitHub as described above. 

An alternative is to use 
```
library("devtools")
install_github("runehaubo/lmerTest", dependencies=TRUE)
```
but that requires you to install all dependent packages from source (which only works if you have the correct compilers installed and set up correctly); installing the pre-compiled packages from CRAN is usually easier.


  

