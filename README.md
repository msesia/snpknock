SNPknock for R
==========================

This package implements efficient algorithms for generating knockoff copies of variables distributed as discrete Markov chains and hidden Markov models, with specific support for genetic data. To learn more, visit https://web.stanford.edu/group/candes/knockoffs/software/snpknock/.

# Installation

## Stable version

The stable version of this package is hosted on [CRAN](https://cran.r-project.org/package=SNPknock). 

To install this package from CRAN, run the following command in your R console:
```r
install.packages("SNPknock")
```

## Development version

You can install the lastest development version by cloning this repository and building the package from source. Alternatively, you can install it directly from your R console using the [devtools](https://CRAN.R-project.org/package=devtools) package.

To install this package with devtools, run the following command in your R console:

```r
library(devtools)
install_bitbucket("msesia/snpknock/R/SNPknock")
```

If you also want install the vignettes along with the package, type instead:

```r
install_bitbucket("msesia/snpknock/R/SNPknock", build_vignette = TRUE)
```

Note that building the vignettes may require additional R packages.

## Resources
For more information and tutorials, visit
https://web.stanford.edu/group/candes/knockoffs/software/snpknock/

## News

To read about the latest changes, visit the [NEWS page](SNPknock/NEWS).

## License

This software is distributed under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) and it comes with ABSOLUTELY NO WARRANTY.