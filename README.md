# ATAComb

Say it like catacomb.

Functions and Rmarkdowns for QC and analysis of scATAC-seq data.

## Requirements

```
Depends:
  GenomicRanges,
  data.table,
  Matrix
Imports:
  assertthat,
  igraph,
  IRanges,
  Matrix,
  parallel,
  RANN
```

## Installation

This package can be installed from Github using the `devtools` package.

You may first need to register your GitHub PAT, as this is a private repository.
```
Sys.setenv(GITHUB_PAT = "your-access-token-here")
devtools::install_github("aifimmunology/ATAComb")
```
## Test Data

TBD

## Tests

Tests for `ATAComb` are implemented using the `testthat` testing suite:  
https://testthat.r-lib.org/

To run tests for this package, download the repository from Github and run `devtools::test()` in R from the base directory of the package.

Extra-stringent, CRAN-level package testing can be performed using `devtools::check()` in R.

## Style and Standards

This package aims to conform to the tidyverse style guide:  
https://style.tidyverse.org/index.html

General information about R package conventions can be found in `R Packages`:  
http://r-pkgs.had.co.nz/


<a id="legal_info"></a>

# Legal Information

<a id="license"></a>

## License

The license for this package is available on Github in the file LICENSE.txt in this repository.

<a id="support"></a>

## Level of Support

We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.

<a id="contrib"></a>

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in the file CONTRIBUTING.md in this repository.

