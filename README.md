# dmhBGM
Bayesian analysis of graphical models of binary and/or ordinal variables using the Double Metropolis Hastings algorithm. 

This is a clone of an earlier version of [bgms](https://github.com/MaartenMarsman/bgms), replacing the pseudolikelihood approach with the Double Metropolis Hastings algorithm. Currently not under active development. If you find a problem with the software, please create an issue on Github.

## Installation

You can install the package from Github using

``` r
if (!requireNamespace("remotes")) { 
  install.packages("remotes")   
}   
remotes::install_github("MaartenMarsman/dmhBGM")
```
