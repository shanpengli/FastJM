
# FastJM

The `FastJM` package implement efficient computation of semi-parametric
joint model of longitudinal and competing risks data.

# Installation

``` r
library(devtools)
devtools::install_github("shanpengli/FastJM", build_vignettes = FALSE)
```

# Example

The `FastJM` package comes with several simulated datasets.

``` r
library(FastJM)
data(ydata)
data(cdata)
fit <- jmcs(ydata = ydata, cdata = cdata, 
long.formula = response ~ time + x1, 
surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
random =  ~ time| ID)
fit
```
