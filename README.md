
# FastJM

The `FastJM` package implement efficient computation of semi-parametric
joint model of longitudinal and competing risks data.

# Installation

``` r
library(devtools)
devtools::install_github("shanpengli/FastJM", build_vignettes = FALSE)
```

# Example

The `FastJM` package comes with several simulated datasets. To fit a
joint model, we use `jmcs` function.

``` r
library(FastJM)
data(ydata)
data(cdata)
fit <- jmcs(ydata = ydata, cdata = cdata, 
long.formula = response ~ time + x1, 
surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
random =  ~ time| ID)
fit
# Call:
#  jmcs(ydata = ydata, cdata = cdata, long.formula = response ~ time + x1, random = ~time | ID, surv.formula = Surv(surv, failure_type) ~ x1 + x2) 
# 
# Data Summary:
# Number of observations: 3067 
# Number of groups: 1000 
# 
# Proportion of competing risks: 
# Risk 1 : 34.9 %
# Risk 2 : 29.8 %
# 
# Numerical intergration:
# Method: pseudo-adaptive Guass-Hermite quadrature
# Number of quadrature points:  6 
# 
# Model Type: joint modeling of longitudinal continuous and competing risks data 
# 
# Model summary:
# Longitudinal process: linear mixed effects model
# Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard
# 
# Loglikelihood:  -8993.048 
# 
# Fixed effects in the longitudinal sub-model:  response ~ time + x1 
# 
#             Estimate      SE   Z value  p-val
# (Intercept)  2.00436 0.03977  50.40076 0.0000
# time         0.98386 0.03134  31.39642 0.0000
# x1          -1.47905 0.05829 -25.37387 0.0000
# 
#         Estimate     SE  Z value  p-val
# sigma^2  0.49167 0.0179 27.46542 0.0000
# 
# Fixed effects in the survival sub-model:  Surv(surv, failure_type) ~ x1 + x2 
# 
#      Estimate      SE   Z value  p-val
# x1_1  0.53967 0.18334   2.94347 0.0032
# x2_1 -1.10429 0.12705  -8.69205 0.0000
# x1_2  0.64837 0.20087   3.22789 0.0012
# x2_2 -1.76851 0.15214 -11.62453 0.0000
# 
# Association parameters:                 
#       Estimate      SE Z value  p-val
# nu1_1  0.94239 0.12086 7.79744 0.0000
# nu1_2  0.34631 0.19589 1.76789 0.0771
# nu2_1  0.95740 0.13516 7.08357 0.0000
# nu2_2  0.00769 0.24067 0.03195 0.9745
# 
# 
# Random effects:                 
#   Formula: ~time | ID 
#                  Estimate      SE  Z value  p-val
# (Intercept)       0.53136 0.03930 13.51907 0.0000
# time              0.25915 0.02259 11.47045 0.0000
# (Intercept):time -0.02614 0.02517 -1.03838 0.2991
```

The `FastJM` package can proceed dynamic prediction given the
longitudinal history information.

``` r
ND <- ydata[ydata$ID %in% c(23, 878, 230, 113), ]
ID <- unique(ND$ID)
NDc <- cdata[cdata$ID  %in% ID, ]
survfit <- survfitjmcs(fit, 
                       ynewdata = ND, 
                       cnewdata = NDc, 
                       u = seq(4, 4.9, by = 0.1), 
                       M = 100)
survfit
```
