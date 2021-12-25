
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
ND <- ydata[ydata$ID %in% c(419, 218), ]
ID <- unique(ND$ID)
NDc <- cdata[cdata$ID  %in% ID, ]
survfit <- survfitjmcs(fit, 
                       ynewdata = ND, 
                       cnewdata = NDc, 
                       u = seq(3, 4.8, by = 0.2), 
                       M = 100)
survfit

# Prediction of Conditional Probabilities of Event
#   based on 100 Monte Carlo samples
# 
# $`218`
# $`218`$`Cumulative incidence probabilities for type 1 failure`
#       times      Mean     Median   95%Lower  95%Upper
# 1  2.441634 0.0000000 0.00000000 0.00000000 0.0000000
# 2  3.000000 0.1042089 0.08501926 0.03167305 0.2361048
# 3  3.200000 0.1265429 0.10504460 0.03963747 0.2839815
# 4  3.400000 0.1583911 0.13562638 0.05111005 0.3528610
# 5  3.600000 0.1893599 0.16617725 0.06106590 0.4168163
# 6  3.800000 0.2158977 0.19320244 0.06977831 0.4685929
# 7  4.000000 0.2320297 0.21009903 0.07502376 0.4987382
# 8  4.200000 0.2537060 0.23138395 0.08186588 0.5374272
# 9  4.400000 0.2567729 0.23435594 0.08281619 0.5427199
# 10 4.600000 0.2832645 0.26198248 0.09068803 0.5865385
# 11 4.800000 0.3001481 0.28039390 0.09835850 0.6121160
# 
# $`218`$`Cumulative incidence probabilities for type 2 failure`
#       times      Mean    Median   95%Lower  95%Upper
# 1  2.441634 0.0000000 0.0000000 0.00000000 0.0000000
# 2  3.000000 0.1313353 0.1234973 0.04242986 0.3261823
# 3  3.200000 0.1604026 0.1517229 0.05344900 0.3873897
# 4  3.400000 0.1946131 0.1861315 0.06727745 0.4539016
# 5  3.600000 0.2115190 0.2033670 0.07463875 0.4844251
# 6  3.800000 0.2386425 0.2317131 0.08740303 0.5298890
# 7  4.000000 0.2546127 0.2500329 0.09537971 0.5546613
# 8  4.200000 0.2804137 0.2728269 0.10906506 0.5915387
# 9  4.400000 0.2906288 0.2831806 0.11478894 0.6051476
# 10 4.600000 0.3208521 0.3140568 0.13262484 0.6428100
# 11 4.800000 0.3443734 0.3334078 0.14814465 0.6684966
# 
# 
# $`419`
# $`419`$`Cumulative incidence probabilities for type 1 failure`
#       times       Mean     Median    95%Lower   95%Upper
# 1  2.432155 0.00000000 0.00000000 0.000000000 0.00000000
# 2  3.000000 0.03054195 0.02650777 0.009682452 0.06045416
# 3  3.200000 0.03858212 0.03358989 0.012290781 0.07604747
# 4  3.400000 0.05134066 0.04492659 0.016484373 0.10049422
# 5  3.600000 0.06488123 0.05705475 0.021002405 0.12609668
# 6  3.800000 0.07742440 0.06837507 0.025250365 0.14949952
# 7  4.000000 0.08571303 0.07591435 0.028094204 0.16477277
# 8  4.200000 0.09778120 0.08697629 0.032288587 0.18673217
# 9  4.400000 0.09958676 0.08864021 0.032921754 0.18998849
# 10 4.600000 0.11682448 0.10441903 0.039049267 0.22063142
# 11 4.800000 0.12878350 0.11543248 0.043377777 0.24151576
# 
# $`419`$`Cumulative incidence probabilities for type 2 failure`
#       times       Mean     Median    95%Lower   95%Upper
# 1  2.432155 0.00000000 0.00000000 0.000000000 0.00000000
# 2  3.000000 0.02109040 0.01930663 0.006954115 0.05474365
# 3  3.200000 0.02644037 0.02425306 0.008762312 0.06825689
# 4  3.400000 0.03321015 0.03053788 0.011074073 0.08515274
# 5  3.600000 0.03685292 0.03393781 0.012335097 0.09411542
# 6  3.800000 0.04324688 0.03994099 0.014582649 0.10960005
# 7  4.000000 0.04727913 0.04374293 0.016015715 0.11924210
# 8  4.200000 0.05428128 0.05037558 0.018534414 0.13575401
# 9  4.400000 0.05723952 0.05318927 0.019610086 0.14264139
# 10 4.600000 0.06653590 0.06206307 0.023022383 0.16401714
# 11 4.800000 0.07480723 0.07000650 0.026133855 0.18250746
```
