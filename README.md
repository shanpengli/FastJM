---
output:
  pdf_document: default
  html_document: default
---

# FastJM

<!-- badges: start -->

[![R-CMD-check](https://github.com/shanpengli/FastJM/workflows/R-CMD-check/badge.svg)](https://github.com/shanpengli/FastJM/actions)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/FastJM)](https://cran.r-project.org/package=FastJM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/FastJM)](https://cran.r-project.org/package=FastJM)
[![CRAN_time_from_release](https://www.r-pkg.org/badges/ago/FastJM)](https://cran.r-project.org/package=FastJM)
[![CRAN_Status_Badge_version_last_release](https://www.r-pkg.org/badges/version-last-release/FastJM)](https://cran.r-project.org/package=FastJM)
<!-- badges: end -->

The `FastJM` package implement efficient computation of semi-parametric
joint model of longitudinal and competing risks data.

# Example

The `FastJM` package comes with several simulated datasets. To fit a
joint model, we use `jmcs` function.

``` r
require(FastJM)
#> Loading required package: FastJM
#> Loading required package: MASS
#> Loading required package: statmod
data(ydata)
data(cdata)
fit <- jmcs(ydata = ydata, cdata = cdata, 
            long.formula = response ~ time + gender + x1 + race, 
            surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
            random =  ~ time| ID)
fit
#> 
#> Call:
#>  jmcs(ydata = ydata, cdata = cdata, long.formula = response ~ time + gender + x1 + race, random = ~time | ID, surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race) 
#> 
#> Data Summary:
#> Number of observations: 3067 
#> Number of groups: 1000 
#> 
#> Proportion of competing risks: 
#> Risk 1 : 34.9 %
#> Risk 2 : 29.8 %
#> 
#> Numerical intergration:
#> Method: pseudo-adaptive Guass-Hermite quadrature
#> Number of quadrature points:  6 
#> 
#> Model Type: joint modeling of longitudinal continuous and competing risks data 
#> 
#> Model summary:
#> Longitudinal process: linear mixed effects model
#> Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard
#> 
#> Loglikelihood:  -8989.389 
#> 
#> Fixed effects in the longitudinal sub-model:  response ~ time + gender + x1 + race 
#> 
#>             Estimate      SE   Z value  p-val
#> (Intercept)  2.01853 0.05704  35.38803 0.0000
#> time         0.98292 0.03147  31.22885 0.0000
#> genderMale  -0.07766 0.05860  -1.32527 0.1851
#> x1          -1.47810 0.05851 -25.26356 0.0000
#> raceWhite    0.04527 0.05911   0.76581 0.4438
#> 
#>         Estimate      SE  Z value  p-val
#> sigma^2  0.49182 0.01793 27.43751 0.0000
#> 
#> Fixed effects in the survival sub-model:  Surv(surv, failure_type) ~ x1 + gender + x2 + race 
#> 
#>              Estimate      SE   Z value  p-val
#> x1_1          0.54672 0.18540   2.94892 0.0032
#> genderMale_1 -0.18781 0.11935  -1.57359 0.1156
#> x2_1         -1.10450 0.12731  -8.67602 0.0000
#> raceWhite_1  -0.10027 0.11802  -0.84960 0.3955
#> x1_2          0.62986 0.20064   3.13927 0.0017
#> genderMale_2  0.10834 0.13065   0.82926 0.4070
#> x2_2         -1.76738 0.15245 -11.59296 0.0000
#> raceWhite_2   0.03194 0.13049   0.24479 0.8066
#> 
#> Association parameters:                 
#>               Estimate      SE Z value  p-val
#> (Intercept)_1  0.93973 0.12160 7.72809 0.0000
#> time_1         0.31691 0.19318 1.64051 0.1009
#> (Intercept)_2  0.96486 0.13646 7.07090 0.0000
#> time_2         0.03772 0.24137 0.15629 0.8758
#> 
#> 
#> Random effects:                 
#>   Formula: ~time | ID 
#>                  Estimate      SE  Z value  p-val
#> (Intercept)       0.52981 0.03933 13.47048 0.0000
#> time              0.25885 0.02262 11.44217 0.0000
#> (Intercept):time -0.02765 0.02529 -1.09330 0.2743
```

The `FastJM` package can make dynamic prediction given the longitudinal
history information. Below is a toy example for competing risks data.
Conditional cumulative incidence probabilities for each failure will be
presented.

``` r
ND <- ydata[ydata$ID %in% c(419, 218), ]
ID <- unique(ND$ID)
NDc <- cdata[cdata$ID  %in% ID, ]
survfit <- survfitjmcs(fit, 
                       ynewdata = ND, 
                       cnewdata = NDc, 
                       u = seq(3, 4.5, by = 0.3), 
                       M = 50,
                       seed = 100)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |=========                                                             |  12%  |                                                                              |==========                                                            |  14%  |                                                                              |===========                                                           |  16%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |================                                                      |  22%  |                                                                              |=================                                                     |  24%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  41%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  45%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  49%  |                                                                              |====================================                                  |  51%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  55%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  59%  |                                                                              |===========================================                           |  61%  |                                                                              |============================================                          |  63%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |==================================================                    |  71%  |                                                                              |===================================================                   |  73%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |===================================================================   |  96%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
survfit
#> 
#> Prediction of Conditional Probabilities of Event
#>  based on 50 Monte Carlo samples
#> 
#> $`218`
#> $`218`$`Cumulative incidence probabilities for type 1 failure`
#>      times      Mean     Median 95%HDLower 95%HDUpper
#> 1 2.441634 0.0000000 0.00000000 0.00000000  0.0000000
#> 2 3.000000 0.1019023 0.09173377 0.02800572  0.2912309
#> 3 3.300000 0.1411337 0.12478303 0.03925318  0.4000652
#> 4 3.600000 0.1951804 0.17273625 0.05517985  0.5465473
#> 5 3.900000 0.2331410 0.20945648 0.06683557  0.6477954
#> 6 4.200000 0.2699894 0.24178899 0.07841441  0.7432997
#> 7 4.500000 0.2918208 0.26053237 0.08534956  0.7977818
#> 
#> $`218`$`Cumulative incidence probabilities for type 2 failure`
#>      times      Mean    Median 95%HDLower 95%HDUpper
#> 1 2.441634 0.0000000 0.0000000 0.00000000  0.0000000
#> 2 3.000000 0.1140959 0.1047456 0.03061918  0.3670499
#> 3 3.300000 0.1665314 0.1549190 0.04647752  0.5120441
#> 4 3.600000 0.1875582 0.1749764 0.05323437  0.5656966
#> 5 3.900000 0.2234971 0.2056436 0.06531201  0.6510426
#> 6 4.200000 0.2534424 0.2310376 0.07590219  0.7156340
#> 7 4.500000 0.2929031 0.2666981 0.09082805  0.7912904
#> 
#> 
#> $`419`
#> $`419`$`Cumulative incidence probabilities for type 1 failure`
#>      times       Mean     Median 95%HDLower 95%HDUpper
#> 1 2.432155 0.00000000 0.00000000 0.00000000 0.00000000
#> 2 3.000000 0.03286548 0.03286250 0.01565922 0.07106902
#> 3 3.300000 0.04747870 0.04762011 0.02278781 0.10154752
#> 4 3.600000 0.06975108 0.07027051 0.03384960 0.14672429
#> 5 3.900000 0.08671654 0.08750940 0.04242203 0.18023420
#> 6 4.200000 0.10497785 0.10599381 0.05181927 0.21525971
#> 7 4.500000 0.11708209 0.11824975 0.05816036 0.23778849
#> 
#> $`419`$`Cumulative incidence probabilities for type 2 failure`
#>      times       Mean     Median  95%HDLower 95%HDUpper
#> 1 2.432155 0.00000000 0.00000000 0.000000000 0.00000000
#> 2 3.000000 0.02269061 0.02021804 0.005899132 0.05233959
#> 3 3.300000 0.03448109 0.03077516 0.009046389 0.07867464
#> 4 3.600000 0.03965429 0.03541735 0.010447487 0.09002123
#> 5 3.900000 0.04945079 0.04421991 0.013145948 0.11104256
#> 6 4.200000 0.05844794 0.05232337 0.015665138 0.12994575
#> 7 4.500000 0.07159601 0.06444799 0.019413217 0.15694409
```

Plot the cumulative incidence function for each failure with the
historical longitudinal observations.

``` r
oldpar <- par(mfrow = c(2, 2), mar = c(5, 4, 4, 4))
plot(survfit, estimator = "both", include.y = TRUE)
```

![](man/figures/README-unnamed-chunk-4-1.png)<!-- -->

``` r
par(oldpar)
```
