
# FastJM

<!-- badges: start -->

[![R-CMD-check](https://github.com/shanpengli/FastJM/workflows/R-CMD-check/badge.svg)](https://github.com/shanpengli/FastJM/actions)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/FastJM)](https://cran.r-project.org/package=FastJM)
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
#>       Estimate      SE Z value  p-val
#> nu1_1  0.93973 0.12160 7.72809 0.0000
#> nu1_2  0.31691 0.19318 1.64051 0.1009
#> nu2_1  0.96486 0.13646 7.07090 0.0000
#> nu2_2  0.03772 0.24137 0.15629 0.8758
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
                       u = seq(3, 4.8, by = 0.2), 
                       M = 100,
                       seed = 100)
survfit
#> 
#> Prediction of Conditional Probabilities of Event
#>  based on 100 Monte Carlo samples
#> 
#> $`218`
#> $`218`$`Cumulative incidence probabilities for type 1 failure`
#>       times      Mean     Median   95%Lower  95%Upper
#> 1  2.441634 0.0000000 0.00000000 0.00000000 0.0000000
#> 2  3.000000 0.1042121 0.09168331 0.02551633 0.2519007
#> 3  3.200000 0.1274778 0.11325876 0.03194108 0.3045305
#> 4  3.400000 0.1610705 0.14534683 0.04182345 0.3776787
#> 5  3.600000 0.1940238 0.17773682 0.05219704 0.4449708
#> 6  3.800000 0.2222930 0.20647921 0.06168556 0.4988277
#> 7  4.000000 0.2394984 0.22306490 0.06779426 0.5299034
#> 8  4.200000 0.2627392 0.24462767 0.07650370 0.5696356
#> 9  4.400000 0.2660452 0.24767362 0.07778946 0.5750665
#> 10 4.600000 0.2946978 0.27582552 0.09027359 0.6197911
#> 11 4.800000 0.3129948 0.29613833 0.09759391 0.6456466
#> 
#> $`218`$`Cumulative incidence probabilities for type 2 failure`
#>       times      Mean    Median   95%Lower  95%Upper
#> 1  2.441634 0.0000000 0.0000000 0.00000000 0.0000000
#> 2  3.000000 0.1206324 0.1084623 0.02732776 0.2649958
#> 3  3.200000 0.1480588 0.1343882 0.03402085 0.3175604
#> 4  3.400000 0.1806585 0.1661039 0.04227422 0.3763843
#> 5  3.600000 0.1969183 0.1820365 0.04647734 0.4040791
#> 6  3.800000 0.2231869 0.2086400 0.05336726 0.4459997
#> 7  4.000000 0.2387306 0.2246734 0.05755825 0.4693403
#> 8  4.200000 0.2638940 0.2481454 0.06453518 0.5034274
#> 9  4.400000 0.2738487 0.2583460 0.06736726 0.5155370
#> 10 4.600000 0.3033470 0.2899113 0.07611484 0.5553568
#> 11 4.800000 0.3262614 0.3132795 0.08311886 0.5887102
#> 
#> 
#> $`419`
#> $`419`$`Cumulative incidence probabilities for type 1 failure`
#>       times       Mean     Median   95%Lower   95%Upper
#> 1  2.432155 0.00000000 0.00000000 0.00000000 0.00000000
#> 2  3.000000 0.02872658 0.02687972 0.01029811 0.05720026
#> 3  3.200000 0.03630935 0.03400285 0.01308744 0.07214694
#> 4  3.400000 0.04834371 0.04526916 0.01758261 0.09574403
#> 5  3.600000 0.06117241 0.05729828 0.02244669 0.12066025
#> 6  3.800000 0.07305485 0.06846044 0.02701788 0.14350806
#> 7  4.000000 0.08089509 0.07582621 0.03007797 0.15848148
#> 8  4.200000 0.09235147 0.08658972 0.03461434 0.18021130
#> 9  4.400000 0.09407310 0.08820711 0.03530296 0.18346131
#> 10 4.600000 0.11050764 0.10361973 0.04199070 0.21432507
#> 11 4.800000 0.12196150 0.11437911 0.04673742 0.23556138
#> 
#> $`419`$`Cumulative incidence probabilities for type 2 failure`
#>       times       Mean     Median    95%Lower   95%Upper
#> 1  2.432155 0.00000000 0.00000000 0.000000000 0.00000000
#> 2  3.000000 0.02141642 0.01879897 0.006000041 0.05236234
#> 3  3.200000 0.02685920 0.02362739 0.007557645 0.06536595
#> 4  3.400000 0.03376509 0.02978113 0.009552346 0.08169752
#> 5  3.600000 0.03749466 0.03312294 0.010641212 0.09041257
#> 6  3.800000 0.04406357 0.03904481 0.012581513 0.10556092
#> 7  4.000000 0.04821486 0.04280435 0.013819254 0.11503185
#> 8  4.200000 0.05542342 0.04936500 0.015990536 0.13128493
#> 9  4.400000 0.05846205 0.05214284 0.016914275 0.13806252
#> 10 4.600000 0.06800699 0.06090433 0.019842673 0.15912501
#> 11 4.800000 0.07648702 0.06876341 0.022495396 0.17740045
```
