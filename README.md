
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
#> Loading required package: statmod
#> Loading required package: MASS
require(survival)
#> Loading required package: survival
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
                       u = seq(3, 4.8, by = 0.2), 
                       method = "GH",
                       obs.time = "time")
survfit
#> 
#> Prediction of Conditional Probabilities of Event
#> based on the pseudo-adaptive Guass-Hermite quadrature rule with 6 quadrature points
#> $`218`
#>       times       CIF1      CIF2
#> 1  2.441634 0.00000000 0.0000000
#> 2  3.000000 0.09629588 0.1110072
#> 3  3.200000 0.11862304 0.1369133
#> 4  3.400000 0.15142590 0.1679708
#> 5  3.600000 0.18413127 0.1839693
#> 6  3.800000 0.21269800 0.2096528
#> 7  4.000000 0.23043413 0.2249182
#> 8  4.200000 0.25459317 0.2500146
#> 9  4.400000 0.25811390 0.2599361
#> 10 4.600000 0.28856883 0.2896654
#> 11 4.800000 0.30829095 0.3134531
#> 
#> $`419`
#>       times       CIF1       CIF2
#> 1  2.432155 0.00000000 0.00000000
#> 2  3.000000 0.02972511 0.02073398
#> 3  3.200000 0.03757608 0.02601222
#> 4  3.400000 0.05003929 0.03270990
#> 5  3.600000 0.06332292 0.03635232
#> 6  3.800000 0.07563241 0.04273814
#> 7  4.000000 0.08376596 0.04677029
#> 8  4.200000 0.09564633 0.05378957
#> 9  4.400000 0.09743720 0.05674168
#> 10 4.600000 0.11449841 0.06602758
#> 11 4.800000 0.12639379 0.07432217
```

To assess the prediction accuracy of the fitted joint model, we may run
`PEjmcs` to calculate the Brier score.

``` r
## evaluate prediction accuracy of fitted joint model using cross-validated Brier Score
PE <- PEjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
             obs.time = "time", method = "GH", 
             quadpoint = NULL, maxiter = 1000, n.cv = 3, 
             survinitial = TRUE)
#> The 1 th validation is done!
#> The 2 th validation is done!
#> The 3 th validation is done!
summary(PE, error = "Brier")
#> 
#> Expected Brier Score at the landmark time of 3 
#> based on 3 fold cross validation
#>   Horizon Time Brier Score 1 Brier Score 2
#> 1          3.6    0.06214519    0.03483533
#> 2          4.0    0.09404153    0.05420762
#> 3          4.4    0.11043488    0.06861870
```

An alternative to assess the prediction accuracy is to run `MAEQjmcs` to
calculate the prediction error by comparing the predicted and empirical
risks stratified on different risk groups based on quantile of the
predicted risks.

``` r
## evaluate prediction accuracy of fitted joint model using cross-validated mean absolute prediction error
MAEQ <- MAEQjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4), 
                 obs.time = "time", method = "GH", 
                 quadpoint = NULL, maxiter = 1000, n.cv = 3, 
                 survinitial = TRUE)
#> The 1 th validation is done!
#> The 2 th validation is done!
#> The 3 th validation is done!
summary(MAEQ, digits = 3)
#> 
#> Sum of absolute error across quintiles of predicted risk scores at the landmark time of 3 
#> based on 3 fold cross validation
#>   Horizon Time  CIF1  CIF2
#> 1          3.6 0.249 0.079
#> 2          4.0 0.292 0.106
#> 3          4.4 0.293 0.147
```

We may also calculate the area under the ROC curve (AUC) to assess the
discrimination measure of joint models.

``` r
## evaluate prediction accuracy of fitted joint model using cross-validated mean AUC
AUC <- AUCjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4),
               obs.time = "time", method = "GH",
               quadpoint = NULL, maxiter = 1000, n.cv = 3)
#> The 1 th validation is done!
#> The 2 th validation is done!
#> The 3 th validation is done!
summary(AUC, digits = 3)
#> 
#> Expected AUC at the landmark time of 3 
#> based on 3 fold cross validation
#>   Horizon Time      AUC1      AUC2
#> 1          3.6 0.7172788 0.6561398
#> 2          4.0 0.6914785 0.6538589
#> 3          4.4 0.6966017 0.6933577
```
