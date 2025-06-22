# FastJM

<!-- badges: start -->

[![R-CMD-check](https://github.com/shanpengli/FastJM/workflows/R-CMD-check/badge.svg)](https://github.com/shanpengli/FastJM/actions)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/FastJM)](https://cran.r-project.org/package=FastJM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/FastJM)](https://cran.r-project.org/package=FastJM)
[![CRAN\_time\_from\_release](https://www.r-pkg.org/badges/ago/FastJM)](https://cran.r-project.org/package=FastJM)
[![CRAN\_Status\_Badge\_version\_last\_release](https://www.r-pkg.org/badges/version-last-release/FastJM)](https://cran.r-project.org/package=FastJM)
[![R-CMD-check](https://github.com/shanpengli/FastJM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shanpengli/FastJM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `FastJM` package implement efficient computation of semi-parametric
joint model of longitudinal and competing risks data.

# Example

The `FastJM` package comes with several simulated datasets. To fit a
joint model, we use `jmcs` function.

    require(FastJM)
    #> Loading required package: FastJM
    #> Loading required package: survival
    #> Loading required package: MASS
    #> Loading required package: statmod
    #> Loading required package: magrittr
    require(survival)
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

The `FastJM` package can make dynamic prediction given the longitudinal
history information. Below is a toy example for competing risks data.
Conditional cumulative incidence probabilities for each failure will be
presented.

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

To assess the prediction accuracy of the fitted joint model, we may run
`PEjmcs` to calculate the Brier score.

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
    #> 1          3.6    0.05888837    0.03483090
    #> 2          4.0    0.08889966    0.05428274
    #> 3          4.4    0.10517866    0.06865793

An alternative to assess the prediction accuracy is to run `MAEQjmcs` to
calculate the prediction error by comparing the predicted and empirical
risks stratified on different risk groups based on quantile of the
predicted risks.

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
    #> 1          3.6 0.083 0.120
    #> 2          4.0 0.178 0.161
    #> 3          4.4 0.191 0.152

We may also calculate the area under the ROC curve (AUC) to assess the
discrimination measure of joint models.

    ## evaluate prediction accuracy of fitted joint model using cross-validated mean AUC
    AUC <- AUCjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4),
                   obs.time = "time", method = "GH",
                   quadpoint = NULL, maxiter = 1000, n.cv = 3, metric = "AUC")
    #> The 1 th validation is done!
    #> The 2 th validation is done!
    #> The 3 th validation is done!
    summary(AUC, digits = 3)
    #> 
    #> Expected AUC at the landmark time of 3 
    #> based on 3 fold cross validation
    #>   Horizon Time      AUC1      AUC2
    #> 1          3.6 0.7366710 0.7097309
    #> 2          4.0 0.7154871 0.6760296
    #> 3          4.4 0.7336741 0.7254964

Alternatively, we can also calculate concordance index (Cindex) as
another discrimination measure.

    ## evaluate prediction accuracy of fitted joint model using cross-validated mean Cindex
    Cindex <- AUCjmcs(fit, seed = 100, landmark.time = 3, horizon.time = c(3.6, 4, 4.4),
                   obs.time = "time", method = "GH",
                   maxiter = 1000, n.cv = 3, metric = "Cindex")
    #> The 1 th validation is done!
    #> The 2 th validation is done!
    #> The 3 th validation is done!
    summary(Cindex, digits = 3)
    #> 
    #> Expected Cindex at the landmark time of 3 
    #> based on 3 fold cross validation
    #>   Horizon Time   Cindex1   Cindex2
    #> 1          3.6 0.6864341 0.6772933
    #> 2          4.0 0.6859882 0.6765425
    #> 3          4.4 0.6862253 0.6757857

To fit a joint model with multiple longitudinal outcomes and competing
risks, we can use the `mvjmcs` function.

    data(mvydata)
    data(mvcdata)
    library(FastJM)
    mvfit <- mvjmcs(ydata = mvydata, cdata = mvcdata,
                  long.formula = list(Y1 ~ X11 + X12 + time,
                                      Y2 ~ X11 + X12 + time),
                  random = list(~ time | ID,
                                ~ 1 | ID),
                  surv.formula = Surv(survtime, cmprsk) ~ X21 + X22, 
                  maxiter = 1000, opt = "optim",
                  tol = 1e-3, print.para = FALSE)
    #> Time difference of 30.46499 secs
    mvfit
    #> 
    #> Call:
    #>  mvjmcs(ydata = mvydata, cdata = mvcdata, long.formula = list(Y1 ~ X11 + X12 + time, Y2 ~ X11 + X12 + time), random = list(~time | ID, ~1 | ID), surv.formula = Surv(survtime, cmprsk) ~ X21 + X22, maxiter = 1000, opt = "optim", tol = 0.001, print.para = FALSE) 
    #> 
    #> Data Summary:
    #> Number of observations: 5645 
    #> Number of groups: 800 
    #> 
    #> Proportion of competing risks: 
    #> Risk 1 : 41.62 %
    #> Risk 2 : 11.25 %
    #> 
    #> Model Type: joint modeling of multivariate longitudinal continuous and competing risks data 
    #> 
    #> Model summary:
    #> Longitudinal process: linear mixed effects model
    #> Event process: cause-specific Cox proportional hazard model with non-parametric baseline hazard
    #> 
    #> Fixed effects in the longitudinal sub-model:  list(Y1 ~ X11 + X12 + time, Y2 ~ X11 + X12 + time) 
    #> 
    #>                  Estimate      SE   Z value  p-val
    #> (Intercept)_bio1  4.97563 0.05380  92.48724 0.0000
    #> X11_bio1          1.46466 0.08028  18.24488 0.0000
    #> X12_bio1          1.99751 0.01426 140.09647 0.0000
    #> time_bio1         0.84039 0.03949  21.28056 0.0000
    #> (Intercept)_bio2  9.97529 0.04927 202.44462 0.0000
    #> X11_bio2          0.97973 0.07333  13.36122 0.0000
    #> X12_bio2          2.00943 0.01309 153.45426 0.0000
    #> time_bio2         0.99380 0.00455 218.61091 0.0000
    #> 
    #>              Estimate      SE    Z value  p-val
    #> sigma^2_bio1  0.49303 0.00018 2716.52548 0.0000
    #> sigma^2_bio2  0.49758 0.00965   51.55103 0.0000
    #> 
    #> Fixed effects in the survival sub-model:  Surv(survtime, cmprsk) ~ X21 + X22 
    #> 
    #>       Estimate      SE  Z value  p-val
    #> X21_1  0.93331 0.13443  6.94271 0.0000
    #> X22_1  0.51069 0.03163 16.14809 0.0000
    #> X21_2 -0.21818 0.24899 -0.87628 0.3809
    #> X22_2  0.48438 0.05920  8.18250 0.0000
    #> 
    #> Association parameters:                 
    #>                   Estimate      SE  Z value  p-val
    #> (Intercept)_1bio1  0.49907 0.07535  6.62377 0.0000
    #> time_1bio1         0.70541 0.08500  8.29933 0.0000
    #> (Intercept)_1bio2 -0.54622 0.07970 -6.85387 0.0000
    #> (Intercept)_2bio1  0.63192 0.13341  4.73654 0.0000
    #> time_2bio1         0.66064 0.16716  3.95222 0.0001
    #> (Intercept)_2bio2 -0.48370 0.15878 -3.04643 0.0023
    #> 
    #> 
    #> Random effects:                 
    #>   bio 1 :  ~time | ID 
    #>   bio 2 :  ~1 | ID 
    #>                       Estimate      SE  Z value  p-val
    #> Intercept1             1.02121 0.06469 15.78559 0.0000
    #> time1                  0.91534 0.05830 15.70184 0.0000
    #> Intercept2             0.88204 0.05324 16.56757 0.0000
    #> Intercept1:time1      -0.09356 0.04533 -2.06390 0.0390
    #> Intercept1:Intercept2  0.04353 0.04051  1.07442 0.2826
    #> time1:Intercept2      -0.06540 0.04225 -1.54769 0.1217

We can extract the components of the model as follows:

    # Longitudinal fixed effects
    fixef(mvfit, process = "Longitudinal")
    #> (Intercept)_bio1         X11_bio1         X12_bio1        time_bio1 
    #>        4.9756250        1.4646558        1.9975109        0.8403906 
    #> (Intercept)_bio2         X11_bio2         X12_bio2        time_bio2 
    #>        9.9752870        0.9797270        2.0094259        0.9938004
    summary(mvfit, process = "Longitudinal")
    #>       Longitudinal   coef     SE 95%Lower 95%Upper p-values
    #> 1 (Intercept)_bio1 4.9756 0.0538   4.8702   5.0811        0
    #> 2         X11_bio1 1.4647 0.0803   1.3073   1.6220        0
    #> 3         X12_bio1 1.9975 0.0143   1.9696   2.0255        0
    #> 4        time_bio1 0.8404 0.0395   0.7630   0.9178        0
    #> 5 (Intercept)_bio2 9.9753 0.0493   9.8787  10.0719        0
    #> 6         X11_bio2 0.9797 0.0733   0.8360   1.1234        0
    #> 7         X12_bio2 2.0094 0.0131   1.9838   2.0351        0
    #> 8        time_bio2 0.9938 0.0045   0.9849   1.0027        0

    # Survival fixed effects
    fixef(mvfit, process = "Event")
    #> $Risk1
    #>     X21_1     X22_1 
    #> 0.9333053 0.5106853 
    #> 
    #> $Risk2
    #>      X21_2      X22_2 
    #> -0.2181847  0.4843826
    summary(mvfit, process = "Event")
    #>             Survival    coef exp(coef) SE(coef) 95%Lower 95%Upper
    #> 1              X21_1  0.9333    2.5429   0.1344   0.6698   1.1968
    #> 2              X22_1  0.5107    1.6664   0.0316   0.4487   0.5727
    #> 3              X21_2 -0.2182    0.8040   0.2490  -0.7062   0.2698
    #> 4              X22_2  0.4844    1.6232   0.0592   0.3684   0.6004
    #> 5  (Intercept)_1bio1  0.4991    1.6472   0.0753   0.3514   0.6468
    #> 6         time_1bio1  0.7054    2.0247   0.0850   0.5388   0.8720
    #> 7  (Intercept)_1bio2 -0.5462    0.5791   0.0797  -0.7024  -0.3900
    #> 8  (Intercept)_2bio1  0.6319    1.8812   0.1334   0.3704   0.8934
    #> 9         time_2bio1  0.6606    1.9360   0.1672   0.3330   0.9883
    #> 10 (Intercept)_2bio2 -0.4837    0.6165   0.1588  -0.7949  -0.1725
    #>    95%exp(Lower) 95%exp(Upper) p-values
    #> 1         1.9539        3.3095   0.0000
    #> 2         1.5663        1.7730   0.0000
    #> 3         0.4935        1.3097   0.3809
    #> 4         1.4454        1.8229   0.0000
    #> 5         1.4211        1.9093   0.0000
    #> 6         1.7140        2.3917   0.0000
    #> 7         0.4954        0.6770   0.0000
    #> 8         1.4484        2.4435   0.0000
    #> 9         1.3952        2.6866   0.0001
    #> 10        0.4516        0.8416   0.0023

    # Random effects for first few subjects
    head(ranef(mvfit))
    #>   (Intercept)_bio1  time_bio1 (Intercept)_bio2
    #> 1        1.2362134 -0.5279470      -1.20320837
    #> 2       -0.5287808 -0.3316951       1.56034227
    #> 3       -1.1591939  0.3289143       0.17095050
    #> 4       -1.4259973 -1.9371873      -0.09556157
    #> 5        0.2387970 -1.9513868       0.02260521
    #> 6       -0.1206852 -0.0225784       0.06433657

Currently, prediction and validation features (e.g., survfitjmcs,
PEjmcs, AUCjmcs) are implemented for models of class jmcs. Extension to
mvjmcs is under active development and will be available later this
year.
