# FastJM
Fast fitting of joint models with longitudinal and competing risks data applying pseudo-adaptive quadrature rules and one-way scan

# Install from github:

    library(devtools)
    install_github("shanpengli/FastJMEigen", build_vignettes = FALSE, ref = "master")

## This package depends on R (>=4.1.0)

# Run the package using a toy example 

    data(ydata)
    data(cdata)
    fit <- jmcs(ydata = ydata, cdata = cdata,
                long.formula = response ~ time + x1,
                random = ~time|ID,
                surv.formula = Surv(surv, failure_type) ~ x1 + x2)
    fit

