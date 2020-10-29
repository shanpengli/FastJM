# FastJM
Fast fitting of joint models with longitudinal and competing risks data applying pseudo-adaptive quadrature rules and one-way scan

# Install from github:
library(devtools)

install_github("shanpengli/FastJM", build_vignettes = TRUE)

Reminder: GSL library must be preinstalled before the installation of FastJM. See INSTALL file.

## This package depends on R (>=3.5.0)

# Run the package using a toy example 
data(ydata)

data(cdata)

fit <- jmcs(ydata = ydata, cdata = cdata,

                    long.formula = response ~ time + x1,
                    
                    surv.formula = Surv(surv, failure_type) ~ x1 + x2,
                    
                    ID = "ID",
                    
                    RE = "time",
                    
                    model = "interslope")

fit


