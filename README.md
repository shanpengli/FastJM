# FastJM
Fast fitting of joint models with longitudinal and competing risks data applying pseudo-adaptive quadrature rules and one-way scan

# Install from github:
library(devtools)

install_github("shanpengli/FastJM", build_vignettes = TRUE)

Reminder: GSL library must be preinstalled before the installation of FastJM. See INSTALL file.

# Run the package using a toy example with filepath
library(FastJM)

set.seed(123)

yfile=system.file("extdata", "simy0.txt", package = "FastJM")

cfile=system.file("extdata", "simc0.txt", package = "FastJM")

mfile=system.file("extdata", "simm0.txt", package = "FastJM")

res2=jmcs(p1=3,yfile,cfile,mfile,point=6,type_file=TRUE, do.trace = F)

res2

# Run the package using a toy example with data frames

## This package depends on R (>=3.5.0)

library(FastJM)

set.seed(123)

data(mydata)

ydata <- mydata$ydata

cdata <- mydata$cdata

mdata <- mydata$mdata

fit <- jmcs(p1 = 3, yfile = ydata, cfile = cdata, mfile = mdata, type_file = F, point = 6, do.trace = F)

