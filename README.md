# FastJM
Fast fitting of joint models with longitudinal and competing risks data applying pseudo-adaptive quadrature rules and one-way scan

# Install from github:
library(devtools)

install_github("shanpengli/FastJM")

# Run the package using a toy example
set.seed(123)

yfile=system.file("extdata", "simy0.txt", package = "FastJM")

cfile=system.file("extdata", "simc0.txt", package = "FastJM")

mfile=system.file("extdata", "simm0.txt", package = "FastJM")

res2=jmcs(p1=3,yfile,cfile,mfile,point=6,type_file=TRUE, do.trace = F)

res2
