# FastJM 1.4.2

* Fix small bugs within C functions.

# FastJM 1.4.1

* Fix small bugs within ```jmcs()```.

# FastJM 1.4.0

* Add the function ```AUCjmcs()``` area under the ROC curve (AUC) to assess the prediction performance of joint models.

# FastJM 1.3.0

* Correct the implementation of dynamic prediction in ```surviftjmcs()```.
* Remove ```plot.surviftjmcs()``` due to theoretical problem.
* Provide two metrics of prediction accuracy of joint model by adding ```PEjmcs()``` and ```MAEQjmcs()```.

# FastJM 1.2.0

* Correct the implementation of dynamic prediction in ```surviftjmcs()``` for the competing risk and add ```summary()``` for providing parameter estimates and SE for both sub-models.

# FastJM 1.1.3

* Correct syntax error on ```jmcs()```.

# FastJM 1.1.2

* Correct testthat error.

# FastJM 1.1.1

* Correct testthat error.

# FastJM 1.1.0

* Provide support for handling categorical variables in both sub-models.

* Provide the ```anova()``` function to compare two fitted joint models.

* Add the ```simulate``` argument in the ```survfitjmcs()``` function to obtain the conditional probabilities using the Gauss-Hermite quadrature rule for numerical integration.

* Adjust the label position of y axis for clarity purposes when ```include.y = TRUE``` in the ```plot.survfitjmcs()``` function.

# FastJM 1.0.1

* Removed unused variables in C functions.

# FastJM 1.0.0

* First version of software with subsequent minor patches. No NEWS file was maintained.
