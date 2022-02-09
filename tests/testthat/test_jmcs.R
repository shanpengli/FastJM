myEps <- if (capabilities("long.double")) .Machine$double.eps else 1e-9

test_that(" joint model (jmcs)",
          {
            fit <- jmcs(ydata = ydata, cdata = cdata, 
                        long.formula = response ~ time + gender + x1 + race, 
                        surv.formula = Surv(surv, failure_type) ~ x1 + gender + x2 + race, 
                        random =  ~ time| ID)
            
            expect_equal(mean(fit$beta), 0.29819263, tolerance = (10 ^ 8) * myEps)
            expect_equal(mean(fit$vcov), 0.000757106178, tolerance = (10 ^ 7) * myEps)
            expect_equal(fit$loglike, -8989.389, tolerance = (10 ^ 12) * myEps)
            
 
          })

