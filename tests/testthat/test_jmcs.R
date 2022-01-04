cat("\nTests for 'joint model (jmcs)'")

myEps <- if (capabilities("long.double")) .Machine$double.eps else 1e-9

test_that(" joint model for ydata and cdata ",
          {
            fit <- jmcs(ydata = ydata, cdata = cdata, 
                        long.formula = response ~ time + x1, 
                        surv.formula = Surv(surv, failure_type) ~ x1 + x2, 
                        random =  ~ time| ID)
            expect_equal(mean(fit$beta), 0.50305914, tolerance = (10 ^ 8) * myEps)
            expect_equal(mean(fit$vcov), 0.00111349097, tolerance = (10 ^ 7) * myEps)
            expect_equal(fit$loglike, -8993.048, tolerance = (10 ^ 12) * myEps)
          })