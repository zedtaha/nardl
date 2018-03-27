context("Calculation of pssbounds")

data(fod)
reg <- nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = FALSE,case=3)

# JM: This one is copied from the output using dput()
# I assume it is correct, but there's no way of being certain about it.
pssoutput <- structure(list(k = 1L, obs = 49L, fstat = 3.99055616887072, tstat = NULL,
                            case = 3, ftest.I0.p10 = 4.19, ftest.I1.p10 = 4.94,
                            ftest.I0.p05 = 5.22, ftest.I1.p05 = 6.07, ftest.I0.p01 = 7.56,
                            ftest.I1.p01 = 8.685),
                       .Names = c("k", "obs", "fstat","tstat", "case","ftest.I0.p10",
                                  "ftest.I1.p10","ftest.I0.p05", "ftest.I1.p05",
                                  "ftest.I0.p01", "ftest.I1.p01"))

test_that("Case III more than 45 obs gives correct results",{

  # This is a placeholder just to help
  expect_equal(pssbounds(case = reg$case,
                         obs = reg$obs,
                         fstat = reg$fstat,
                         k = reg$k),
               pssoutput)
})
