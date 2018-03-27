context("Calculation of pssbounds")

data(fod)
reg <- nardl(food~inf,p=4,q=4,fod,ic="aic",maxlags = FALSE,graph = FALSE,case=3)

test_that("Case III more than 45 obs gives correct results",{

  # This is a placeholder just to help
  expect_equal(pssbounds(case = reg$case,
                         obs = reg$obs,
                         fstat = reg$fstat,
                         k = reg$k),
    1)
})
