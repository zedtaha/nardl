context("Helper functions for nardl")

# Test cases
mat <- matrix(1:10, ncol = 2)
vec <- 1:6

test_that("lagm works correctly",{
  # Test whether the lag works correctly
  expect_equal(lagm(mat,2),
               matrix(c(NA,1,2,3,4,NA,6,7,8,9,
                        NA,NA,1,2,3,NA,NA,6,7,8),
                      ncol = 4,dimnames = list(NULL,c("_1","_1","_2","_2"))))
  # Test whether error is given
  expect_error(lagm(vec,2))
})
