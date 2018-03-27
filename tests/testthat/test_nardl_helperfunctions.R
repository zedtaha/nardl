context("Helper functions for nardl")

# Test cases
mat <- matrix(1:10, ncol = 2)
vec <- 1:6
mat2 <- mat
colnames(mat2) <- c("A","B")

test_that("lagm works correctly",{
  # Test whether the lag works correctly
  expect_equal(lagm(mat,2),
               matrix(c(NA,1,2,3,4,NA,6,7,8,9,
                        NA,NA,1,2,3,NA,NA,6,7,8),
                      ncol = 4,dimnames = list(NULL,c("1_1","2_1","1_2","2_2"))))
  # Test whether error is given
  expect_error(lagm(vec,2))
  expect_equal(lagm(mat2,2),
               matrix(c(NA,1,2,3,4,NA,6,7,8,9,
                        NA,NA,1,2,3,NA,NA,6,7,8),
                      ncol = 4,dimnames = list(NULL,c("A_1","B_1","A_2","B_2"))))
})
