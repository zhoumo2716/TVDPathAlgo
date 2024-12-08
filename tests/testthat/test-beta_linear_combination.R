library(testthat)
library(TVDPathAlgo)

test_that("test_linear_combination validates interpolated beta solution", {
  # Input signal
  y <- c(1, 3, 4, 6, 8, 6)
  
  # Run path_algo to generate results
  result <- path_algo(y, method = "tri")
  
  # Choose an intermediate lambda value
  lambda_test <- (result$lambda_values[2] + result$lambda_values[3]) / 2
  
  # Test linear combination
  valid_interpolation <- test_beta_linear_combination(result, lambda_test, y)
  
  # Check if the interpolation test passes
  expect_true(valid_interpolation, "Linear combination test for beta failed!")
})
