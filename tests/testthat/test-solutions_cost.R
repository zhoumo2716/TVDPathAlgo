library(TVDPathAlgo)
test_that("solutions_cost handles invalid alpha input", {
  # Simulate a simple signal
  y <- c(1, 2, 3, 4)
  
  # Generate a result using path_algo
  result <- path_algo(y, lambda_0 = 100, method = "tri")
  
  # Test for invalid alpha (less than 0)
  expect_error(
    solutions_cost(result, y, alpha = -0.1),
    regexp = "Error: 'alpha' must be a numeric value between 0 and 1"
  )
  
  # Test for invalid alpha (greater than 1)
  expect_error(
    solutions_cost(result, y, alpha = 1.1),
    regexp = "Error: 'alpha' must be a numeric value between 0 and 1"
  )
  
  # Test for non-numeric alpha
  expect_error(
    solutions_cost(result, y, alpha = "invalid"),
    regexp = "Error: 'alpha' must be a numeric value between 0 and 1"
  )
  
  # Test for missing alpha
  expect_silent(
    solutions_cost(result, y)  # Should default to alpha = 0.5
  )
  
  # Test for edge cases: alpha = 0 and alpha = 1
  expect_silent(
    solutions_cost(result, y, alpha = 0)
  )
  expect_silent(
    solutions_cost(result, y, alpha = 1)
  )
})
