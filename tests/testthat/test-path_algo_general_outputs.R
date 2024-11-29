library(TVDPathAlgo)
test_that("path_algo handles invalid inputs and errors correctly", {
  # Test 1: y is empty
  expect_error(
    path_algo(y = numeric(0)),
    "Signal 'y' is required but missing.",
    info = "Expected an error when 'y' is empty."
  )
  
  # Test 2: y is a string (invalid type)
  expect_error(
    path_algo(y = "invalid_input"),
    info = "Expected an error when 'y' is a string."
  )
  
  # Test 3: lambda_0 is a string (invalid type)
  expect_error(
    path_algo(y = c(1, 3, 4, 6, 8, 6), lambda_0 = "invalid_input"),
    info = "Expected an error when 'lambda_0' is a string."
  )
  
  # Test 4: method is invalid
  expect_error(
    path_algo(y = c(1, 3, 4, 6, 8, 6), method = "invalid_method"),
    info = "Expected an error when 'method' is invalid."
  )
  
  # Test 5: valid input should return a valid output
  result <- path_algo(y = c(1, 3, 4, 6, 8, 6), lambda_0 = 100, method = "tri")
  expect_type(result, "list")
  expect_type(result$lambda_index, "double")
  expect_type(result$lambda_values, "double")
  expect_type(result$beta_values, "list")
  expect_type(result$knot_locations, "list")
  expect_type(result$loop_number, "double")
})
