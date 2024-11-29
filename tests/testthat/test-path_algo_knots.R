test_that("Knots for lambda_test are consistent with path_algo output", {
  # Input signal
  y <- c(1, 3, 4, 6, 8, 6)
  
  # Run path_algo
  lambda_0 <- 100
  result <- path_algo(y, lambda_0 = lambda_0)
  
  # Extract values from result
  algorithm_lambda_values <- result$lambda_values
  algorithm_knot_locations <- result$knot_locations
  
  # Define test lambdas
  lambda_test <- c(50, 4, 2.8, 2.2, 1, 0.8, 0.58, 0.3)
  
  for (i in seq_along(lambda_test)) {
    # Compute beta for test lambda
    beta_test <- tvdenoising(y, lambda_test[i])
    
    # Compute knots for test lambda
    change_points_test <- which(diff(beta_test) != 0)
    change_locations_test <- if (length(change_points_test) == 0) 0 else change_points_test + 1
    
    # Find the nearest lambda
    nearest_index <- which.min(abs(algorithm_lambda_values - lambda_test[i]))
    nearest_lambda <- algorithm_lambda_values[nearest_index]
    
    # Determine the second lambda based on nearest lambda
    if (nearest_lambda > lambda_test[i]) {
      second_index <- max(which(algorithm_lambda_values <= lambda_test[i]))
    } else {
      second_index <- min(which(algorithm_lambda_values >= lambda_test[i]))
    }
    
    # Retrieve knots for nearest and second nearest lambdas
    nearest_knots <- algorithm_knot_locations[[nearest_index]]
    second_knots <- algorithm_knot_locations[[second_index]]
    
    # Test if all knots are within the combined set of nearest and second knots
    expect_true(
      all(change_locations_test %in% c(nearest_knots, second_knots)),
      info = paste("Mismatch detected for lambda_test =", lambda_test[i])
    )
  }
})
