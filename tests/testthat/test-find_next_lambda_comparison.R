test_that("find_next_lambda and find_next_lambda_tri return the same results", {
  # Input signal and parameters
  y <- c(1, 3, 4, 6, 8, 6)
  current_lambda <- 100
  B <- numeric(0)  # Initial boundary set is empty
  S <- integer()   # Initial direction set is empty
  
  # Run find_next_lambda
  result_general <- find_next_lambda(y, current_lambda, B, S)
  
  # Run find_next_lambda_tri
  result_tri <- find_next_lambda_tri(y, current_lambda, B, S)
  
  # Compare results
  expect_equal(result_general$next_lambda, result_tri$next_lambda,
               info = "The next_lambda values do not match.")
  expect_equal(result_general$B_new, result_tri$B_new,
               info = "The B_new values do not match.")
  expect_equal(result_general$S_new, result_tri$S_new,
               info = "The S_new values do not match.")
})
