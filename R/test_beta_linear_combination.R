#' Test Linear Interpolation for Intermediate Lambda
#'
#' This function tests whether the denoised solution (`beta`) for an intermediate 
#' lambda value can be accurately constructed as a weighted linear combination 
#' of the solutions at the bounding lambda values. The test compares the 
#' objective function outcomes for the interpolated and actual solutions.
#'
#' @param result A list returned by the `path_algo` function, containing 
#'               lambda values, beta solutions, and knot locations etc.
#' @param lambda_test A numeric value specifying the intermediate lambda to test.
#' @param y A numeric vector representing the original noisy signal.
#' 
#' @return A logical value indicating whether the test passed. Returns `TRUE` 
#'         if the objective values of the interpolated and actual solutions 
#'         differ by less than a specified tolerance.
#' 
#' @examples
#' # Example usage
#' y <- c(1, 3, 4, 6, 8, 6)
#' result <- path_algo(y, method = "tri")
#' lambda_test <- (result$lambda_values[2] + result$lambda_values[3]) / 2
#' test_beta_linear_combination(result, lambda_test, y)
#'
#' @export
test_beta_linear_combination <- function(result, lambda_test, y) {
  # Extract lambda values and beta solutions
  lambda_values <- result$lambda_values
  beta_values <- result$beta_values
  
  # Find the bounding lambdas
  k <- max(which(lambda_values >= lambda_test))
  lambda_k <- lambda_values[k]
  lambda_k1 <- lambda_values[k + 1]
  
  beta_k <- beta_values[[k]]
  beta_k1 <- beta_values[[k + 1]]
  
  # Compute weights
  w_k <- (lambda_k1 - lambda_test) / (lambda_k1 - lambda_k)
  w_k1 <- (lambda_test - lambda_k) / (lambda_k1 - lambda_k)
  
  # Construct beta_lambda
  beta_lambda <- w_k * beta_k + w_k1 * beta_k1
  
  # Compute the objective for beta_lambda
  compute_objective <- function(beta, y, lambda) {
    fidelity <- sum((y - beta)^2) / 2
    regularization <- lambda * sum(abs(diff(beta)))
    return(fidelity + regularization)
  }
  
  obj_interpolated <- compute_objective(beta_lambda, y, lambda_test)
  
  # Compute the actual beta_lambda by solving TV denoising
  beta_actual <- tvdenoising::tvdenoising(y, lambda_test)
  obj_actual <- compute_objective(beta_actual, y, lambda_test)
  
  # Compare the objectives
  cat("Objective for interpolated beta:", obj_interpolated, "\n")
  cat("Objective for actual beta:", obj_actual, "\n")
  
  # Return comparison
  return(abs(obj_interpolated - obj_actual) < 1e-6)  # Tolerance for floating-point comparison
}
