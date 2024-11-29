#' Path Algorithm for Regularization
#'
#' This function implements a path-following algorithm for solving a regularization problem.
#' It iteratively finds the next lambda value, computes the solution, and tracks knot locations.
#'
#' @param y A numeric vector representing the input signal.
#' @param lambda_0 A numeric value specifying the initial lambda value. Defaults to \eqn{\max(y) \times 100}.
#' @param method A character string specifying the method. Must be `"tri"` (tri-diagonal) or `"general"`.
#'        Defaults to `"tri"`.
#'
#' @return A list containing:
#' \describe{
#'   \item{lambda_index}{A numeric vector of indices corresponding to lambda values.}
#'   \item{lambda_values}{A numeric vector of lambda values at each iteration.}
#'   \item{beta_values}{A list of vectors, each containing the beta solution for a lambda value.}
#'   \item{knot_locations}{A list of vectors, each containing the knot locations for a lambda value.}
#'   \item{loop}{An integer indicating the total number of iterations.}
#' }
#'
#' @importFrom tvdenoising tvdenoising
#' 
#' @examples
#' # Example usage with dummy data
#' y <- c(1, 2, 3, 4)
#' result <- path_algo(y, lambda_0 = 100, method = "tri")
#' print(result)
#'
#' @export
path_algo <- function(y, lambda_0, method) {
  if (missing(y) || length(y) == 0) {
    stop("Error: Signal 'y' is required but missing or empty.")
  }
  
  if (!is.numeric(y)) {
    stop("Error: Signal 'y' must be a numeric vector.")
  }
  
  if (missing(lambda_0)) {
    lambda_0 <- max(y) * 100  
  }
  
  if (!is.numeric(lambda_0) || length(lambda_0) != 1 || lambda_0 < 0) {
      stop("Error: 'lambda_0' must be a single non-negative numeric value.")
  }
  
  if (missing(method)) {
    method <- "tri"
  } else {
    if (!(method %in% c("tri", "general"))) {
      stop("Error: Argument 'method' is not valid. Please input 'tri' or 'general'.")
    }
  }
  
  # Initialize vectors and lists to store results
  lambda_index <- c(1)
  lambda_values <- c(lambda_0)
  beta_values <- list(tvdenoising::tvdenoising(y, lambda_0))  # Use fully qualified name for clarity
  knot_locations <- list(0)  # Use a list for knot locations
  
  B <- numeric(0)  # Initial boundary set is empty
  S <- integer()  # Initial direction set is empty
  
  loop <- 1
  current_lambda <- lambda_0
  while (current_lambda > 0) {
    # Find the next lambda and update B and S using the selected method
    if (method == "tri") {
      result <- find_next_lambda_tri(y, current_lambda, B, S)
    } else if (method == "general") {
      result <- find_next_lambda(y, current_lambda, B, S)
    }
    
    next_lambda <- result$next_lambda
    B <- result$B_new
    S <- result$S_new
    
    # Solve for beta
    beta <- tvdenoising::tvdenoising(y, next_lambda)
    beta_values[[loop + 1]] <- beta  # Append to list
    
    # Record the knot locations
    change_point <- which(diff(beta) != 0)
    knot_location <- change_point + 1
    knot_locations[[loop + 1]] <- knot_location  # Append to list
    
    # Record the new lambda and update the loop
    lambda_values <- c(lambda_values, next_lambda)
    current_lambda <- next_lambda
    lambda_index <- c(lambda_index, length(lambda_index) + 1)
    loop <- loop + 1
  }
  
  return(list(
    lambda_index = lambda_index,
    lambda_values = lambda_values,
    beta_values = beta_values,
    knot_locations = knot_locations,
    loop_number = loop
  ))
}
