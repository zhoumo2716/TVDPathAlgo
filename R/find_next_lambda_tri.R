#' Find the Next Lambda Value for Regularization Path using Tri-diagonal matrix solver
#'
#' This function computes the next lambda value and updates the boundary set (`B`) 
#' and direction set (`S`) for a given regularization path problem. The method 
#' involves solving an optimization problem using difference matrices and boundary conditions.
#' This differs from find_next_lambda function in terms that it uses Tri-diagonal matrix solver for between steps
#'
#' @param y A numeric vector representing the signal data.
#' @param current_lambda A numeric value indicating the current lambda value in the path.
#' @param B A numeric vector of indices representing the current boundary set. Can be empty.
#' @param S A numeric vector of direction values corresponding to the boundary set. Can be empty.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{next_lambda}{The next lambda value along the regularization path. Returns `0` if no valid lambda is found.}
#'   \item{B_new}{An updated boundary set with newly added indices.}
#'   \item{S_new}{An updated direction set with values corresponding to the updated boundary set.}
#' }
#'
#' @examples
#' # Example signal data
#' y <- c(1, 2, 3, 4, 5)
#' current_lambda <- 10
#' B <- c()
#' S <- c()
#'
#' # Find the next lambda, updated B and S
#' result <- find_next_lambda(y, current_lambda, B, S)
#' print(result)
#'
#' @export



#find_next_lambda function v2 Tri-diagonal matrix solver 
find_next_lambda_tri <- function(y, current_lambda, B, S) {
  
  # Define the difference matrix D for the y-signal
  n <- length(y)
  D <- matrix(0, n - 1, n)
  original_indices <- 1:(n - 1)
  remaining_indices <- setdiff(original_indices, B)
  
  # Exit condition: Check if all original indices are already in B
  if (length(remaining_indices) == 0) {
    #cat("All indices are already in B. Exiting...\n")
    return(list(next_lambda = 0, B_new = B, S_new = S))
  }
  
  for (i in 1:(n - 1)) {
    D[i, i] <- -1
    D[i, i + 1] <- 1
  }
  
  # Handle the case when B and S are empty (starting with a high lambda)
  if (is.null(B) || length(B) == 0) {
    # Full difference matrix as no boundary has been hit
    D_minus_B <- D
    M_inv <- tridiagonal_inverse_matrix(D %*% t(D))  # Using the tridiagonal solver
    D_B <- NULL
  } else if (length(remaining_indices) < 2) {  # If D_{-B} * D_{-B}^T is smaller than 2x2
    # Extract submatrices based on the current boundary set B
    D_minus_B <- D[-sort(B), , drop = FALSE]
    M_inv <- solve(D_minus_B %*% t(D_minus_B))  # Directly solve for small matrix
    D_B <- D[sort(B), , drop = FALSE]
  } else {
    # General case with remaining indices > 2
    D_minus_B <- D[-sort(B), , drop = FALSE]
    M_inv <- tridiagonal_inverse_matrix(D_minus_B %*% t(D_minus_B))  # Efficient tridiagonal inverse
    D_B <- D[sort(B), , drop = FALSE]
  }
  
  
  
  # Debugging: Print matrices and current B, S
  #cat("Current lambda:", current_lambda, "\n")
  #cat("Current B:", B, "\n")
  #cat("Current S:", S, "\n")
  #cat("D matrix:\n")
  #cat("D original indices:", original_indices, "\n")
  #print(D)
  #cat("D_{-B} indices:",  remaining_indices, "\n")
  #cat("D_{-B} matrix:\n")
  #print(D_minus_B)
  #if (is.null(D_B) || nrow(D_B) == 0) {
    #cat("D_{B} matrix is empty\n")
  #} else {
    #cat("D_{B} matrix:\n")
    #print(D_B)
  #}
  #cat("Current M_inv: \n")
  #print(M_inv)
  
  
  # Compute a and b
  a <- M_inv %*% D_minus_B %*% y
  if (length(B) > 0) {
    b <- M_inv %*% D_minus_B %*% t(D_B) %*% S
  } else {
    b <- rep(1, length(a))  # Default value if no boundary is active
  }
  
  #cat("Current a: \n")
  #print(a)
  #cat("Current b: \n")
  #print(b)
  #cat("length of a: \n")
  #print(length(a))
  
  # Compute t_i for all coordinates, considering both directions (±1)
  valid_pos <- (b + 1) != 0
  valid_neg <- (b - 1) != 0
  
  # Initialize t_values_pos and t_values_neg to Inf by default
  t_values_pos <- rep(Inf, length(a))
  t_values_neg <- rep(Inf, length(a))
  
  # Populate t_values_pos and t_values_neg with valid computations, set invalid ones to Inf
  t_values_pos[valid_pos] <- a[valid_pos] / (b[valid_pos] + 1)  # Valid pos computations
  t_values_neg[valid_neg] <- a[valid_neg] / (b[valid_neg] - 1)  # Valid neg computations
  
  # Initialize t_values to store the final result
  t_values <- numeric(length(a))
  
  # Populate t_values with the valid values in order
  for (i in 1:length(a)) {
    if (t_values_pos[i] > 0 & t_values_pos[i] < current_lambda) {
      t_values[i] <- t_values_pos[i]
    } else if (t_values_neg[i] > 0 & t_values_neg[i] < current_lambda) {
      t_values[i] <- t_values_neg[i]
    }
  }
  
  
  #cat("valid pos: \n")
  #print(valid_pos)
  #cat("Current t_pos: \n")
  #print(t_values_pos)
  #cat("Current t_neg: \n")
  #print(t_values_neg)
  #cat("Current t: \n")
  #print(t_values)
  
  # Find the maximum positive t_i (next lambda)
  if (length(t_values) == 0) {
    next_lambda <- 0
    new_indices <- integer(0)
  } else {
    # Create a named vector to map t_values to their corresponding original indices
    named_t_values <- setNames(t_values, remaining_indices)  # remaining_indices should be precomputed
    next_lambda <- max(named_t_values[named_t_values < current_lambda])
    #new_indices <- as.integer(names(named_t_values[named_t_values == next_lambda]))
    new_indices <- as.integer(names(named_t_values[abs(named_t_values - next_lambda) < 1e-8]))  # Tolerance for floating-point
  } 
  
  #cat("Named t_values:\n")
  #for (name in names(named_t_values)) {
    #cat("Index:", name, "-> Value:", named_t_values[name], "\n")
  #}
  #cat("will be added indices are: ", new_indices, "\n")
  
  
  # Update boundary set B and direction set S
  if (length(new_indices) > 0) {
    B_new <- c(B, new_indices)  # Append all new indices
    #S_new <- c(S, sign(t_values[new_indices]))  # Append signs for all new indices
    S_new <- c(S, ifelse(named_t_values[as.character(new_indices)] > 0, 1, -1))
  } else {
    B_new <- B  # No update
    S_new <- S  # No update
  }
  
  #cat("t_values:", named_t_values, "\n")
  #cat("new_indices:", new_indices, "\n")
  #cat("t_values[new_indices]:", named_t_values[as.character(new_indices)], "\n")
  
  B_new <- as.vector(B_new)
  S_new <- as.vector(S_new)
  
  #cat("Next lambda:", next_lambda, "B_new:", paste(B_new, collapse = ", "), "S_new:", paste(S_new, collapse = ", "), "\n")
  
  
  # Return results as a list
  return(list(next_lambda = next_lambda, B_new = B_new, S_new = S_new))
}
