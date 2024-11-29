#' Compute the Inverse of a Tridiagonal Matrix
#'
#' This function calculates the inverse of a tridiagonal matrix using forward 
#' elimination and back substitution. It solves the equation \eqn{Ax = e_i} 
#' for each column vector \eqn{e_i} (standard basis vectors) to construct the 
#' inverse matrix.
#'
#' @param A A numeric square matrix of dimension \eqn{n \times n}, where \eqn{A} 
#' must be tridiagonal. The matrix should have non-zero entries only on the main 
#' diagonal, superdiagonal, and subdiagonal.
#'
#' @return A numeric matrix of the same dimension as \eqn{A}, representing the inverse 
#' of the input tridiagonal matrix.
#'
#' @examples
#' # Example tridiagonal matrix
#' A <- matrix(c(4, 1, 0, 
#'               1, 4, 1,
#'               0, 1, 4), byrow = TRUE, nrow = 3)
#'
#' # Compute the inverse
#' A_inv <- tridiagonal_inverse_matrix(A)
#' print(A_inv)
#'
#' # Verify the result
#' all.equal(A %*% A_inv, diag(3)) # Should return TRUE
#'
#' @export
tridiagonal_inverse_matrix <- function(A) {
  n <- nrow(A)  # Dimension of the matrix
  A_inv <- matrix(0, n, n)  # Placeholder for the inverse matrix
  
  # Extract diagonals
  a <- c(0, A[row(A) == col(A) + 1])   # Subdiagonal (length n-1, padded with 0 at the start)
  b <- A[row(A) == col(A)]             # Main diagonal (length n)
  c <- c(A[row(A) == col(A) - 1], 0)   # Superdiagonal (length n-1, padded with 0 at the end)
  
  # Solve Ax = e_i for each column vector e_i
  for (i in 1:n) {
    # Create right-hand side (standard basis vector e_i)
    e_i <- numeric(n)
    e_i[i] <- 1
    
    # Forward elimination
    cp <- c  # Copy of superdiagonal
    bp <- b  # Copy of main diagonal
    dp <- e_i  # Copy of right-hand side
    for (j in 2:n) {
      w <- a[j] / bp[j - 1]  # Eliminate subdiagonal entries
      bp[j] <- bp[j] - w * cp[j - 1]  # Update main diagonal
      dp[j] <- dp[j] - w * dp[j - 1]  # Update right-hand side
    }
    
    # Back substitution
    x <- numeric(n)
    x[n] <- dp[n] / bp[n]
    for (j in (n - 1):1) {
      x[j] <- (dp[j] - cp[j] * x[j + 1]) / bp[j]
    }
    
    # Store the solution as a column in the inverse matrix
    A_inv[, i] <- x
  }
  
  return(A_inv)
}

