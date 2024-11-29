test_that("tridiagonal_inverse_matrix computes the same result as solve()", {
  # Helper function to create a random tridiagonal matrix
  create_tridiagonal_matrix <- function(n, main_diag, sub_diag, super_diag) {
    A <- matrix(0, n, n)
    diag(A) <- main_diag
    diag(A[-1, ]) <- sub_diag
    diag(A[, -1]) <- super_diag
    return(A)
  }
  
  # Test cases
  test_cases <- list(
    # Case 1: Same diagonal values
    list(n = 5, main_diag = rep(4, 5), sub_diag = rep(-1, 4), super_diag = rep(-1, 4)),
    
    # Case 2: Different diagonal values
    list(n = 6, main_diag = c(3, 4, 5, 6, 7, 8), sub_diag = c(-1, -2, -1, -3, -2), super_diag = c(1, 2, 3, 2, 1)),
    
    # Case 3: Larger matrix with different values
    list(n = 10, main_diag = 2:11, sub_diag = rep(-1, 9), super_diag = rep(1, 9))
  )
  
  for (case in test_cases) {
    # Generate the matrix
    A <- create_tridiagonal_matrix(case$n, case$main_diag, case$sub_diag, case$super_diag)
    
    # Compute the inverse using tridiagonal_inverse_matrix
    A_inv_custom <- tridiagonal_inverse_matrix(A)
    
    # Compute the inverse using solve()
    A_inv_solve <- solve(A)
    
    # Check if the custom result matches solve()
    expect_true(
      all.equal(A_inv_custom, A_inv_solve, tolerance = 1e-8),
      info = paste("Failed for matrix:\n", paste(capture.output(print(A)), collapse = "\n"))
    )
  }
})
