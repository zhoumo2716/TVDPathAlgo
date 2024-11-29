
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview

<!-- badges: start -->

<!-- badges: end -->

`TVDPathAlgo` is an R package designed to provide an efficient
path-following algorithm for solving Total Variation Denoising (TVD)
problems. It computes solutions for a sequence of regularization
parameters (lambda) while tracking the knot locations (change points) in
the resulting piecewise constant signal.

This main function in this package is `path_algo`. There are also other
functions involved in this packgae.

## Installation

Before installing TVDPathAlgo, ensure that you have the tvdenoising
package installed, as it is a key dependency. Both can be installed from
GitHub using the pak package:

``` r
# install.packages("pak")
pak::pkg_install("glmgen/tvdenoising")
pak::pkg_install("zhoumo2716/TVDPathAlgo")
```

## Example

This is a basic example:

``` r
library(TVDPathAlgo)

# Input signal
y <- c(1, 3, 4, 6, 8, 6)

# Solve using the path-following algorithm
result <- path_algo(y, lambda_0 = 100, method = "tri")

# View results
print(result$lambda_values)      # Sequence of lambda values
#> [1] 100.0000000   3.0000000   2.5000000   2.0000000   0.6666667   0.0000000
print(result$knot_locations)     # Knot locations for each lambda
#> [[1]]
#> [1] 0
#> 
#> [[2]]
#> [1] 3 4
#> 
#> [[3]]
#> [1] 3 4
#> 
#> [[4]]
#> [1] 2 3 4
#> 
#> [[5]]
#> [1] 2 3 4 5
#> 
#> [[6]]
#> [1] 2 3 4 5 6
```

## Functions

### `path_algo`

**Description**: Solves the Total Variation Denoising problem for a
given signal using a path-following algorithm.

**Arguments**:  
- `y`: Numeric vector (signal).  
- `lambda_0`: Non-negative numeric value for the initial lambda.  
- `method`: `"tri"` (tridiagonal optimization) or `"general"`.

**Returns**:  
- `lambda_index`: Indices corresponding to the lambda values.  
- `lambda_values`: Sequence of lambda values.  
- `knot_locations`: Knot locations for each lambda.  
- `beta_values`: Solutions at each lambda.  
- `loop_number`: Total number of iterations.

### `find_next_lambda`

**Description**: Finds the next regularization parameter (`lambda`) in
the sequence for the “general” method.

**Arguments**:  
- `y`: Numeric vector (signal).  
- `current_lambda`: Current value of the regularization parameter.  
- `B`: Current boundary set.  
- `S`: Current direction set.

**Returns**:  
- `next_lambda`: Next lambda value in the sequence.  
- `B_new`: Updated boundary set.  
- `S_new`: Updated direction set.

### `find_next_lambda_tri`

**Description**: Finds the next regularization parameter (`lambda`) in
the sequence for the “tri” (tridiagonal optimization) method. It solves
matrix inverse computation using `tridiagonal_inverse_matrix`.

**Arguments**:  
- `y`: Numeric vector (signal).  
- `current_lambda`: Current value of the regularization parameter.  
- `B`: Current boundary set.  
- `S`: Current direction set.

**Returns**:  
- `next_lambda`: Next lambda value in the sequence.  
- `B_new`: Updated boundary set.  
- `S_new`: Updated direction set.

### `tridiagonal_inverse_matrix`

**Description**: Computes the inverse of a tridiagonal matrix.

**Arguments**:  
- `A`: Tridiagonal matrix.

**Returns**:  
- Inverted matrix `A`.
