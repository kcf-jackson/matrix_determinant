# Test.
A <- matrix(rnorm(100), 10, 10);


# Naive implementation (cofactor expansion)
source("determinant.R")
system.time(d1 <- compute_determinant(A))
system.time(d2 <- det(A))
all.equal(d1, d2)


# Rcpp translation of the naive implementation
Rcpp::sourceCpp("determinant.cpp")
system.time(d3 <- Rcpp_compute_determinant(A))
all.equal(d3, d2)
# significant speedup but doesn't scale
# A 11 x 11 matrix would take about 10 times as long due to the nature
# of cofactor expansion.
A <- matrix(rnorm(121), 11, 11)
system.time(d2 <- det(A))
system.time(d3 <- Rcpp_compute_determinant(A))
all.equal(d3, d2)


# Try naive block-based implementation (without parallelisation)
# Idea; if we split a matrix into 4 blocks [M] = [A B; C D],
# if A is invertible, then det(M) = det(A) * det(D - C A_inverse B})
A <- matrix(rnorm(1000000), 1000, 1000)
system.time(d2 <- det(A))
system.time(d3 <- compute_determinant_by_block(A))
all.equal(d3, d2)
# Can handle much larger matrices (d = 1000 in about a second)


# Limit of base::det function
# The base package uses LAPACK subroutines (which is based on 
# LU decomposition)
c(10, 20, 50, 100, 1000, 5000) %>% 
  purrr::map(
    function(x) {
      A <- matrix(rnorm(i^2), i, i)
      system.time(det(A))
      # system.time(determinant(A, T)) #log-determinant
    }
  )
# The performance starts to suffer when d = 5000


# Randomized linear-complexity algorithm
A <- matrix(rnorm(1000000), 1000, 1000)
pdA <- t(A) %*% A
system.time(d2 <- determinant(pdA)$modulus)
system.time(d4 <- log_det(pdA, m = 30, n = 15))
print(d2); print(d4)
