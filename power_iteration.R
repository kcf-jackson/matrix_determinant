# An algorithm to find the greatest absolute eigenvalue of a matrix.
# Assumptions:
# 1. A has a strictly dominant eigenvalue.
# 2. The starting vector has a non-zero coefficients in the eigenspace decomposition.
power_iter <- function(A, init_b0, tol = 1e-5) {
  if (missing(init_b0)) 
    init_b0 <- rnorm(ncol(A))
  old_b <- init_b0
  
  count <- 0
  err <- 1e5
  while (err > tol) {
    new_b <- A %*% old_b
    new_b <- new_b / l2norm(new_b)
    err <- l2norm(new_b - old_b)
    old_b <- new_b
    
    count <- count + 1
    if (count %% 1000 == 0) 
      cat(count, "iterations performed.\n")
  }
  
  cat("Converged in", count, "iterations.\n")
  list(
    eigenvector = new_b,
    eigenvalue = as.numeric(
      t(Conj(new_b)) %*% A %*% new_b / sum(Conj(new_b) * new_b)
    )
  )
}
l2norm <- function(vec0) {
  sqrt(sum(vec0^2))
}