# Inverse iteration (aka inverse power iteration) for finding eigenvalues 
# and eigenvectors near an approximate eigenvalue position.
inverse_iter <- function(A, init_b0, init_lambda0, tol = 1e-5){
  d <- ncol(A)
  if (missing(init_b0)) init_b0 <- rnorm(d)
  if (missing(init_lambda0)) init_lambda0 <- rnorm(1)
  lambda <- init_lambda0
  
  old_b <- init_b0 / l2norm(init_b0)
  new_b <- solve(A - lambda * diag(d), old_b)
  new_b <- new_b / l2norm(new_b)
  err <- l2norm(new_b - old_b)
  
  count <- 0
  while(err > tol) {
    old_b <- new_b
    new_b <- solve(A - lambda * diag(d), old_b)
    new_b <- new_b / l2norm(new_b)
    err <- l2norm(new_b - old_b)
    count <- count + 1
  }  
  
  cat("Converged in", count, "iterations.\n")
  list(
    eigenvector = new_b,
    eigenvalue = as.numeric(
      t(Conj(new_b)) %*% A %*% new_b / sum(Conj(new_b) * new_b)
    )
  )
}

