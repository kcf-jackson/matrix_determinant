# Exercise. Program a function to compute the determinant of a matrix.
compute_determinant <- function(A) {
  if (ncol(A) != nrow(A)) {
    stop("Only square matrix has determinant.")
  }
  
  if (nrow(A) == 2) {
    return(A[1,1] * A[2,2] - A[1,2] * A[2,1])
  } else {
    res <- 0
    for (i in 1:ncol(A)) {
      sign <- ifelse(i %% 2 == 1, 1, -1)
      res <- res + sign * A[1,i] * compute_determinant(A[-1, -i])
    }
  }
  res
}
