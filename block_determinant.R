# Exercise. Program a function to compute the determinant of a matrix
# using block-based approach.
library(magrittr)
compute_determinant_by_block <- function(A, min_block_size = 10, 
                                         engine = Rcpp_compute_determinant,
                                         ...) {
  if (length(A) == 1) return(A)
  if (ncol(A) != nrow(A)) 
    stop("Only square matrix has determinant.")
  if (nrow(A) <= min_block_size) return(engine(A, ...))
  
  midpt <- floor(ncol(A) / 2)
  num_col <- ncol(A)
  block_matrices <- list(
    A1 = A[1:midpt, 1:midpt],
    B1 = A[1:midpt, (midpt + 1):num_col],
    C1 = A[(midpt + 1):num_col, 1:midpt],
    D1 = A[(midpt + 1):num_col, (midpt + 1):num_col]
  )
  if (matrixcalc::is.singular.matrix(block_matrices$A1)) 
    stop("Block (1,1) is not invertible.")
  
  attach(block_matrices)
  M1 <- A1
  M2 <- D1 - C1 %*% solve(A1) %*% B1
  detach(block_matrices)
  
  block_matrices <- list(M1 = M1, M2 = M2)
  block_matrices %>%
    lapply(compute_determinant_by_block) %>%
    do.call(c, .) %>% prod()
}
