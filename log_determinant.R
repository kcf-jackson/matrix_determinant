# Randomized linear-time algorithm
init_fun <- function(delta) {
  return(
    function(x) {
      log(1 - 0.5 * ((1 - 2*delta) * x + 1))
    }
  )
}

# Compute log determinant for positive definite matrices with sigma_max < 1.
# delta takes values from 0 to 0.5 exclusively. 
log_det_pd <- function(B, m, n, delta) {
  d <- nrow(B)
  A <- diag(d) - B
  G <- 0
  
  cheb_coeff <- pracma::chebCoeff(init_fun(delta), -1, 1, n)
  for (i in 1:m) {
    v <- sample(c(0,1), d, replace = T)
    u <- cheb_coeff[1] * v
    if (n > 1) {
      w0 <- v
      w1 <- A %*% v
      u <- u + cheb_coeff[2] * A %*% v
      for (j in 2:n) {
        w2 <- 2 * A %*% w1 - w0
        u <- u + cheb_coeff[j+1] * w2
        w0 <- w1
        w1 <- w2
      }
    }
    G <- G + sum(v * u) / m
  }
  G
}
log_det_pd.guarantee <- function(epsilon, delta) {
  cat("m needs to be greater than", 54 * epsilon^{-2} * log(2 / delta), "\n")
  UP <- log( (20 / epsilon) * (sqrt(2 / delta - 1) - 1) * (log(2/delta) / log(1 / (1 - delta))) )
  DOWN <- log( (sqrt(2 - delta) + sqrt(delta)) / (sqrt(2 + delta) - sqrt(delta)) )
  cat("n needs to be greater than", UP / DOWN, "\n")
}


log_det <- function(C, m, n, sigma_min, sigma_max) {
  MORE WORK HERE FOR FINDING THE MINIMUM AND MAXIMUM EIGENVALUES ESTIMATE
  # if (missing(sigma_min)) {
  #   sigma_min <- rARPACK::eigs(C , 1, which = "SM")$values
  #   if (sigma_min <= 0) {
  #     stop("sigma_min has to be greater than 0.")
  #   } 
  # }
  # if (missing(sigma_max)) {
  #   sigma_max <- rARPACK::eigs(C , 1, which = "LM")$values
  #   # sigma_max <- sqrt(L1_norm(C) * Linf_norm(C))
  # }
  d <- nrow(C)
  B <- t(C) %*% C / (sigma_min^2 + sigma_max^2)
  delta <- sigma_min^2 / (sigma_min^2 + sigma_max^2)
  G <- log_det_pd(B, m, n, delta)
  0.5 * (G + d * log(sigma_min^2 + sigma_max^2))
}
L1_norm <- function(C) {
  sum(abs(C))
}
Linf_norm <- function(C) {
  max(abs(C))
}
