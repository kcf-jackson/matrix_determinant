# #Future investigation
# #Reference: http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter10.pdf
#
# arnoldi_iter <- function(A, n) {
#   arnoldi_vectors <- rnorm(ncol(A)) %>% normalise() %>% matrix()
#   for (j in 2:n) {
#     arnoldi_vectors %<>% add_new_vec(A, .)
#   }
#   arnoldi_vectors
# }
# add_new_vec <- function(A, M) {
#   k <- ncol(M)
#   qk <- A %*% M[,k]
#   for (j in 1:(k-1)) {
#     qj <- M[,j]
#     h <- sum(qj * qk)
#     qk <- qk - h * qj
#   }
#   cbind(M, normalise(qk))
# }
# 
# 
# # Doesn't give results due to numerical issues
# naive_lanczos <- function(A, m) {
#   beta <- 0
#   last_v <- 0
#   current_v <- rnorm(ncol(A)) %>% normalise()
#   for (j in 1:(m-1)) {
#     w <- A %*% current_v
#     alpha <- sum(w * current_v)
#     w <- (w - alpha * current_v - beta * last_v)
#     beta <- l2norm(w)
#     new_v <- w / beta
#     
#     last_v <- current_v
#     current_v <- new_v
#   }
#   w <- A %*% new_v
#   alpha <- sum(w * new_v)
#   return(
#     list(eigenvector = w, eigenvalue = alpha)
#   )
# }
# normalise <- function(v0) {
#   v0 / l2norm(v0)
# }