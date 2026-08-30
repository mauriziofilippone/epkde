# BayesKDE: Internal utility functions
# Filippone & Sanguinetti (2011), CSDA 55:3104-3122

# -------------------------------------------------------------------------
# Log-normalization constant of the Gamma distribution
#   Gamma(x | shape, rate)  ->  log Z = shape * log(rate) - log Gamma(shape)
# -------------------------------------------------------------------------
gamma_log_norm <- function(shape, rate) {
  shape * log(rate) - lgamma(shape)
}

# -------------------------------------------------------------------------
# Log-normalization constant of the (unnormalized) Wishart distribution
#   W(Lambda | W0, nu)  parameterized with PRECISION matrix W0
#   Uses a Cholesky decomposition for numerical stability.
# -------------------------------------------------------------------------
wishart_log_norm <- function(prec_matrix, nu) {
  d   <- nrow(prec_matrix)
  # log|prec_matrix| via Cholesky factor
  U   <- chol(prec_matrix)
  log_det <- 2 * sum(log(diag(U)))

  -nu * log_det / 2 -
    nu * d / 2 * log(2) -
    d * (d - 1) / 4 * log(pi) -
    sum(lgamma((nu + 1 - seq_len(d)) / 2))
}

# -------------------------------------------------------------------------
# Squared pairwise Euclidean distance matrix
#   Returns an n x n matrix D where D[i,j] = ||x[i,] - x[j,]||^2
#   Uses the identity ||a-b||^2 = ||a||^2 + ||b||^2 - 2 a'b for speed.
# -------------------------------------------------------------------------
#' @importFrom stats dist
#' @keywords internal
sq_dist_matrix <- function(x) {
  as.matrix(dist(x))^2
}

# -------------------------------------------------------------------------
# Squared pairwise distances between two sets of points x (n x d) and y (m x d)
#   Returns an n x m matrix D where D[i,j] = ||x[i,] - y[j,]||^2
# -------------------------------------------------------------------------
sq_dist_cross <- function(x, y) {
  row_norms_x <- rowSums(x^2)
  row_norms_y <- rowSums(y^2)
  outer(row_norms_x, row_norms_y, "+") - 2 * tcrossprod(x, y)
}

# -------------------------------------------------------------------------
# Numerically stable log-sum:  log(sum(exp(log_vals)))
# -------------------------------------------------------------------------
log_sum_exp <- function(log_vals) {
  m <- max(log_vals[is.finite(log_vals)])
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(log_vals - m)))
}

# -------------------------------------------------------------------------
# Extract n and d from a data matrix or vector
# -------------------------------------------------------------------------
get_nd <- function(x) {
  if (is.matrix(x)) {
    list(n = nrow(x), d = ncol(x))
  } else {
    list(n = length(x), d = 1L)
  }
}
