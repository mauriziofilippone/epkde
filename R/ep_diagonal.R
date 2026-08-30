# BayesKDE: EP for diagonal kernel precision matrix
# Filippone & Sanguinetti (2011), Appendix.
#
# Model:  KDE with Gaussian kernel K(x | lambda) = N(x; 0, diag(lambda)^{-1})
#         where lambda = (lambda_1, ..., lambda_d) is a vector of per-dimension
#         precisions.  Each lambda_k has an independent Gamma prior.
#         Likelihood: leave-one-out cross-validated log-likelihood (eq. 4)
#
# EP approximates the posterior with independent Gamma marginals:
#   q(lambda_k) = Gamma(lambda_k | post_shape[k], post_rate[k])

#' EP inference for diagonal KDE bandwidth
#'
#' Runs Expectation Propagation to approximate the posterior distribution of
#' the diagonal kernel precision \eqn{\boldsymbol{\lambda} = (\lambda_1,
#' \ldots, \lambda_d)} in multivariate Gaussian KDE.  Each dimension is
#' assigned an independent Gamma prior.
#'
#' @param x Numeric matrix of training data (\eqn{n \times d}).
#' @param prior_shape Numeric vector of length \eqn{d}; shape parameters
#'   \eqn{a_{0k}} for the Gamma prior on each \eqn{\lambda_k}.
#'   Use \code{rep(0, d)} for an improper flat prior.
#' @param prior_rate  Numeric vector of length \eqn{d}; rate parameters
#'   \eqn{b_{0k}} for the Gamma prior on each \eqn{\lambda_k}.
#'   Use \code{rep(0, d)} for an improper flat prior.
#' @param max_iter Maximum number of EP outer iterations. Default 100.
#' @param tol Convergence tolerance. Default \code{1e-6}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{post_shape}}{Length-\eqn{d} vector of posterior Gamma shapes.}
#'     \item{\code{post_rate}}{Length-\eqn{d} vector of posterior Gamma rates.}
#'     \item{\code{log_evidence}}{Approximate log model evidence.}
#'     \item{\code{log_norm_const}}{Log normalization constant of the posterior.}
#'     \item{\code{shape_factors}}{\eqn{n \times d} matrix of EP factor shapes.}
#'     \item{\code{rate_factors}}{\eqn{n \times d} matrix of EP factor rates.}
#'     \item{\code{log_norm_factors}}{Length-\eqn{n} vector of EP factor log
#'       normalization constants.}
#'   }
#'
#' @references Filippone, M. & Sanguinetti, G. (2011). Approximate inference
#'   of the bandwidth in multivariate kernel density estimation.
#'   \emph{Computational Statistics & Data Analysis}, 55(12), 3104-3122.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 2), ncol = 2)
#' fit <- ep_kde_diagonal(x, prior_shape = rep(1, 2), prior_rate = rep(1, 2))
#' cat("Posterior mean precision:", fit$post_shape / fit$post_rate, "\n")
#'
#' @export
ep_kde_diagonal <- function(x, prior_shape, prior_rate,
                             max_iter = 100, tol = 1e-6) {

  if (!is.matrix(x)) stop("x must be a matrix for the diagonal case.")
  n <- nrow(x)
  d <- ncol(x)
  stopifnot(length(prior_shape) == d, length(prior_rate) == d)

  # --- Precompute per-dimension half squared distances ---
  # half_sq_dists[r, j, k] = (x[r,k] - x[j,k])^2 / 2
  # Store as array of d matrices, each n x n
  half_sq_dists <- array(0, dim = c(n, n, d))
  for (k in seq_len(d)) {
    mat                 <- as.matrix(dist(x[, k]))^2 / 2
    diag(mat)           <- 1   # diagonal excluded in EP sum; value is dummy
    half_sq_dists[,,k] <- mat
  }

  # --- Initialise EP factors ---
  # shape_factors[j, k], rate_factors[j, k]: parameters of j-th factor for dim k
  init_shape <- rep(1,    d)
  init_rate  <- rep(0.01, d)

  shape_factors    <- matrix(0, n, d)
  rate_factors     <- matrix(0, n, d)
  for (k in seq_len(d)) {
    shape_factors[, k] <- 1 - prior_shape[k] / n + init_shape[k] / n
    rate_factors[, k]  <- -prior_rate[k]  / n + init_rate[k]  / n
  }
  log_norm_factors <- rep(0, n)

  # Log normalization of the prior
  log_prior_norm <- 0
  nonzero_k <- which(prior_shape != 0)
  for (k in nonzero_k)
    log_prior_norm <- log_prior_norm + gamma_log_norm(prior_shape[k], prior_rate[k])

  post_shape  <- init_shape
  post_rate   <- init_rate
  log_norm_const <- 0

  # --- EP outer loop ---
  for (iter in seq_len(max_iter)) {
    post_shape_old <- post_shape
    post_rate_old  <- post_rate

    for (j in seq_len(n)) {
      # 1. Cavity distribution parameters (vector over dimensions)
      cavity_shape <- colSums(shape_factors) + prior_shape - shape_factors[j, ] - (n - 1)
      cavity_rate  <- colSums(rate_factors)  + prior_rate  - rate_factors[j, ]
      log_cavity_norm <- log_prior_norm + sum(log_norm_factors) - log_norm_factors[j]

      # Cavity rate negative in any dimension → EP ill-defined for this observation.
      # Skip the factor update and move on.
      neg_dims <- which(cavity_rate < 0)
      if (length(neg_dims) > 0) {
        warning("EP cavity rate negative in dimension(s) ",
                paste(neg_dims, collapse = ", "), " at observation ", j,
                "; skipping factor update. Consider larger prior_rate values.")
        next
      }

      # 2. Tilted distribution parameters
      #    For the diagonal case the tilted distribution is a product of d
      #    Gamma mixtures, one per dimension.
      tilt_shape <- cavity_shape + 0.5                   # length-d vector
      # tilt_rates is n x d:  tilt_rates[r, k] = cavity_rate[k] + half_sq_dists[r,j,k]
      tilt_rates <- matrix(0, n, d)
      for (k in seq_len(d))
        tilt_rates[, k] <- cavity_rate[k] + half_sq_dists[, j, k]

      # 3. Unnormalized log weights: sum over dimensions of -log Z_Gamma
      log_weights <- rep(0, n)
      for (k in seq_len(d))
        log_weights <- log_weights - gamma_log_norm(tilt_shape[k], tilt_rates[, k])

      if (any(is.nan(log_weights))) {
        warning("EP did not converge for observation ", j,
                ". Try increasing prior_shape or prior_rate.")
        return(NULL)
      }

      weights      <- exp(log_weights)
      weights[j]   <- 0
      norm_weights <- weights / sum(weights)

      log_norm_const_j <- log_cavity_norm - log(n - 1) -
                          d / 2 * log(2 * pi) + log(sum(weights))

      # 4. Match Gamma moments per dimension
      # moment1[k] = E[lambda_k], moment2[k] = Var[lambda_k]
      moment1 <- colSums(norm_weights * t(tilt_shape / t(tilt_rates)))
      moment2 <- colSums(norm_weights *
                   (t(tilt_shape / t(tilt_rates^2)) +
                    t((tilt_shape / t(tilt_rates))^2))) - moment1^2

      new_shape    <- moment1^2 / moment2
      new_rate     <- moment1   / moment2
      log_new_norm <- log_norm_const_j + sum(gamma_log_norm(new_shape, new_rate))

      # 5. Update j-th EP factor
      log_norm_factors[j] <- log_new_norm - log_cavity_norm
      shape_factors[j, ]  <- new_shape - cavity_shape + 1
      rate_factors[j, ]   <- new_rate  - cavity_rate
    }

    if (sum(abs(post_shape - new_shape) + abs(post_rate - new_rate)) < tol) break
    post_shape     <- new_shape
    post_rate      <- new_rate
    log_norm_const <- log_new_norm
  }

  post_shape     <- new_shape
  post_rate      <- new_rate
  log_norm_const <- log_new_norm

  list(
    post_shape       = post_shape,
    post_rate        = post_rate,
    log_norm_const   = log_norm_const,
    log_evidence     = log_norm_const - sum(gamma_log_norm(post_shape, post_rate)),
    shape_factors    = shape_factors,
    rate_factors     = rate_factors,
    log_norm_factors = log_norm_factors
  )
}
