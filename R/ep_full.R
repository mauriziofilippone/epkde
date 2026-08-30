# BayesKDE: EP for a full kernel precision matrix (Lambda)
# Filippone & Sanguinetti (2011), Sections 3 & 5.
#
# Model:  KDE with Gaussian kernel K(x | Lambda) = N(x; 0, Lambda^{-1})
#         where Lambda is a d x d positive-definite precision matrix.
#         Prior: p(Lambda) = Wishart(Lambda | Lambda_0^{-1}, nu_0)
#                (nu_0 = 0 for an improper flat prior)
#         Likelihood: leave-one-out cross-validated log-likelihood (eq. 4)
#
# EP approximates the posterior with a Wishart distribution:
#   q(Lambda) = Wishart(Lambda | Sigma_new^{-1}, nu_new)
#
# The Woodbury / Sherman-Morrison identity (eq. 21 in the paper) is used to
# update the cavity precision matrix efficiently without storing n^2 d x d matrices.

#' EP inference for full KDE bandwidth (precision matrix)
#'
#' Runs Expectation Propagation to approximate the posterior distribution of
#' the full \eqn{d \times d} kernel precision matrix \eqn{\Lambda} in
#' multivariate Gaussian KDE.  The approximate posterior is a Wishart
#' distribution \eqn{q(\Lambda) = \mathcal{W}(\Lambda \mid \Sigma^{-1}, \nu)}.
#'
#' @param x Numeric matrix of training data (\eqn{n \times d}).
#' @param prior_cov \eqn{d \times d} prior covariance matrix
#'   \eqn{\Sigma_0 = \Lambda_0^{-1}}.  Pass \code{NULL} or a large-diagonal
#'   matrix for a diffuse prior.  Ignored when \code{prior_nu = 0}.
#' @param prior_nu Wishart degrees-of-freedom parameter \eqn{\nu_0} for the
#'   prior.  Use 0 for an improper flat prior.
#' @param max_iter Maximum number of EP outer iterations. Default 100.
#' @param tol Convergence tolerance on \eqn{\nu}. Default \code{1e-6}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{post_prec}}{Posterior Wishart precision parameter \eqn{\Lambda}
#'       (i.e., \eqn{E[\Lambda] / \nu}).}
#'     \item{\code{post_nu}}{Posterior Wishart degrees-of-freedom \eqn{\nu}.}
#'     \item{\code{post_mean}}{Posterior mean precision matrix
#'       \eqn{E[\Lambda] = \nu \Lambda_{new}}.}
#'     \item{\code{bandwidth_matrix}}{Posterior mean bandwidth matrix
#'       (inverse of \code{post_mean}).}
#'     \item{\code{log_evidence}}{Approximate log model evidence.}
#'     \item{\code{log_norm_const}}{Log normalization constant of the posterior.}
#'   }
#'   Returns \code{NULL} with a warning if EP fails to converge.
#'
#' @references Filippone, M. & Sanguinetti, G. (2011). Approximate inference
#'   of the bandwidth in multivariate kernel density estimation.
#'   \emph{Computational Statistics & Data Analysis}, 55(12), 3104-3122.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(50 * 2), ncol = 2)
#' fit <- ep_kde_full(x, prior_cov = NULL, prior_nu = 0)
#' cat("Posterior nu:", fit$post_nu, "\n")
#' print(fit$bandwidth_matrix)
#'
#' @export
ep_kde_full <- function(x, prior_cov = NULL, prior_nu = 0,
                         max_iter = 100, tol = 1e-6) {

  if (!is.matrix(x)) stop("x must be a matrix for the full precision case.")
  n <- nrow(x)
  d <- ncol(x)

  Id <- diag(d)

  # --- Prior covariance setup ---
  # For a diffuse (near-flat) prior we use a large covariance matrix
  if (is.null(prior_cov) || prior_nu == 0) {
    prior_cov <- Id * 1e3
    prior_nu  <- 0
  }
  log_prior_norm <- if (prior_nu == 0) 0
                    else wishart_log_norm(solve(prior_cov), prior_nu)

  # --- Initialise EP factors ---
  # Each factor contributes an approximate Wishart:
  #   f_j(Lambda) ~ Wishart(Lambda | sigma_factors[,,j]^{-1}, nu_factors[j])
  # Sigma parameterization: sum of sigma_factors^{-1} gives the posterior precision
  init_cov <- Id * 1e3   # diffuse starting posterior
  init_nu  <- d + 2

  # Initialise so that sum over factors + prior = starting posterior
  sigma_factors <- array(0, dim = c(d, d, n))
  nu_factors    <- rep(0, n)
  log_norm_factors <- rep(0, n)

  # sigma_factors[,,j] = inv( inv(init_cov)/n - inv(prior_cov)/n )
  init_sigma_j <- init_cov / n - prior_cov / n
  for (j in seq_len(n)) sigma_factors[,, j] <- init_sigma_j
  nu_factors[] <- init_nu / n - prior_nu / n + d + 1

  post_cov <- init_cov
  post_nu  <- init_nu
  log_norm_const <- 0

  # Pre-allocate Woodbury precision update matrices (one per data point)
  woodbury_mats <- array(0, dim = c(d, d, n))
  log_z         <- rep(0, n)

  # --- EP outer loop ---
  for (iter in seq_len(max_iter)) {
    post_nu_old <- post_nu

    for (j in seq_len(n)) {
      # 1. Cavity covariance matrix (sum of prior + all factors except j)
      #    cavity_cov = inv( inv(prior_cov) + sum_{i != j} inv(sigma_factors[,,i]) )
      sum_sigma <- prior_cov
      for (r in seq_len(n)) sum_sigma <- sum_sigma + sigma_factors[,, r]
      sum_sigma <- sum_sigma - sigma_factors[,, j]

      # Cavity covariance not positive definite → EP ill-defined for this observation.
      # Skip the factor update and move on.
      min_eig <- min(eigen(sum_sigma, only.values = TRUE)$values)
      if (min_eig < 0) {
        warning("EP cavity covariance not positive definite at observation ", j,
                "; skipping factor update. Consider a smaller prior_cov (tighter prior).")
        next
      }

      cavity_cov  <- sum_sigma
      cavity_nu   <- prior_nu + sum(nu_factors) - nu_factors[j] - (d + 1) * (n - 1)
      log_cavity_norm <- log_prior_norm + sum(log_norm_factors) - log_norm_factors[j]

      # 2. Woodbury update: for each r, compute
      #    U_r = Q - Q a_r a_r' Q / (1 + a_r' Q a_r),  a_r = x[j,] - x[r,]
      #    which equals the precision of the cavity marginalized over x[j,] - x[r,]
      #    (Sherman-Morrison rank-1 update, eq. 21 in paper)
      eig          <- eigen(cavity_cov)
      min_eig_val  <- min(eig$values)
      if (min_eig_val < 0) {
        eig$values   <- eig$values - (min_eig_val - 1e-3)
      }
      # Q = inv(cavity_cov) via eigendecomposition
      cavity_prec  <- eig$vectors %*% (t(eig$vectors) / eig$values)
      cavity_prec  <- (cavity_prec + t(cavity_prec)) / 2   # enforce symmetry

      diff_vecs <- t(x) - x[j, ]        # d x n matrix: diff_vecs[,r] = x[r,] - x[j,]
      prec_diffs <- cavity_prec %*% diff_vecs  # d x n: Q * a_r for each r
      denom      <- 1 + colSums(diff_vecs * prec_diffs)  # 1 + a_r' Q a_r, length n

      for (r in seq_len(n)) {
        woodbury_mats[,, r] <- cavity_prec -
                               (prec_diffs[, r] %*% t(prec_diffs[, r])) / denom[r]
      }
      woodbury_mats[,, j] <- Id   # diagonal term excluded

      tilt_alpha <- cavity_nu + 1   # degrees of freedom for Wishart weights

      # 3. Compute unnormalized log weights (log of Wishart normalization constants)
      for (r in seq_len(n))
        log_z[r] <- -wishart_log_norm(woodbury_mats[,, r], tilt_alpha)

      log_z[j] <- -Inf  # exclude self-term
      min_logz  <- min(log_z[is.finite(log_z)])
      weights   <- exp(log_z - min_logz)
      weights[j] <- 0
      norm_weights <- weights / sum(weights)

      if (any(is.nan(norm_weights))) {
        warning("EP (full) failed to converge. Returning improper result.")
        return(list(
          post_prec       = Id,
          post_nu         = d + 1,
          post_mean       = Id * (d + 1),
          bandwidth_matrix = solve(Id * (d + 1)),
          log_norm_const  = 0,
          log_evidence    = -Inf
        ))
      }

      log_norm_const_j <- log_cavity_norm - log(n - 1) - d / 2 * log(2 * pi) +
                          min_logz + log(sum(weights))

      # 4. Match Wishart moments (eqs. 13-17 in paper)
      # E1 = E[Lambda] = tilt_alpha * sum_r w_r U_r
      E1 <- matrix(0, d, d)
      for (r in seq_len(n))
        E1 <- E1 + tilt_alpha * norm_weights[r] * woodbury_mats[,, r]

      # Frobenius-norm-based variance (scalar)  (eq. 14)
      rho <- 0
      for (r in seq_len(n))
        rho <- rho + norm_weights[r] *
               ((tilt_alpha + tilt_alpha^2) * sum(woodbury_mats[,, r]^2) +
                tilt_alpha * sum(diag(woodbury_mats[,, r]))^2)
      rho <- rho - sum(E1^2)

      # New Wishart parameters (eqs. 15-17)
      new_nu    <- (sum(E1^2) + sum(diag(E1))^2) / rho
      new_sigma <- solve(E1 / new_nu)   # new posterior covariance parameter
      log_new_norm <- log_norm_const_j + wishart_log_norm(solve(new_sigma), new_nu)

      # 5. Update the j-th EP factor parameters (eqs. 18-20)
      sigma_factors[,, j] <- new_sigma - cavity_cov
      nu_factors[j]       <- new_nu - cavity_nu + (d + 1)
      log_norm_factors[j] <- log_new_norm - log_cavity_norm
    }

    if (abs(post_nu - new_nu) < tol) break
    post_cov       <- new_sigma
    post_nu        <- new_nu
    log_norm_const <- log_new_norm
  }

  post_cov       <- new_sigma
  post_nu        <- new_nu
  log_norm_const <- log_new_norm
  post_prec      <- solve(post_cov)

  list(
    post_prec        = post_prec,
    post_nu          = post_nu,
    post_mean        = post_nu * post_prec,       # E[Lambda] = nu * Lambda_new
    bandwidth_matrix = post_cov / post_nu,        # E[Sigma] = Sigma_new / nu
    log_norm_const   = log_norm_const,
    log_evidence     = log_norm_const - wishart_log_norm(post_prec, post_nu)
  )
}
