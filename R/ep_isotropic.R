# BayesKDE: EP for isotropic (scalar) kernel precision
# Implements the algorithm of Filippone & Sanguinetti (2011), Section 3 & Appendix.
#
# Model:  KDE with Gaussian kernel K(x | lambda) = N(x; 0, lambda^{-1} I)
#         Prior:  p(lambda) = Gamma(lambda | prior_shape, prior_rate)
#                 (use prior_shape = prior_rate = 0 for an improper flat prior)
#         Likelihood: leave-one-out cross-validated log-likelihood (eq. 4)
#
# EP approximates the posterior p(lambda | X) with a Gamma distribution
# q(lambda) = Gamma(lambda | post_shape, post_rate)
# by iteratively refining one approximate likelihood factor per data point.

#' EP inference for isotropic KDE bandwidth
#'
#' Runs Expectation Propagation to approximate the posterior distribution of
#' the scalar kernel precision \eqn{\lambda} (inverse bandwidth squared) in
#' multivariate Gaussian KDE with an isotropic kernel.  The approximate
#' posterior is \eqn{\Gamma(\lambda \mid a, b)} and the log model evidence is
#' also returned.
#'
#' @param x Numeric matrix of training data (\eqn{n \times d}) or a numeric
#'   vector for the univariate case.
#' @param prior_shape Non-negative scalar; shape parameter \eqn{a_0} of the
#'   Gamma prior on \eqn{\lambda}.  Use 0 for an improper flat prior.
#' @param prior_rate Non-negative scalar; rate parameter \eqn{b_0} of the
#'   Gamma prior on \eqn{\lambda}.  Use 0 for an improper flat prior.
#' @param max_iter Maximum number of EP outer iterations. Default 100.
#' @param tol Convergence tolerance on the change in posterior parameters.
#'   Default \code{1e-9}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{post_shape}}{Posterior Gamma shape \eqn{a}.}
#'     \item{\code{post_rate}}{Posterior Gamma rate \eqn{b}.}
#'     \item{\code{log_evidence}}{Approximate log model evidence.}
#'     \item{\code{log_norm_const}}{Log normalization constant of the
#'       approximate posterior (for expert use).}
#'     \item{\code{shape_factors}}{Length-\eqn{n} vector of EP factor shape
#'       parameters (one per data point).}
#'     \item{\code{rate_factors}}{Length-\eqn{n} vector of EP factor rate
#'       parameters (one per data point).}
#'     \item{\code{log_norm_factors}}{Length-\eqn{n} vector of EP factor log
#'       normalization constants (one per data point).}
#'   }
#'
#' @references Filippone, M. & Sanguinetti, G. (2011). Approximate inference
#'   of the bandwidth in multivariate kernel density estimation.
#'   \emph{Computational Statistics & Data Analysis}, 55(12), 3104-3122.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(100 * 2), ncol = 2)   # 100 bivariate observations
#' fit <- ep_kde_isotropic(x, prior_shape = 1, prior_rate = 1)
#' cat("Posterior mean precision:", fit$post_shape / fit$post_rate, "\n")
#' cat("Log model evidence:", fit$log_evidence, "\n")
#'
#' @export
ep_kde_isotropic <- function(x, prior_shape = 0, prior_rate = 0,
                              max_iter = 100, tol = 1e-9) {

  dims <- get_nd(x)
  n <- dims$n
  d <- dims$d
  if (!is.matrix(x)) x <- matrix(x, ncol = 1)

  # --- Precompute half squared pairwise distances (n x n) ---
  # half_sq_dists[r, j] = ||x_r - x_j||^2 / 2
  # Diagonal set to 1 (a dummy value; diagonal terms are excluded in the sum)
  sq_dists   <- sq_dist_matrix(x)
  half_sq_dists <- sq_dists / 2
  diag(half_sq_dists) <- 1

  # --- Initialise EP factors ---
  # Each data point j contributes one approximate factor
  # f_j(lambda) ~= Gamma(lambda | .) with parameters (shape_factors[j], rate_factors[j])
  # Initialised so that the product of all factors equals the chosen starting posterior.
  init_shape <- 1.0     # starting posterior shape
  init_rate  <- 0.01    # starting posterior rate

  shape_factors    <- rep((1 - prior_shape / n + init_shape / n), n)
  rate_factors     <- rep((-prior_rate  / n + init_rate  / n),    n)
  log_norm_factors <- rep(0, n)

  # Log normalization of the prior (0 for the improper flat prior)
  log_prior_norm <- if (prior_shape == 0 && prior_rate == 0) 0
                    else gamma_log_norm(prior_shape, prior_rate)

  # Running posterior parameters (product of prior and all factors)
  post_shape  <- init_shape
  post_rate   <- init_rate
  log_norm_const <- 0



  # --- EP outer loop ---
  for (iter in seq_len(max_iter)) {
    post_shape_old <- post_shape
    post_rate_old  <- post_rate

    for (j in seq_len(n)) {
      # 1. Compute cavity distribution (posterior with factor j removed)
      #    q^{-j}(lambda) = Gamma(lambda | cavity_shape, cavity_rate)
      cavity_shape <- sum(shape_factors) + prior_shape - shape_factors[j] - (n - 1)
      cavity_rate  <- sum(rate_factors)  + prior_rate  - rate_factors[j]
      log_cavity_norm <- log_prior_norm +
                         sum(log_norm_factors) - log_norm_factors[j]

      # Cavity rate negative → EP approximation ill-defined for this observation.
      # Skip the factor update (leave it at its current value) and move on.
      if (cavity_rate < 0) {
        warning("EP cavity rate negative at observation ", j,
                "; skipping factor update. Consider a larger prior_rate.")
        next
      }
      if (cavity_shape <= 0) next

      # 2. Form the tilted distribution f_j(lambda) * q^{-j}(lambda)
      #    For fixed j, the leave-one-out likelihood for observation j is a
      #    mixture over r != j of Gamma(lambda | tilt_shape, tilt_rates[r]).
      tilt_shape <- cavity_shape + d / 2          # scalar
      tilt_rates <- cavity_rate + half_sq_dists[, j]  # length-n vector

      # 3. Compute unnormalized mixture weights (in log space for stability)
      log_weights    <- -gamma_log_norm(tilt_shape, tilt_rates)
      if (any(is.nan(log_weights))) {
        warning("EP did not converge for observation ", j,
                ". Try increasing prior_shape or prior_rate.")
        return(NULL)
      }

      # Centering trick (matches ep.hyper in original ep.r):
      # subtract midrange of non-self log-weights before exponentiating
      # to prevent underflow when weights span a large dynamic range.
      log_weights[j] <- -Inf
      a_tmp          <- mean(range(log_weights[seq_len(n) != j]))

      # Exclude the diagonal term (r == j is not part of the LOO sum)
      weights        <- exp(log_weights - a_tmp)
      weights[j]     <- 0
      norm_weights   <- weights / sum(weights)

      # Log normalization of the tilted distribution (eq. 12 in paper)
      log_norm_const_j <- log_cavity_norm - log(n - 1) -
                          d / 2 * log(2 * pi) + a_tmp + log(sum(weights))

      # 4. Match moments of a Gamma to the tilted distribution
      #    E[lambda] and Var[lambda] under the weighted mixture
      moment1 <- sum(norm_weights * tilt_shape / tilt_rates)
      moment2 <- sum(norm_weights *
                     (tilt_shape / tilt_rates^2 + (tilt_shape / tilt_rates)^2)) -
                 moment1^2

      # New Gamma parameters matching mean and variance
      new_shape <- moment1^2 / moment2
      new_rate  <- moment1   / moment2
      log_new_norm <- log_norm_const_j + gamma_log_norm(new_shape, new_rate)

      # 5. Update the j-th EP factor (by dividing out the cavity)
      log_norm_factors[j] <- log_new_norm - log_cavity_norm
      shape_factors[j]    <- new_shape + n - (sum(shape_factors) + prior_shape - shape_factors[j])
      rate_factors[j]     <- new_rate  - (sum(rate_factors)  + prior_rate  - rate_factors[j])
    }

    # Check convergence on posterior parameters
    if ((abs(post_shape - new_shape) + abs(post_rate - new_rate)) < tol) break

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
    log_evidence     = log_norm_const - gamma_log_norm(post_shape, post_rate),
    shape_factors    = shape_factors,
    rate_factors     = rate_factors,
    log_norm_factors = log_norm_factors,
    prior_shape      = prior_shape,
    prior_rate       = prior_rate
  )
}


#' Online EP update for isotropic KDE bandwidth
#'
#' Updates a previously computed EP approximation (from \code{ep_kde_isotropic})
#' to incorporate a new batch of \eqn{m} observations, without re-processing
#' the original \eqn{n} data points from scratch.
#'
#' @param x_old Numeric matrix (\eqn{n \times d}) of the original training data.
#' @param x_new Numeric matrix (\eqn{m \times d}) of the new observations to
#'   incorporate.
#' @param prior_shape,prior_rate Prior parameters (same values used in the
#'   original \code{ep_kde_isotropic} call).
#' @param fit_old The list returned by \code{ep_kde_isotropic} on \code{x_old}.
#' @param max_iter,tol Convergence control (see \code{ep_kde_isotropic}).
#'
#' @return A named list with the same structure as \code{ep_kde_isotropic},
#'   but with \code{shape_factors}, \code{rate_factors}, and
#'   \code{log_norm_factors} of length \eqn{n + m}.
#'
#' @references Filippone, M. & Sanguinetti, G. (2011). Approximate inference
#'   of the bandwidth in multivariate kernel density estimation.
#'   \emph{Computational Statistics & Data Analysis}, 55(12), 3104-3122.
#'
#' @export
ep_kde_online <- function(x_old, x_new, prior_shape, prior_rate, fit_old,
                           max_iter = 100, tol = 1e-9) {

  dims_old <- get_nd(x_old)
  dims_new <- get_nd(x_new)
  n <- dims_old$n
  m <- dims_new$n
  d <- dims_old$d

  if (!is.matrix(x_old)) x_old <- matrix(x_old, ncol = 1)
  if (!is.matrix(x_new)) x_new <- matrix(x_new, ncol = 1)

  # Retrieve EP factors from the previous fit
  shape_factors_old    <- fit_old$shape_factors
  rate_factors_old     <- fit_old$rate_factors
  log_norm_factors_old <- fit_old$log_norm_factors

  # Current posterior from old fit serves as updated prior for the new data
  updated_prior_shape <- fit_old$prior_shape + sum(shape_factors_old) - n
  updated_prior_rate  <- fit_old$prior_rate  + sum(rate_factors_old)

  log_prior_norm <- if (prior_shape == 0 && prior_rate == 0) 0
                    else gamma_log_norm(prior_shape, prior_rate)

  # --- Step 1: Refine the n old factors to account for the n+m model ---
  # Precompute distances from old points to new points and among new points
  all_points     <- rbind(x_new, x_old)         # (m+n) x d
  sq_dists_all   <- sq_dist_cross(all_points, x_new)   # (m+n) x m
  # Rows 1..m are new-to-new; rows (m+1)..(m+n) are old-to-new
  half_sq_old_to_new <- sq_dists_all[(m + 1):(m + n), ] / 2  # n x m

  shape_factors    <- shape_factors_old
  rate_factors     <- rate_factors_old
  log_norm_factors <- log_norm_factors_old

  # Posterior from old factors (starting point for refinement)
  post_shape <- updated_prior_shape
  post_rate  <- updated_prior_rate

  for (iter in seq_len(max_iter)) {
    for (j in seq_len(n)) {
      cavity_shape <- sum(shape_factors) + prior_shape - shape_factors[j] - (n - 1)
      cavity_rate  <- sum(rate_factors)  + prior_rate  - rate_factors[j]
      log_cavity_norm <- log_prior_norm +
                         sum(log_norm_factors) - log_norm_factors[j]

      # Cavity rate negative → skip factor update for this observation.
      if (cavity_rate < 0) {
        warning("Online EP cavity rate negative at old observation ", j,
                "; skipping factor update. Consider a larger prior_rate.")
        next
      }
      if (cavity_shape <= 0) next

      # Tilted distribution: mixture over m new points + collapsed old terms
      tilt_shape_new   <- cavity_shape + d / 2                    # for each new point
      tilt_rates_new   <- cavity_rate + half_sq_old_to_new[j, ]  # length m

      # Collapsed representation of old data contribution
      tilt_shape_old <- cavity_shape + shape_factors_old[j] - 1
      tilt_rate_old  <- cavity_rate  + rate_factors_old[j]

      alpha_vec  <- c(rep(tilt_shape_new, m), tilt_shape_old)
      beta_vec   <- c(tilt_rates_new,         tilt_rate_old)

      log_weights <- -gamma_log_norm(alpha_vec, beta_vec)
      # Weight the collapsed old term by its original factor normalization
      log_weights[m + 1] <- log_weights[m + 1] + log(n - 1) +
                             log_norm_factors_old[j] + d / 2 * log(2 * pi)

      if (any(is.nan(log_weights))) {
        warning("Online EP did not converge for old observation ", j, ".")
        return(NULL)
      }

      weights      <- exp(log_weights)
      norm_weights <- weights / sum(weights)

      log_norm_const_j <- log_cavity_norm - log(n + m - 1) -
                          d / 2 * log(2 * pi) + log(sum(weights))

      moment1 <- sum(norm_weights * alpha_vec / beta_vec)
      moment2 <- sum(norm_weights *
                     (alpha_vec / beta_vec^2 + (alpha_vec / beta_vec)^2)) -
                 moment1^2

      new_shape    <- moment1^2 / moment2
      new_rate     <- moment1   / moment2
      log_new_norm <- log_norm_const_j + gamma_log_norm(new_shape, new_rate)

      log_norm_factors[j] <- log_new_norm - log_cavity_norm
      shape_factors[j]    <- new_shape + n - (sum(shape_factors) + prior_shape - shape_factors[j])
      rate_factors[j]     <- new_rate  - (sum(rate_factors)  + prior_rate  - rate_factors[j])
    }

    if ((abs(post_shape - new_shape) + abs(post_rate - new_rate)) < tol) break
    post_shape <- new_shape
    post_rate  <- new_rate
  }

  # Save refined old-factor posterior as updated prior for new factors
  new_prior_shape <- new_shape
  new_prior_rate  <- new_rate
  log_new_prior_norm <- log_new_norm

  saved_shape_factors    <- shape_factors
  saved_rate_factors     <- rate_factors
  saved_log_norm_factors <- log_norm_factors

  # --- Step 2: Run EP for the m new factors ---
  # Use full (n+m) x m distance matrix: rows 1..m = new-to-new, rows (m+1)..(m+n) = old-to-new
  half_sq_dists_full <- sq_dists_all / 2   # (n+m) x m
  # Zero out self-distances for new points (diagonal of top-left m x m block)
  for (k in seq_len(m)) half_sq_dists_full[k, k] <- 1

  shape_factors_new    <- rep(1,    m)
  rate_factors_new     <- rep(0.01, m)
  log_norm_factors_new <- rep(0,    m)

  for (iter in seq_len(max_iter)) {
    for (j in seq_len(m)) {
      cavity_shape <- sum(shape_factors_new) + new_prior_shape -
                      shape_factors_new[j] - (m - 1)
      cavity_rate  <- sum(rate_factors_new) + new_prior_rate -
                      rate_factors_new[j]
      log_cavity_norm <- log_new_prior_norm +
                         sum(log_norm_factors_new) - log_norm_factors_new[j]

      # Cavity rate negative → skip factor update for this observation.
      if (cavity_rate < 0) {
        warning("Online EP cavity rate negative at new observation ", j,
                "; skipping factor update. Consider a larger prior_rate.")
        next
      }
      if (cavity_shape <= 0) next

      tilt_shape <- cavity_shape + d / 2
      tilt_rates <- cavity_rate + half_sq_dists_full[, j]  # length n+m

      log_weights    <- -gamma_log_norm(tilt_shape, tilt_rates)
      weights        <- exp(log_weights)
      weights[j]     <- 0   # zero out self (new point j is row j)
      norm_weights   <- weights / sum(weights)

      log_norm_const_j <- log_cavity_norm - log(n + m - 1) -
                          d / 2 * log(2 * pi) + log(sum(weights))

      moment1 <- sum(norm_weights * tilt_shape / tilt_rates)
      moment2 <- sum(norm_weights *
                     (tilt_shape / tilt_rates^2 + (tilt_shape / tilt_rates)^2)) -
                 moment1^2

      new_shape    <- moment1^2 / moment2
      new_rate     <- moment1   / moment2
      log_new_norm <- log_norm_const_j + gamma_log_norm(new_shape, new_rate)

      log_norm_factors_new[j] <- log_new_norm - log_cavity_norm
      shape_factors_new[j]    <- new_shape + m -
                                 (sum(shape_factors_new) + new_prior_shape - shape_factors_new[j])
      rate_factors_new[j]     <- new_rate  -
                                 (sum(rate_factors_new) + new_prior_rate - rate_factors_new[j])
    }

    if ((abs(post_shape - new_shape) + abs(post_rate - new_rate)) < tol) break
    post_shape <- new_shape
    post_rate  <- new_rate
  }

  post_shape     <- new_shape
  post_rate      <- new_rate
  log_norm_const <- log_new_norm

  list(
    post_shape       = post_shape,
    post_rate        = post_rate,
    log_norm_const   = log_norm_const,
    log_evidence     = log_norm_const - gamma_log_norm(post_shape, post_rate),
    shape_factors    = c(saved_shape_factors,    shape_factors_new),
    rate_factors     = c(saved_rate_factors,     rate_factors_new),
    log_norm_factors = c(saved_log_norm_factors, log_norm_factors_new),
    prior_shape      = prior_shape,
    prior_rate       = prior_rate
  )
}
