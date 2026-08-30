# BayesKDE: KDE density evaluation and model comparison utilities

#' Evaluate kernel density estimate at test points
#'
#' Given training data and a fitted precision parameter (from one of the
#' \code{ep_kde_*} functions), evaluates the KDE
#' \deqn{\hat{p}(\mathbf{x}) = \frac{1}{n} \sum_{j=1}^n
#'   \mathcal{N}(\mathbf{x}; \mathbf{x}_j, \Lambda^{-1})}
#' at a set of test points.
#'
#' @param x_test  Numeric matrix (\eqn{m \times d}) of test locations, or a
#'   numeric vector for the univariate case.
#' @param x_train Numeric matrix (\eqn{n \times d}) of training data used to
#'   fit the KDE.
#' @param precision_matrix \eqn{d \times d} positive-definite precision matrix
#'   \eqn{\Lambda}.  For the isotropic fit use \code{fit$post_shape /
#'   fit$post_rate * diag(d)}.  For the diagonal fit use
#'   \code{diag(fit$post_shape / fit$post_rate)}.  For the full fit use
#'   \code{fit$post_mean}.
#' @param log_scale Logical; if \code{TRUE} return log-densities.
#'   Default \code{FALSE}.
#'
#' @return Numeric vector of length \eqn{m} with the (log-)density values at
#'   \code{x_test}.
#'
#' @examples
#' set.seed(1)
#' x_train <- matrix(rnorm(100 * 2), ncol = 2)
#' fit <- ep_kde_isotropic(x_train, prior_shape = 1, prior_rate = 1)
#' Lambda <- fit$post_shape / fit$post_rate * diag(2)
#' x_test <- matrix(c(0, 0, 1, 1), ncol = 2)
#' p_hat  <- kde_predict(x_test, x_train, Lambda)
#'
#' @export
kde_predict <- function(x_test, x_train, precision_matrix, log_scale = FALSE) {

  if (!is.matrix(x_test))  x_test  <- matrix(x_test,  ncol = 1)
  if (!is.matrix(x_train)) x_train <- matrix(x_train, ncol = 1)
  if (!is.matrix(precision_matrix)) precision_matrix <- matrix(precision_matrix)

  m  <- nrow(x_test)
  n  <- nrow(x_train)
  d  <- ncol(x_train)

  stopifnot(ncol(x_test) == d,
            nrow(precision_matrix) == d, ncol(precision_matrix) == d)

  # Log normalization constant of the Gaussian kernel: d/2 * log(det(Lambda)/(2pi)^d)
  chol_prec     <- chol(precision_matrix)
  log_det_prec  <- 2 * sum(log(diag(chol_prec)))
  log_kernel_norm <- 0.5 * (log_det_prec - d * log(2 * pi))

  # Squared Mahalanobis distances: (x_test[i,] - x_train[j,])' Lambda (x_test[i,] - x_train[j,])
  # Result: m x n matrix
  # Lower Cholesky factor L = t(U): ||delta %*% L||^2 = delta Lambda t(delta) (correct Mahalanobis)
  L       <- t(chol_prec)
  sq_maha <- sq_dist_cross(x_test %*% L, x_train %*% L)

  # log p(x_test[i]) = log(1/n) + log_sum_exp over j of (log_kernel_norm - sq_maha[i,j]/2)
  log_contrib <- log_kernel_norm - sq_maha / 2   # m x n
  log_density  <- apply(log_contrib, 1, log_sum_exp) - log(n)

  if (log_scale) log_density else exp(log_density)
}


#' Extract the log model evidence for model comparison
#'
#' A convenience wrapper that returns the log model evidence from a fitted
#' \code{ep_kde_*} object.  The evidence can be used to compute Bayes factors
#' comparing, e.g., isotropic vs. diagonal vs. full precision structures.
#'
#' @param fit A list returned by \code{ep_kde_isotropic}, \code{ep_kde_diagonal},
#'   or \code{ep_kde_full}.
#'
#' @return A single numeric value: the approximate log model evidence.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 2), ncol = 2)
#' fit_iso  <- ep_kde_isotropic(x, prior_shape = 1, prior_rate = 1)
#' fit_diag <- ep_kde_diagonal(x,  prior_shape = rep(1, 2), prior_rate = rep(1, 2))
#' fit_full <- ep_kde_full(x, prior_cov = NULL, prior_nu = 0)
#'
#' cat("Log evidence  isotropic:", model_evidence(fit_iso),  "\n")
#' cat("Log evidence  diagonal: ", model_evidence(fit_diag), "\n")
#' cat("Log evidence  full:     ", model_evidence(fit_full), "\n")
#'
#' # Bayes factor: full vs. isotropic (on log scale)
#' cat("Log BF (full vs. iso):", model_evidence(fit_full) - model_evidence(fit_iso), "\n")
#'
#' @export
model_evidence <- function(fit) {
  if (is.null(fit$log_evidence))
    stop("The supplied object does not contain a 'log_evidence' element.")
  fit$log_evidence
}
