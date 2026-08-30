# BayesKDE – test-basic.R
#
# Fast tests (no skip_on_cran).  These run on every devtools::test() call and
# cover:
#   1. Internal utility correctness
#   2. Input coercion (vector vs. matrix)
#   3. 1-D model equivalence (isotropic == diagonal for d = 1)
#   4. kde_predict: log_scale consistency & approximate integration to 1
#   5. Regression snapshots (fixed seed, baked-in numerics)
#   6. Online update composes correctly (two-step ~ one-shot offline)
#   7. Prior sensitivity (informative prior shifts the posterior)
#   8. model_evidence() wrapper

library(testthat)
library(BayesKDE)

# ── shared fixtures ────────────────────────────────────────────────────────────
set.seed(99)
n <- 60; d <- 2
x60 <- matrix(rnorm(n * d), ncol = d)


# ══════════════════════════════════════════════════════════════════════════════
# 1.  Internal utility functions
# ══════════════════════════════════════════════════════════════════════════════

test_that("gamma_log_norm matches analytic formula", {
  # log Z(Gamma(shape, rate)) = shape * log(rate) - lgamma(shape)
  expect_equal(gamma_log_norm(2, 3),
               2 * log(3) - lgamma(2),
               tolerance = 1e-12)
  # Check a few more (a, b) pairs
  expect_equal(gamma_log_norm(0.5, 2),
               0.5 * log(2) - lgamma(0.5),
               tolerance = 1e-12)
  expect_equal(gamma_log_norm(10, 0.1),
               10 * log(0.1) - lgamma(10),
               tolerance = 1e-12)
})

test_that("sq_dist_matrix is symmetric with zero diagonal", {
  x3 <- matrix(c(0,0, 1,1, 2,2), ncol = 2, byrow = TRUE)
  D  <- sq_dist_matrix(x3)
  expect_true(all(diag(D) == 0))
  expect_equal(D, t(D))
  # Distance between (0,0) and (1,1): (1-0)^2 + (1-0)^2 = 2
  expect_equal(D[1, 2], 2)
})

test_that("sq_dist_cross matches sq_dist_matrix when x == y", {
  x3 <- matrix(c(0,0, 1,1, 2,2), ncol = 2, byrow = TRUE)
  D_mat   <- sq_dist_matrix(x3)
  D_cross <- sq_dist_cross(x3, x3)
  expect_equal(D_mat, D_cross, tolerance = 1e-12)
})

test_that("log_sum_exp matches naive computation", {
  v <- c(1, 2, 3)
  expect_equal(log_sum_exp(v),
               log(sum(exp(v))),
               tolerance = 1e-12)
  # Overflow-proof: huge values
  v_big <- c(1000, 1001, 1002)
  expect_true(is.finite(log_sum_exp(v_big)))
  expect_equal(log_sum_exp(v_big),
               log_sum_exp(v_big - 1000) + 1000,
               tolerance = 1e-10)
  # All -Inf returns -Inf
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)
})


# ══════════════════════════════════════════════════════════════════════════════
# 2.  Input coercion: vector input == single-column matrix
# ══════════════════════════════════════════════════════════════════════════════

test_that("ep_kde_isotropic accepts a plain vector and matches matrix result", {
  set.seed(5)
  xv <- rnorm(50)
  fv <- ep_kde_isotropic(xv, prior_shape = 1, prior_rate = 1)
  fm <- ep_kde_isotropic(matrix(xv, ncol = 1), prior_shape = 1, prior_rate = 1)
  expect_equal(fv$post_shape, fm$post_shape, tolerance = 1e-10)
  expect_equal(fv$post_rate,  fm$post_rate,  tolerance = 1e-10)
})

test_that("kde_predict accepts a plain vector for x_test", {
  set.seed(5)
  xv   <- rnorm(50)
  fit  <- ep_kde_isotropic(xv, prior_shape = 1, prior_rate = 1)
  Lam  <- fit$post_shape / fit$post_rate * diag(1)
  xtest_vec <- rnorm(5)
  p_vec <- kde_predict(xtest_vec, matrix(xv, ncol = 1), Lam)
  p_mat <- kde_predict(matrix(xtest_vec, ncol = 1), matrix(xv, ncol = 1), Lam)
  expect_equal(p_vec, p_mat, tolerance = 1e-12)
})


# ══════════════════════════════════════════════════════════════════════════════
# 3.  1-D model equivalence: isotropic == diagonal when d = 1
# ══════════════════════════════════════════════════════════════════════════════

test_that("isotropic and diagonal give identical posteriors for d = 1", {
  set.seed(1)
  x1  <- matrix(rnorm(40), ncol = 1)
  fi1 <- ep_kde_isotropic(x1, prior_shape = 1, prior_rate = 1)
  fd1 <- ep_kde_diagonal(x1,  prior_shape = 1, prior_rate = 1)

  expect_equal(fi1$post_shape / fi1$post_rate,
               fd1$post_shape / fd1$post_rate,
               tolerance = 1e-5,
               label = "posterior mean precision")

  expect_equal(fi1$log_evidence, fd1$log_evidence,
               tolerance = 1e-4,
               label = "log evidence")
})


# ══════════════════════════════════════════════════════════════════════════════
# 4.  kde_predict correctness
# ══════════════════════════════════════════════════════════════════════════════

test_that("kde_predict log_scale=TRUE equals log of log_scale=FALSE", {
  fit <- ep_kde_isotropic(x60, prior_shape = 1, prior_rate = 1)
  Lam <- fit$post_shape / fit$post_rate * diag(d)
  x_test <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)

  p  <- kde_predict(x_test, x60, Lam, log_scale = FALSE)
  lp <- kde_predict(x_test, x60, Lam, log_scale = TRUE)

  expect_equal(log(p), lp, tolerance = 1e-12)
})

test_that("kde_predict densities integrate to ~1 (1-D, fine grid)", {
  set.seed(7)
  x_tr <- matrix(rnorm(200), ncol = 1)
  fit  <- ep_kde_isotropic(x_tr, prior_shape = 1, prior_rate = 1)
  Lam  <- fit$post_shape / fit$post_rate * diag(1)

  grid <- matrix(seq(-5, 5, length.out = 2000), ncol = 1)
  dens <- kde_predict(grid, x_tr, Lam)
  dx   <- 10 / 1999  # grid spacing

  expect_equal(sum(dens) * dx, 1, tolerance = 0.01)
})

test_that("kde_predict output has correct length and is strictly positive", {
  m_test <- 15
  fit  <- ep_kde_isotropic(x60, prior_shape = 1, prior_rate = 1)
  Lam  <- fit$post_shape / fit$post_rate * diag(d)
  xts  <- matrix(rnorm(m_test * d), ncol = d)
  dens <- kde_predict(xts, x60, Lam)

  expect_length(dens, m_test)
  expect_true(all(dens > 0))
})


# ══════════════════════════════════════════════════════════════════════════════
# 5.  Regression snapshots (seed 99, n=60, d=2)
#     Expected values computed once and baked in.
# ══════════════════════════════════════════════════════════════════════════════

test_that("ep_kde_isotropic regression: fixed numeric outputs (seed 99)", {
  fit <- ep_kde_isotropic(x60, prior_shape = 1, prior_rate = 1)

  expect_equal(fit$post_shape, 17.55584, tolerance = 1e-3)
  expect_equal(fit$post_rate,  6.966547, tolerance = 1e-3)
  expect_equal(fit$log_evidence, -173.1099, tolerance = 0.05)
  expect_true(is.finite(fit$log_evidence))
})

test_that("ep_kde_diagonal regression: fixed numeric outputs (seed 99)", {
  fit <- ep_kde_diagonal(x60,
                         prior_shape = rep(1, d),
                         prior_rate  = rep(1, d))

  expect_equal(fit$post_shape[1], 8.804598, tolerance = 1e-3)
  expect_equal(fit$post_shape[2], 9.240947, tolerance = 1e-3)
  expect_equal(fit$post_rate[1],  3.426444, tolerance = 1e-3)
  expect_equal(fit$post_rate[2],  4.076326, tolerance = 1e-3)
  expect_equal(fit$log_evidence, -174.4443, tolerance = 0.05)
})

test_that("ep_kde_full regression: post_nu and finite log_evidence (seed 99)", {
  fit <- ep_kde_full(x60, prior_cov = NULL, prior_nu = 0)

  expect_equal(fit$post_nu, 59.6139, tolerance = 0.1)
  expect_true(is.finite(fit$log_evidence))
  # Posterior mean should be positive definite
  eigs <- eigen(fit$post_mean, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})


# ══════════════════════════════════════════════════════════════════════════════
# 6.  Online update: two-step composition ~ one-shot offline
# ══════════════════════════════════════════════════════════════════════════════

test_that("online EP two-step update is close to offline fit", {
  set.seed(42)
  xall <- matrix(rnorm(60 * d), ncol = d)
  x_a  <- xall[1:30,  ]
  x_b  <- xall[31:50, ]
  x_c  <- xall[51:60, ]

  fit_a    <- ep_kde_isotropic(x_a, prior_shape = 1, prior_rate = 1)
  fit_ab   <- ep_kde_online(x_a, x_b, prior_shape = 1, prior_rate = 1,
                             fit_old = fit_a)
  fit_abc  <- ep_kde_online(rbind(x_a, x_b), x_c,
                             prior_shape = 1, prior_rate = 1,
                             fit_old = fit_ab)
  fit_all  <- ep_kde_isotropic(xall, prior_shape = 1, prior_rate = 1)

  mean_2step   <- fit_abc$post_shape / fit_abc$post_rate
  mean_offline <- fit_all$post_shape / fit_all$post_rate

  # EP is approximate, so allow generous tolerance (< 20 % relative error)
  expect_lt(abs(mean_2step - mean_offline) / mean_offline, 0.20)
})


# ══════════════════════════════════════════════════════════════════════════════
# 7.  Prior sensitivity
# ══════════════════════════════════════════════════════════════════════════════

test_that("tight informative prior pulls posterior away from data-driven value", {
  set.seed(3)
  xp <- matrix(rnorm(80 * d), ncol = d)

  fit_flat  <- ep_kde_isotropic(xp, prior_shape = 0,   prior_rate = 0)
  fit_tight <- ep_kde_isotropic(xp, prior_shape = 500, prior_rate = 1)

  mean_flat  <- fit_flat$post_shape  / fit_flat$post_rate
  mean_tight <- fit_tight$post_shape / fit_tight$post_rate

  # Tight prior centred at 500 should pull the posterior mean much higher
  expect_gt(mean_tight, mean_flat * 5)
})


# ══════════════════════════════════════════════════════════════════════════════
# 8.  model_evidence() wrapper
# ══════════════════════════════════════════════════════════════════════════════

test_that("model_evidence extracts log_evidence from all three model types", {
  fi <- ep_kde_isotropic(x60, prior_shape = 1,        prior_rate = 1)
  fd <- ep_kde_diagonal( x60, prior_shape = rep(1, d), prior_rate = rep(1, d))
  ff <- ep_kde_full(     x60, prior_cov = NULL,        prior_nu = 0)

  expect_equal(model_evidence(fi), fi$log_evidence)
  expect_equal(model_evidence(fd), fd$log_evidence)
  expect_equal(model_evidence(ff), ff$log_evidence)

  expect_error(model_evidence(list(a = 1)), "log_evidence")
})
