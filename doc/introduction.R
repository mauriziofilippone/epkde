## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(epkde)
set.seed(42)

## ----data---------------------------------------------------------------------
# Generate data from a bivariate Gaussian mixture
n <- 150; d <- 2
x <- rbind(matrix(rnorm(n/2 * d, mean = c(-2,  0)), ncol = d),
           matrix(rnorm(n/2 * d, mean = c( 2,  0)), ncol = d))

## ----iso----------------------------------------------------------------------
fit_iso <- ep_kde_isotropic(x, prior_shape = 1, prior_rate = 1)
cat("Posterior mean precision (lambda):", fit_iso$post_shape / fit_iso$post_rate, "\n")
cat("Log model evidence:", fit_iso$log_evidence, "\n")

## ----diag---------------------------------------------------------------------
fit_diag <- ep_kde_diagonal(x,
                             prior_shape = rep(1, d),
                             prior_rate  = rep(1, d))
cat("Posterior mean precision per dimension:",
    fit_diag$post_shape / fit_diag$post_rate, "\n")
cat("Log model evidence:", fit_diag$log_evidence, "\n")

## ----full---------------------------------------------------------------------
fit_full <- ep_kde_full(x, prior_nu = 0)
cat("Posterior nu:", fit_full$post_nu, "\n")
cat("Posterior mean precision matrix:\n")
print(fit_full$post_mean)
cat("Log model evidence:", fit_full$log_evidence, "\n")

## ----bf-----------------------------------------------------------------------
ev <- c(isotropic = model_evidence(fit_iso),
        diagonal  = model_evidence(fit_diag),
        full      = model_evidence(fit_full))
cat("Log evidences:\n"); print(ev)
cat("\nLog Bayes factor (diagonal vs. isotropic):",
    ev["diagonal"] - ev["isotropic"], "\n")

## ----predict------------------------------------------------------------------
# Grid for evaluation
gr    <- seq(-6, 6, length.out = 50)
grid  <- as.matrix(expand.grid(gr, gr))

# Use the diagonal fit: Lambda = diag(post_shape / post_rate)
Lambda_diag <- diag(fit_diag$post_shape / fit_diag$post_rate)
p_hat       <- kde_predict(grid, x, Lambda_diag)

# Contour plot
contour(gr, gr, matrix(p_hat, 50),
        main = "Bayesian KDE (diagonal bandwidth)",
        xlab = expression(x[1]), ylab = expression(x[2]))
points(x, pch = 20, cex = 0.4)

## ----online-------------------------------------------------------------------
n_init <- 100
x_init <- x[1:n_init, ]
x_more <- x[(n_init + 1):n, ]

fit_init   <- ep_kde_isotropic(x_init, prior_shape = 1, prior_rate = 1)
fit_update <- ep_kde_online(x_init, x_more,
                             prior_shape = 1, prior_rate = 1,
                             fit_old = fit_init)
fit_all    <- ep_kde_isotropic(x, prior_shape = 1, prior_rate = 1)

cat("Online posterior mean precision :", fit_update$post_shape / fit_update$post_rate, "\n")
cat("Offline (all-at-once) mean prec :", fit_all$post_shape    / fit_all$post_rate,    "\n")

