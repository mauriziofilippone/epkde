# epkde

**Bayesian Bandwidth Selection for Multivariate KDE via Expectation Propagation**

`epkde` implements the approximate Bayesian method for bandwidth selection in
multivariate kernel density estimation (KDE) from Filippone & Sanguinetti
(2011). The method uses the Expectation Propagation (EP) algorithm to
approximate the posterior distribution of the kernel precision matrix under a
leave-one-out cross-validated likelihood.

Three covariance structures are supported:

- **Isotropic** — scalar precision (fastest)
- **Diagonal** — per-dimension precisions
- **Full** — unconstrained precision matrix

An online variant is also provided for the isotropic case, allowing incremental
updates as new data arrive.

## Installation

```r
# From GitHub
devtools::install_github("mauriziofilippone/epkde")
```

## Quick start

```r
library(epkde)
set.seed(1)

# Fit isotropic EP bandwidth to bivariate data
x <- matrix(rnorm(200 * 2), ncol = 2)
fit <- ep_kde_isotropic(x, prior_shape = 1, prior_rate = 1)

# Posterior mean precision
lambda_hat <- fit$post_shape / fit$post_rate
cat("Posterior mean precision:", lambda_hat, "\n")
cat("Log model evidence:", fit$log_evidence, "\n")

# Evaluate the KDE at test points
x_test <- matrix(c(0, 0, 1, 1), ncol = 2)
Lambda  <- lambda_hat * diag(2)
p_hat   <- kde_predict(x_test, x, Lambda)

# Compare models via Bayes factors
fit_diag <- ep_kde_diagonal(x, prior_shape = rep(1, 2), prior_rate = rep(1, 2))
cat("Log BF (diagonal vs isotropic):",
    model_evidence(fit_diag) - model_evidence(fit), "\n")
```

## Functions

| Function | Description |
|---|---|
| `ep_kde_isotropic()` | EP for scalar (isotropic) precision |
| `ep_kde_diagonal()` | EP for diagonal precision matrix |
| `ep_kde_full()` | EP for full precision matrix |
| `ep_kde_online()` | Online EP update for isotropic precision |
| `kde_predict()` | Evaluate the KDE at test points |
| `model_evidence()` | Extract log model evidence for model comparison |

## Reference

Filippone, M. & Sanguinetti, G. (2011). Approximate inference of the bandwidth
in multivariate kernel density estimation. *Computational Statistics & Data
Analysis*, 55(12), 3104–3122. <https://doi.org/10.1016/j.csda.2011.05.023>

## License

GPL-3 © Maurizio Filippone
