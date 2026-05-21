# Correcting Measurement Error and Zero Inflation in Functional Covariates for Scalar-on-Function Quantile Regression

Welcome! This repository is the code companion for our paper **"Correcting
Measurement Error and Zero Inflation in Functional Covariates for
Scalar-on-Function Quantile Regression"**, by Caihong Qin, Lan Xue, Ufuk
Beyaztas, Roger S. Zoh, Mark Benden, Jeff Goldsmith, and Carmen D. Tekwe.

The paper was published in *Statistics in Medicine*, with DOI:
https://doi.org/10.1002/sim.70598.

The paper studies scalar-on-function quantile regression when the functional
covariate is observed through repeated measurements that contain both
measurement error and excess zeros. The main idea is to first recover the latent
functional covariate by correcting zero inflation and measurement error, and
then plug the recovered curve into joint quantile regression across multiple
quantile levels.

## Example

The example below generates sample data, applies the proposed BE-ZIME
correction, and then shows how to use the recovered curves in joint quantile
regression.

Install the required packages once before running the example:

```r
install.packages(c("MASS", "lme4", "quantregGrowth"))
```

```r
source("zime_functions.R")

# Generate latent curves, zero-inflated repeated measurements, and a scalar response.
sample_data <- generate_zime_sample(
  n = 100,                  # sample size
  grid_length = 100,        # number of observed time points
  n_replicates = 7,         # number of repeated measurements per subject
  n_segments = 2,           # number of zero-inflation intervals
  first_segment_prob = 0.4, # nonzero probability in the first segment
  pattern = "alternating",  # piecewise pattern for nonzero probabilities
  seed = 2026               # random seed for reproducibility
)

# Apply BE-ZIME to recover latent curves and zero-inflation probabilities.
zime_fit <- estimate_be_zime(
  w = sample_data$w,              # observed zero-inflated replicated measurements
  grid = sample_data$grid,        # observation grid
  basis_df = 6,                   # number of full B-spline basis functions
  basis_degree = 3,               # cubic B-spline basis
  n_segments = 2                  # working number of zero-inflation intervals
)

zime_fit$converged
zime_fit$iterations

# Summarize correction accuracy in the generated example.
zero_inflation_mse <- mean((zime_fit$zero_inflation_hat - sample_data$zero_inflation)^2)
x_mse <- mean((zime_fit$x_hat - sample_data$x)^2)

zero_inflation_mse
x_mse

# Use the recovered curves in joint quantile regression when quantregGrowth is installed.
if (requireNamespace("quantregGrowth", quietly = TRUE)) {
  qfit <- fit_joint_quantile_regression(
    y = sample_data$y,             # scalar outcome
    x_hat = zime_fit$x_hat,        # recovered latent functional covariate
    z = sample_data$z,             # additional error-free scalar covariates
    basis_object = zime_fit$basis, # basis used to compute functional scores
    taus = c(0.25, 0.5, 0.75)      # quantile levels fit jointly
  )

  qfit$check_loss
  qfit$beta_grid
}
```

To use your own data, prepare `w` as an array with dimensions
`n_subjects x n_grid_points x n_replicates`. The entries should be nonnegative
measurements of the same underlying functional covariate. Choose `n_segments`
as the working number of piecewise-constant zero-inflation intervals. The true
and estimated nonzero probabilities are returned as
`n_subjects x n_segments` matrices.

## Details

The repository includes:

- `zime_functions.R`: one sample data-generation function, one correction
  function, and a function for joint quantile regression using
  `quantregGrowth::gcrq`.

Required R packages:

- `MASS`, `splines`, and `lme4` for the BE-ZIME correction.
- `quantregGrowth` for the joint quantile regression step.

Main functions:

- `generate_zime_sample()`: generates latent curves, zero-inflated replicated
  measurements, scalar responses, and scalar covariates.
- `estimate_be_zime()`: builds the B-spline basis and corrects both zero
  inflation and measurement error.
- `fit_joint_quantile_regression()`: fits the second-stage joint quantile
  regression with `quantregGrowth::gcrq`, matching the approach used in the
  simulation code.
