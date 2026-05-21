# Generate latent curves, zero-inflated replicates, scalar covariates, and response.
generate_zime_sample <- function(
  n = 100,
  grid_length = 100,
  n_replicates = 7,
  n_fourier_basis = 50,
  p_nonzero = NULL,
  n_segments = 1,
  nonzero_prob = 0.7,
  variation = 0,
  first_segment_prob = 0.4,
  pattern = c("constant", "alternating"),
  eta_beta = 0,
  eta_eps = 1,
  sigma = 0.1,
  gamma = c(0.5, 0.6),
  seed = NULL
) {
  pattern <- match.arg(pattern)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  grid <- seq(0, 1, length.out = grid_length)
  if (is.null(p_nonzero)) {
    if (n_segments == 1) {
      lower <- max(nonzero_prob - variation, 0)
      upper <- min(nonzero_prob + variation, 1)
      p_nonzero <- matrix(runif(n, lower, upper), nrow = 1)
    } else if (pattern == "constant") {
      p_nonzero <- matrix(rep(nonzero_prob, n_segments * n), nrow = n_segments)
    } else if (pattern == "alternating") {
      segment_prob <- rep(c(first_segment_prob, 1 - first_segment_prob), length.out = n_segments)
      p_nonzero <- matrix(rep(segment_prob, n), nrow = n_segments)
    }
  }
  if (!is.matrix(p_nonzero) || ncol(p_nonzero) != n) {
    stop("p_nonzero must be an n_segments x n matrix.", call. = FALSE)
  }

  n_segments <- nrow(p_nonzero)
  basis <- matrix(
    unlist(lapply(seq_len(grid_length), function(i) {
      c(1, sqrt(2) * cos(seq_len(n_fourier_basis - 1) * pi * grid[i]))
    })),
    ncol = grid_length
  )
  loadings <- (-1)^(seq_len(n_fourier_basis) + 1) * seq_len(n_fourier_basis)^(-1)
  x <- 5 + matrix(runif(n * n_fourier_basis, -sqrt(3), sqrt(3)), ncol = n_fourier_basis) %*%
    (basis * loadings)
  x <- pmax(x, 1e-8)

  if (n_segments == 1) {
    segments <- list(seq_len(grid_length))
  } else {
    segments <- split(seq_len(grid_length), cut(seq_len(grid_length), breaks = n_segments, labels = FALSE))
  }
  w <- array(0, dim = c(n, grid_length, n_replicates))
  for (i in seq_len(n)) {
    p_i <- numeric(grid_length)
    for (m in seq_along(segments)) {
      p_i[segments[[m]]] <- p_nonzero[m, i]
    }
    for (j in seq_len(n_replicates)) {
      w[i, , j] <- stats::rpois(grid_length, x[i, ]) *
        stats::rbinom(grid_length, size = 1, prob = p_i)
    }
  }

  z_continuous <- stats::rnorm(n, mean = 1, sd = 0.5)
  z_binary <- stats::rbinom(n, size = 1, prob = 0.6)
  z <- cbind(z_continuous = z_continuous, z_binary = z_binary)
  beta0 <- 0.5 * sin(pi * grid)
  h <- function(tau) 0.2 * (tau - 0.5)
  u <- stats::runif(n)
  beta_u <- (1 + eta_beta * h(u)) %o% beta0
  mu <- rowMeans(x * beta_u) + gamma[1] * z_continuous + gamma[2] * z_binary
  sigma_u <- pmax(sigma * (1 + eta_eps * h(u)), 1e-8)
  y <- mu + sigma_u * stats::qnorm(u)

  list(
    x = x,
    w = w,
    y = y,
    z = z,
    grid = grid,
    p_nonzero = t(p_nonzero),
    zero_inflation = 1 - t(p_nonzero),
    u = u,
    beta_used = beta_u,
    beta0 = beta0
  )
}

# Estimate latent curves from projected repeated measurements using lme4.
.estimate_latent_curves_with_basis <- function(w, basis_object) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Install package 'lme4' to run the BE-ZIME correction.", call. = FALSE)
  }

  scores <- apply(w, MARGIN = c(1, 3), function(x) basis_object$projection %*% x)
  n <- dim(scores)[2]
  n_basis <- dim(scores)[1]
  n_replicates <- dim(scores)[3]
  x_scores <- matrix(NA_real_, nrow = n, ncol = n_basis)

  for (k in seq_len(n_basis)) {
    wk <- scores[k, , ]
    data_k <- data.frame(
      wk = as.vector(wk),
      seqn = rep(seq_len(n), n_replicates)
    )
    model <- lme4::lmer(wk ~ (1 | seqn), data = data_k)
    if (lme4::isSingular(model)) {
      x_scores[, k] <- rowMeans(wk)
    } else {
      x_scores[, k] <- stats::predict(model)[seq_len(n)]
    }
  }

  x_hat <- x_scores %*% t(basis_object$basis)
  list(x_hat = x_hat, x_scores = x_scores, basis = basis_object)
}

# Estimate latent curves after correcting both zero inflation and measurement error.
estimate_be_zime <- function(
  w,
  grid = NULL,
  basis_df = 6,
  basis_degree = 3,
  n_segments = 1,
  max_iter = 10,
  tolerance = 1e-3,
  step_length = 0.1,
  min_prob = 1e-6
) {
  n <- dim(w)[1]
  grid_length <- dim(w)[2]
  n_replicates <- dim(w)[3]
  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = grid_length)
  }
  basis <- as.matrix(splines::bs(grid, df = basis_df, degree = basis_degree, intercept = TRUE))
  basis_object <- list(
    basis = basis,
    projection = MASS::ginv(crossprod(basis)) %*% t(basis),
    grid = grid,
    df = basis_df,
    degree = basis_degree
  )

  if (n_segments == 1) {
    segments <- list(seq_len(grid_length))
  } else {
    segments <- split(seq_len(grid_length), cut(seq_len(grid_length), breaks = n_segments, labels = FALSE))
  }

  p_nonzero <- matrix(NA_real_, nrow = n, ncol = n_segments)
  for (m in seq_along(segments)) {
    p_nonzero[, m] <- rowMeans(apply(w[, segments[[m]], , drop = FALSE] > 0, c(1, 3), mean))
  }
  p_nonzero <- matrix(pmin(1 - min_prob, pmax(min_prob, p_nonzero)), nrow = n, ncol = n_segments)

  converged <- FALSE
  x_fit <- NULL
  for (iter in seq_len(max_iter)) {
    w_star <- w
    for (j in seq_len(n_replicates)) {
      for (m in seq_along(segments)) {
        w_star[, segments[[m]], j] <- w[, segments[[m]], j] / p_nonzero[, m]
      }
    }

    x_fit <- .estimate_latent_curves_with_basis(w = w_star, basis_object = basis_object)
    x_hat <- pmax(x_fit$x_hat, 1e-8)
    lambda <- array(rep(as.vector(x_hat), n_replicates), dim = c(n, grid_length, n_replicates))
    density <- pmax(stats::dpois(as.vector(w), as.vector(lambda)), 1e-50)
    density <- array(density, dim = c(n, grid_length, n_replicates))
    observed_positive <- w > 0

    p_new <- p_nonzero
    for (i in seq_len(n)) {
      for (m in seq_along(segments)) {
        zero_prob <- 1 - p_nonzero[i, m]
        idx <- segments[[m]]
        gi <- as.vector(density[i, idx, ])
        ypos <- as.numeric(observed_positive[i, idx, ])

        for (inner in seq_len(max_iter)) {
          denom <- zero_prob + (1 - zero_prob) * gi
          deriv1 <- (1 - ypos) * (1 - gi) / denom - ypos / (1 - zero_prob)
          deriv2 <- -(1 - ypos) * ((1 - gi) / denom)^2 - ypos / (1 - zero_prob)^2
          step <- sum(deriv1) / sum(deriv2)
          if (!is.finite(step)) {
            break
          }
          zero_next <- zero_prob - step_length * step
          zero_next <- pmin(1 - min_prob, pmax(min_prob, zero_next))
          if (abs(zero_next - zero_prob) < tolerance) {
            zero_prob <- zero_next
            break
          }
          zero_prob <- zero_next
        }
        p_new[i, m] <- 1 - zero_prob
      }
    }

    p_new <- matrix(pmin(1 - min_prob, pmax(min_prob, p_new)), nrow = n, ncol = n_segments)
    if (mean(abs(p_new - p_nonzero)) < tolerance) {
      p_nonzero <- p_new
      converged <- TRUE
      break
    }
    p_nonzero <- p_new
  }

  list(
    x_hat = x_fit$x_hat,
    x_scores = x_fit$x_scores,
    p_nonzero_hat = p_nonzero,
    zero_inflation_hat = 1 - p_nonzero,
    basis = basis_object,
    iterations = iter,
    converged = converged
  )
}

# Fit joint scalar-on-function quantile regression using quantregGrowth::gcrq.
fit_joint_quantile_regression <- function(y, x_hat, z = NULL, basis_object, taus = c(0.25, 0.5, 0.75)) {
  if (!requireNamespace("quantregGrowth", quietly = TRUE)) {
    stop("Install package 'quantregGrowth' to fit joint quantile regression.", call. = FALSE)
  }

  x_scores <- as.matrix(x_hat) %*% t(basis_object$projection)
  if (is.null(z)) {
    model <- quantregGrowth::gcrq(y ~ as.matrix(x_scores), tau = taus)
  } else {
    z <- as.matrix(z)
    model <- quantregGrowth::gcrq(y ~ as.matrix(x_scores) + z, tau = taus)
  }

  residuals <- as.matrix(model$residuals)
  losses <- sapply(seq_along(taus), function(j) {
    mean(residuals[, j] * (taus[j] - (residuals[, j] < 0)))
  })

  coefficients <- as.matrix(model$coefficients)
  n_basis <- ncol(x_scores)
  beta_scores <- coefficients[seq_len(n_basis) + 1, , drop = FALSE]
  beta_grid <- t(basis_object$projection) %*% beta_scores * ncol(x_hat)

  list(
    model = model,
    x_scores = x_scores,
    coefficients = coefficients,
    beta_grid = beta_grid,
    check_loss = losses,
    taus = taus
  )
}
