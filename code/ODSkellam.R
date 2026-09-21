# ODSkellam: core functions for a signed, over-dispersed Skellam model
#
# The script is deliberately dependency-free. It provides the model layer only:
# parameter construction, cell moments, independent cell simulation, and direct
# simulation of an aggregate future reserve. Fitting and triangle selection are
# left to the user.
#
# All quantities are expressed in the chosen integer monetary unit. The
# dispersion parameter must therefore be interpreted in that same unit.

ods_validate_parameters <- function(alpha, beta_plus, beta_minus, phi) {
  if (!is.numeric(alpha) || length(alpha) == 0L ||
      any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop("alpha must contain finite, strictly positive values.")
  }
  if (!is.numeric(beta_plus) || length(beta_plus) == 0L ||
      any(!is.finite(beta_plus)) || any(beta_plus < 0)) {
    stop("beta_plus must contain finite, non-negative values.")
  }
  if (!is.numeric(beta_minus) || length(beta_minus) != length(beta_plus) ||
      any(!is.finite(beta_minus)) || any(beta_minus < 0)) {
    stop("beta_minus must match beta_plus and contain finite, non-negative values.")
  }
  if (length(phi) != 1L || !is.numeric(phi) || !is.finite(phi) || phi <= 1) {
    stop("phi must be a finite scalar greater than one.")
  }
  invisible(TRUE)
}

ods_parameters <- function(alpha, beta_plus, beta_minus, phi) {
  ods_validate_parameters(alpha, beta_plus, beta_minus, phi)

  lambda_plus <- alpha %o% beta_plus
  lambda_minus <- alpha %o% beta_minus

  structure(
    list(
      alpha = alpha,
      beta_plus = beta_plus,
      beta_minus = beta_minus,
      phi = phi,
      lambda_plus = lambda_plus,
      lambda_minus = lambda_minus,
      size_plus = lambda_plus / (phi - 1),
      size_minus = lambda_minus / (phi - 1),
      prob = 1 / phi
    ),
    class = "ods_parameters"
  )
}

ods_moments <- function(alpha, beta_plus, beta_minus, phi) {
  parameters <- ods_parameters(alpha, beta_plus, beta_minus, phi)

  list(
    mean = parameters$lambda_plus - parameters$lambda_minus,
    variance = phi * (parameters$lambda_plus + parameters$lambda_minus)
  )
}

ods_from_cl <- function(ultimate, development_factors) {
  if (!is.numeric(ultimate) || length(ultimate) == 0L ||
      any(!is.finite(ultimate)) || any(ultimate <= 0)) {
    stop("ultimate must contain finite, strictly positive values.")
  }
  if (!is.numeric(development_factors) || length(development_factors) == 0L ||
      any(!is.finite(development_factors)) || any(development_factors == 0)) {
    stop("development_factors must be finite and non-zero.")
  }

  factor_products <- rev(cumprod(rev(development_factors)))
  development_share <- c(
    1 / prod(development_factors),
    (development_factors - 1) / factor_products
  )
  if (any(!is.finite(development_share))) {
    stop("The development factors do not produce finite ODS shares.")
  }

  list(
    alpha = ultimate,
    beta_plus = pmax(development_share, 0),
    beta_minus = pmax(-development_share, 0),
    development_share = development_share
  )
}

.ods_rnbinom <- function(n, lambda, phi) {
  draws <- matrix(0, nrow = n, ncol = length(lambda))
  active <- which(lambda > 0)

  if (length(active) > 0L) {
    draws[, active] <- matrix(
      rnbinom(
        n * length(active),
        size = rep(lambda[active] / (phi - 1), each = n),
        prob = 1 / phi
      ),
      nrow = n,
      ncol = length(active)
    )
  }
  draws
}

r_ods <- function(alpha, beta_plus, beta_minus, phi, n = 1L) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(n)
  parameters <- ods_parameters(alpha, beta_plus, beta_minus, phi)

  plus <- .ods_rnbinom(n, as.vector(parameters$lambda_plus), phi)
  minus <- .ods_rnbinom(n, as.vector(parameters$lambda_minus), phi)
  draws <- plus - minus

  if (n == 1L) {
    return(matrix(
      draws[1L, ],
      nrow = length(alpha),
      ncol = length(beta_plus),
      dimnames = dimnames(parameters$lambda_plus)
    ))
  }

  array(
    t(draws),
    dim = c(length(alpha), length(beta_plus), n),
    dimnames = list(
      accident_period = names(alpha),
      development_period = names(beta_plus),
      replication = seq_len(n)
    )
  )
}

.ods_validate_future_mask <- function(future_mask, dimensions) {
  if (!is.matrix(future_mask) || !is.logical(future_mask) ||
      !identical(dim(future_mask), dimensions)) {
    stop("future_mask must be a logical matrix matching the ODS cell dimensions.")
  }
  invisible(TRUE)
}

ods_reserve_moments <- function(alpha, beta_plus, beta_minus, phi, future_mask) {
  parameters <- ods_parameters(alpha, beta_plus, beta_minus, phi)
  .ods_validate_future_mask(future_mask, dim(parameters$lambda_plus))

  lambda_plus <- sum(parameters$lambda_plus[future_mask])
  lambda_minus <- sum(parameters$lambda_minus[future_mask])

  list(
    mean = lambda_plus - lambda_minus,
    variance = phi * (lambda_plus + lambda_minus),
    positive_mean = lambda_plus,
    negative_mean = lambda_minus
  )
}

r_ods_reserve <- function(n, alpha, beta_plus, beta_minus, phi, future_mask) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(n)
  parameters <- ods_parameters(alpha, beta_plus, beta_minus, phi)
  .ods_validate_future_mask(future_mask, dim(parameters$lambda_plus))

  lambda_plus <- sum(parameters$lambda_plus[future_mask])
  lambda_minus <- sum(parameters$lambda_minus[future_mask])
  plus <- .ods_rnbinom(n, lambda_plus, phi)
  minus <- .ods_rnbinom(n, lambda_minus, phi)

  as.numeric(plus[, 1L] - minus[, 1L])
}
