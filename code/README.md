# ODSkellam core functions

`ODSkellam.R` is a small, dependency-free implementation of the model layer of
the over-dispersed Skellam (ODS) reserving model. It is intended to be sourced
by users who want to build their own fitting, calibration, and reserving
workflow.

The script does not select a chain-ladder model, estimate dispersion, construct
prediction intervals, or run a simulation study.

## Model parameterisation

For accident period $i$ and development period $j$, define two non-negative
component intensities

$$
\lambda_{ij}^{+}=\alpha_i\beta_j^{+},
\qquad
\lambda_{ij}^{-}=\alpha_i\beta_j^{-}.
$$

The incremental ODS cell is the difference of two independent negative-binomial
variables,

$$
X_{ij}=Y_{ij}^{+}-Y_{ij}^{-},
\qquad
Y_{ij}^{\pm}\sim\mathrm{NB}\left(
\frac{\lambda_{ij}^{\pm}}{\phi-1},
\frac{1}{\phi}\right),
$$

where $\phi > 1$. Consequently,

$$
\mathbb{E}[X_{ij}]=\lambda_{ij}^{+}-\lambda_{ij}^{-},
\qquad
\mathrm{Var}(X_{ij})=\phi(\lambda_{ij}^{+}+\lambda_{ij}^{-}).
$$

All simulated values are integers in the monetary unit chosen by the user.
The value of `phi` is therefore unit-dependent.

## Quick start

```r
source("ODSkellam.R")

alpha <- c(100, 120)
beta_plus <- c(0.8, 0.3)
beta_minus <- c(0.0, 0.1)
phi <- 2

ods_moments(alpha, beta_plus, beta_minus, phi)

set.seed(2026)
cell_draw <- r_ods(alpha, beta_plus, beta_minus, phi)
cell_draw
```

`alpha` contains accident-period scales. `beta_plus` and `beta_minus` contain
development-period component shares and must have the same length. Both signed
components may be positive in the general model; the CL constructor below
uses the one-sided decomposition implied by each signed development share.

## Functions

### `ods_validate_parameters()`

Checks that `alpha` is strictly positive, both beta vectors are finite and
non-negative, and `phi` is finite and strictly greater than one.

### `ods_parameters()`

Builds the cell-level component intensities and negative-binomial parameters.
The returned list contains `lambda_plus`, `lambda_minus`, `size_plus`,
`size_minus`, and `prob`, together with the supplied model parameters.

### `ods_from_cl()`

Maps completed chain-ladder ultimates and selected development factors to the
ODS parameterisation:

```r
model <- ods_from_cl(
  ultimate = completed_ultimates,
  development_factors = selected_factors
)
```

`development_factors` must contain only the factors connecting adjacent
development periods. Some chain-ladder implementations return an additional
terminal factor equal to `1`; remove it before calling this function. The
returned list contains `alpha`, `beta_plus`, `beta_minus`, and
`development_share`.

### `ods_moments()`

Returns cell-level `mean` and `variance` matrices.

### `r_ods()`

Simulates independent ODS cells.

```r
one_draw <- r_ods(alpha, beta_plus, beta_minus, phi)
many_draws <- r_ods(alpha, beta_plus, beta_minus, phi, n = 1000L)
```

With `n = 1`, the result is an accident-period by development-period matrix.
With `n > 1`, the result is a three-dimensional array whose third dimension is
`replication`.

### `ods_reserve_moments()`

Returns the mean and variance of the aggregate over cells selected by a logical
`future_mask` matrix. The mask must have the same dimensions as the ODS cell
surface.

### `r_ods_reserve()`

Simulates the aggregate reserve directly using the equivalent sum of the two
negative-binomial component totals:

```r
reserve_draws <- r_ods_reserve(
  n = 10000L,
  alpha = model$alpha,
  beta_plus = model$beta_plus,
  beta_minus = model$beta_minus,
  phi = phi,
  future_mask = future_mask
)
```

This is faster than simulating every future cell when only the aggregate reserve
is required.

## Example with `ChainLadder::RAA`

The `RAA` triangle is included in the `ChainLadder` package. It is a positive
cumulative triangle, so this example produces the one-sided ODP special case of
ODS. An incurred triangle with negative development can produce non-zero
`beta_minus` values.

```r
library(ChainLadder)
source("odskellam.R")

# Select a conventional Mack chain-ladder completion.
mack <- MackChainLadder(
  RAA,
  est.sigma = "Mack",
  mse.method = "Mack"
)

# MackChainLadder() returns a trailing terminal factor equal to 1.
development_factors <- as.numeric(
  mack$f[seq_len(ncol(RAA) - 1L)]
)
ultimate <- as.numeric(mack$FullTriangle[, ncol(RAA)])

model <- ods_from_cl(
  ultimate = ultimate,
  development_factors = development_factors
)

# For a 10 by 10 triangle, cells with i + j > 11 are future cells when
# rows and columns are numbered from 1.
future_mask <- outer(
  seq_len(nrow(RAA)),
  seq_len(ncol(RAA)),
  "+"
) > nrow(RAA) + 1L

phi <- 2 # Illustrative only; calibrate phi for an applied analysis.

reserve_moments <- ods_reserve_moments(
  alpha = model$alpha,
  beta_plus = model$beta_plus,
  beta_minus = model$beta_minus,
  phi = phi,
  future_mask = future_mask
)
reserve_moments

set.seed(2026)
reserve_draws <- r_ods_reserve(
  n = 10000L,
  alpha = model$alpha,
  beta_plus = model$beta_plus,
  beta_minus = model$beta_minus,
  phi = phi,
  future_mask = future_mask
)

mean(reserve_draws)
quantile(reserve_draws, probs = c(0.025, 0.5, 0.975))

# A single simulated full cell surface, including observed and future positions.
simulated_cells <- r_ods(
  alpha = model$alpha,
  beta_plus = model$beta_plus,
  beta_minus = model$beta_minus,
  phi = phi
)
dim(simulated_cells)
```

The fitted chain-ladder quantities provide the deterministic mean anchor only.
The example's value `phi = 2` is deliberately illustrative; an applied workflow
should estimate or calibrate dispersion separately and document the monetary
unit used.