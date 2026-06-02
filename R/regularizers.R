#' Regularizers
#'
#' Regularizers keep the optimized weights "regular" (close to uniform, not
#' extreme). Each is an S7 class implementing [prox()].
#'
#' @param limit Optional upper cap (> 1) on how far a weight may stray from
#'   uniform (`EntropyRegularizer`, `KLRegularizer`); default `Inf`.
#' @param prior Positive numeric reference distribution (`KLRegularizer`).
#' @param base Numeric base/target weights for `SumSquaresRegularizer`; an empty
#'   vector (default) shrinks toward zero.
#' @name regularizers
NULL

#' @describeIn regularizers Identity regularizer: no penalty, prox is a no-op.
#' @export
ZeroRegularizer <- S7::new_class("ZeroRegularizer", package = "raking")

S7::method(prox, ZeroRegularizer) <- function(object, v, lam, ...) {
  v
}

# evaluate: an unpenalized model contributes nothing to the objective.
S7::method(evaluate, ZeroRegularizer) <- function(object, v) 0

#' @describeIn regularizers Negative-entropy regularizer. Pushes weights toward
#'   uniform. Optional `limit` (> 1) caps how far any weight can stray from
#'   uniform; the default `Inf` means no cap.
#' @export
EntropyRegularizer <- S7::new_class(
  "EntropyRegularizer",
  package = "raking",
  properties = list(
    limit = S7::new_property(S7::class_numeric, default = Inf)
  ),
  validator = function(self) {
    if (length(self@limit) != 1) return("`limit` must be a single number")
    if (is.finite(self@limit) && self@limit <= 1) {
      return("`limit` must be greater than 1")
    }
    NULL
  }
)

S7::method(prox, EntropyRegularizer) <- function(object, v, lam, ...) {
  what <- entropy_prox(v, lam)
  n <- length(v)
  pmin(pmax(what, 1 / (object@limit * n)), object@limit / n)
}

# evaluate: negative entropy relative to uniform b = 1/n, i.e. sum(w log(w/b)).
S7::method(evaluate, EntropyRegularizer) <- function(object, v) {
  n <- length(v)
  sum(v * log(v * n))
}

#' @describeIn regularizers KL regularizer: pulls weights toward a given
#'   `prior` distribution rather than uniform.
#' @export
KLRegularizer <- S7::new_class(
  "KLRegularizer",
  package = "raking",
  properties = list(
    prior = S7::class_numeric,
    limit = S7::new_property(S7::class_numeric, default = Inf)
  ),
  validator = function(self) {
    if (any(self@prior <= 0)) return("`prior` values must all be positive")
    if (length(self@limit) != 1) return("`limit` must be a single number")
    if (is.finite(self@limit) && self@limit <= 1) {
      return("`limit` must be greater than 1")
    }
    NULL
  }
)

S7::method(prox, KLRegularizer) <- function(object, v, lam, ...) {
  what <- entropy_prox(v + lam * log(object@prior), lam)
  n <- length(v)
  pmin(pmax(what, 1 / (object@limit * n)), object@limit / n)
}

# evaluate: KL divergence of w from the prior, sum(w log(w/prior)).
S7::method(evaluate, KLRegularizer) <- function(object, v) {
  sum(v * log(v / object@prior))
}

#' @describeIn regularizers Sum-of-squares regularizer: squared deviation from
#'   the base weights `base` (doc 6.3). With the default empty `base` it shrinks
#'   toward zero; `fit()` sets `base` to the (uniform or design) base weights.
#' @export
SumSquaresRegularizer <- S7::new_class(
  "SumSquaresRegularizer",
  package = "raking",
  properties = list(
    base = S7::new_property(S7::class_numeric, default = numeric(0))
  )
)

S7::method(prox, SumSquaresRegularizer) <- function(object, v, lam, ...) {
  if (length(object@base) == 0) {
    v / (1 + 2 * lam)
  } else {
    (v + 2 * lam * object@base) / (1 + 2 * lam)
  }
}

# evaluate: sum of squared deviations from base (or from zero if no base set).
S7::method(evaluate, SumSquaresRegularizer) <- function(object, v) {
  if (length(object@base) == 0) sum(v^2) else sum((v - object@base)^2)
}
