#' Loss functions
#'
#' Each loss is an S7 class describing one constraint: a `target` to hit and a
#' proximal operator the solver calls each iteration. We build them up one at a
#' time in this file.
#'
#' @param target Numeric vector of target values the constraint should match.
#' @param diag_weight Numeric scalar or vector of per-element weights
#'   (`LeastSquaresLoss`; default 1).
#' @param lower,upper Numeric vectors giving the lower/upper offsets of the
#'   allowed band around `target` (`InequalityLoss`).
#' @param scale Numeric scale factor on the KL term (`KLLoss`; default 0.5).
#' @param inner Inner soft loss object that a `BandedLoss` wraps and clips into
#'   the band.
#' @name losses
NULL

#' @describeIn losses Exact-match loss. Drives the weighted statistic to equal
#'   the target exactly; its proximal operator simply returns the target.
#' @export
EqualityLoss <- S7::new_class(
  "EqualityLoss",
  package = "raking",
  properties = list(
    target = S7::class_numeric
  )
)

# prox for EqualityLoss: exact matching snaps straight to the target, so it
# ignores the current vector and the step size.
S7::method(prox, EqualityLoss) <- function(object, v, lam, ...) {
  object@target
}

# evaluate: hard equality constraint -- 0 when matched, +Inf otherwise. Used by
# the solver's best-iterate objective tracking.
S7::method(evaluate, EqualityLoss) <- function(object, v) {
  if (isTRUE(all.equal(as.numeric(v), object@target, tolerance = 1e-8))) {
    0
  } else {
    Inf
  }
}

#' @describeIn losses Least-squares loss. Penalizes squared deviation from the
#'   target, optionally weighted per element by `diag_weight`.
#' @export
LeastSquaresLoss <- S7::new_class(
  "LeastSquaresLoss",
  package = "raking",
  properties = list(
    target = S7::class_numeric,
    diag_weight = S7::new_property(S7::class_numeric, default = 1)
  )
)

# prox: minimizes weighted squared error plus the proximal term. lam is the
# step size (the paper's lam = 1/rho).
S7::method(prox, LeastSquaresLoss) <- function(object, v, lam, ...) {
  dw2 <- object@diag_weight^2
  (dw2 * object@target + v / lam) / (dw2 + 1 / lam)
}

# evaluate: the actual objective value, sum of weighted squared deviations.
S7::method(evaluate, LeastSquaresLoss) <- function(object, v) {
  sum((object@diag_weight * (v - object@target))^2)
}

#' @describeIn losses Inequality (range) loss. Allows the weighted statistic to
#'   sit anywhere in `[target + lower, target + upper]` at no cost; its prox
#'   clips into that band.
#' @export
InequalityLoss <- S7::new_class(
  "InequalityLoss",
  package = "raking",
  properties = list(
    target = S7::class_numeric,
    lower = S7::class_numeric,
    upper = S7::class_numeric
  ),
  validator = function(self) {
    if (length(self@lower) != length(self@upper)) {
      return("`lower` and `upper` must have the same length")
    }
    if (any(self@lower > self@upper)) {
      return("each `lower` must be <= the matching `upper`")
    }
    NULL
  }
)

# prox: clip (v - target) into [lower, upper], then shift back by target.
S7::method(prox, InequalityLoss) <- function(object, v, lam, ...) {
  object@target + pmin(pmax(v - object@target, object@lower), object@upper)
}

# evaluate: hard range constraint -- 0 inside the band, +Inf outside. The
# solver uses this both for best-iterate tracking and the tol>0 convergence
# check (doc 6.4 item 2).
S7::method(evaluate, InequalityLoss) <- function(object, v) {
  d <- v - object@target
  if (all(d >= object@lower - 1e-9) && all(d <= object@upper + 1e-9)) {
    0
  } else {
    Inf
  }
}

# Internal helper: the "entropy prox" via the Lambert-W function. Shared by
# KLLoss here and (later) the entropy regularizer. lamW::lambertW0 is the real
# branch of Lambert W. Not exported.
entropy_prox <- function(f, lam) {
  lam * lamW::lambertW0(exp(f / lam - 1) / lam)
}

#' @describeIn losses Kullback-Leibler loss. Matches a target distribution in
#'   KL divergence; `scale` weights the term (default 0.5, as in regrake).
#' @export
KLLoss <- S7::new_class(
  "KLLoss",
  package = "raking",
  properties = list(
    target = S7::class_numeric,
    scale = S7::new_property(S7::class_numeric, default = 0.5)
  )
)

# prox: shift by lam*scale*log(target), then apply the entropy prox.
S7::method(prox, KLLoss) <- function(object, v, lam, ...) {
  s <- object@scale
  entropy_prox(v + lam * s * log(object@target), lam * s)
}

# evaluate: scale * sum of KL divergence terms x*log(x/target) - x + target.
S7::method(evaluate, KLLoss) <- function(object, v) {
  object@scale * sum(v * log(v / object@target) - v + object@target)
}

#' @describeIn losses L1 (absolute-deviation) loss. Robust matching; its
#'   proximal operator is soft-thresholding around the target.
#' @export
L1Loss <- S7::new_class(
  "L1Loss",
  package = "raking",
  properties = list(target = S7::class_numeric)
)

# prox: shift to target, soft-threshold by lam, shift back.
S7::method(prox, L1Loss) <- function(object, v, lam, ...) {
  z <- v - object@target
  object@target + sign(z) * pmax(abs(z) - lam, 0)
}

S7::method(evaluate, L1Loss) <- function(object, v) {
  sum(abs(v - object@target))
}

#' @describeIn losses Penalized-band loss. Wraps a soft `inner` loss (KL/L2/L1)
#'   and a hard band `[lower, upper]`: its prox runs the inner soft prox (which
#'   pulls toward the target) and then clips the result into the band. This is
#'   the loss `fit()` builds for a margin with `tol > 0`: the soft loss prefers
#'   the exact target, while the band is a hard wall the achieved value cannot
#'   cross. (A `tol = 0` margin uses `EqualityLoss` instead, since a zero-width
#'   band is just the target point.)
#' @export
BandedLoss <- S7::new_class(
  "BandedLoss",
  package = "raking",
  properties = list(
    inner  = S7::class_any,        # an S7 loss carrying the same target
    target = S7::class_numeric,
    lower  = S7::class_numeric,    # absolute lower band edge (target - tol)
    upper  = S7::class_numeric     # absolute upper band edge (target + tol)
  ),
  validator = function(self) {
    if (any(self@lower > self@upper)) {
      return("each `lower` must be <= the matching `upper`")
    }
    NULL
  }
)

# prox: soft pull toward the target, then clip into the hard band.
S7::method(prox, BandedLoss) <- function(object, v, lam, ...) {
  s <- prox(object@inner, v, lam)
  pmin(pmax(s, object@lower), object@upper)
}

# evaluate: the inner soft loss value (the band itself is enforced by the prox).
S7::method(evaluate, BandedLoss) <- function(object, v) {
  evaluate(object@inner, v)
}
