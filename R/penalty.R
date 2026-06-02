# --- Penalty constructors --------------------------------------------------
# User-facing constructors for the penalty() pipe step (design doc 4.4).
# Each returns an S7 regularizer object the solver already knows how to use.

#' Penalty constructors
#'
#' Build a weight penalty r(w) for a raking model. Each returns an S7
#' regularizer object, used as the `reg` argument to `penalty()`.
#'
#' @name penalties
NULL

#' @describeIn penalties Entropy penalty: pulls weights toward uniform.
#' @export
entropy <- function() {
  EntropyRegularizer()
}

#' @describeIn penalties KL penalty: pulls weights toward a positive reference
#'   distribution `prior`.
#' @param prior A positive numeric vector (the reference distribution).
#' @export
kl <- function(prior) {
  KLRegularizer(prior = prior)
}

#' @describeIn penalties Sum-of-squares penalty: ridge-style shrinkage.
#' @export
sum_squares <- function() {
  SumSquaresRegularizer()
}
