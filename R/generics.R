#' Proximal operator
#'
#' The proximal operator that every loss and regularizer must implement. The
#' ADMM solver calls `prox()` on each object, once per iteration.
#'
#' @param object A loss or regularizer object (an S7 object).
#' @param v Numeric vector the operator acts on.
#' @param lam Proximal step-size parameter (the paper's `lam`). The solver
#'   passes `1 / rho` for losses and `lambda / rho` for the regularizer.
#' @param ... Additional arguments used by specific methods.
#' @return A numeric vector the same length as `v`.
#' @export
prox <- S7::new_generic("prox", "object", function(object, v, lam, ...) {
  S7::S7_dispatch()
})

#' Evaluate an objective term
#'
#' Returns the scalar value a loss contributes to the objective function, given
#' the achieved values `v`. The solver uses this to track the best solution.
#'
#' @param object A loss object (an S7 object).
#' @param v Numeric vector of achieved (weighted) values.
#' @return A length-one numeric value.
#' @export
evaluate <- S7::new_generic("evaluate", "object", function(object, v) {
  S7::S7_dispatch()
})
