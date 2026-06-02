# --- Solver control settings (design doc 4.7) ------------------------------

#' Solver control settings
#'
#' @param maxiter Maximum ADMM iterations.
#' @param rho ADMM penalty parameter.
#' @param margin_tol Margin-based convergence tolerance. Internally converted to
#'   a raw ADMM tolerance scaled by problem size: `eps = margin_tol / sqrt(m + 2n)`,
#'   giving roughly constant margin accuracy regardless of size. Set to `NULL`
#'   to use raw `eps_abs`/`eps_rel`. Supplying `eps_abs`/`eps_rel` without
#'   `margin_tol` sets it to `NULL` automatically.
#' @param eps_abs,eps_rel Raw ADMM residual tolerances (used when
#'   `margin_tol = NULL`).
#' @param verbose Print residuals while solving.
#' @return A `rake_control` list of settings.
#' @examples
#' rake_control(maxiter = 1000)
#' # use raw residual tolerances instead of margin_tol:
#' rake_control(eps_abs = 1e-7, eps_rel = 1e-7)
#' @export
rake_control <- function(maxiter = 5000, rho = 50, margin_tol = 1e-4,
                         eps_abs = 1e-5, eps_rel = 1e-5, verbose = FALSE) {
  named <- names(match.call())[-1]
  if (!("margin_tol" %in% named) &&
      ("eps_abs" %in% named || "eps_rel" %in% named)) {
    margin_tol <- NULL
  }
  structure(
    list(maxiter = maxiter, rho = rho, margin_tol = margin_tol,
         eps_abs = eps_abs, eps_rel = eps_rel, verbose = verbose),
    class = "rake_control"
  )
}

# Merge a control list (rake_control or plain list) with defaults and resolve
# the margin_tol size-scaling against problem dimensions m, n. Internal.
resolve_control <- function(control, m, n) {
  if (is.null(control)) control <- list()
  get <- function(field, default) {
    if (field %in% names(control)) control[[field]] else default
  }
  margin_tol <- get("margin_tol", 1e-4)
  # A plain list giving raw eps but no margin_tol means "use raw eps".
  has_eps <- "eps_abs" %in% names(control) || "eps_rel" %in% names(control)
  if (!("margin_tol" %in% names(control)) && has_eps) margin_tol <- NULL

  eps_abs <- get("eps_abs", 1e-5)
  eps_rel <- get("eps_rel", 1e-5)
  if (!is.null(margin_tol)) {
    eps_abs <- eps_rel <- margin_tol / sqrt(m + 2 * n)
  }
  maxiter <- get("maxiter", 5000)
  if (!is.numeric(maxiter) || length(maxiter) != 1 || maxiter < 1) {
    stop("`maxiter` must be a positive integer", call. = FALSE)
  }
  rho <- get("rho", 50)
  if (!is.numeric(rho) || length(rho) != 1 || rho <= 0) {
    stop("`rho` must be a positive number", call. = FALSE)
  }
  list(
    maxiter = maxiter,
    rho     = rho,
    eps_abs = eps_abs,
    eps_rel = eps_rel,
    verbose = get("verbose", FALSE)
  )
}
