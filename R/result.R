# --- Result extractors (design doc 4.8) ------------------------------------

#' @importFrom stats weights
#' @include fit.R
NULL

# Kish effective sample size: (sum w)^2 / sum(w^2).
kish_ess <- function(w) sum(w)^2 / sum(w^2)

# weights(fit): the calibrated weights, in input row order. Registered on the
# stats::weights S3 generic (like print, via S7::methods_register in .onLoad).
S7::method(weights, raking_fit) <- function(object, ...) {
  object@weights
}

# summary(fit): balance table + convergence + weight diagnostics (doc 4.8, 7);
# returns the fit invisibly.
S7::method(summary, raking_fit) <- function(object, ...) {
  d <- object@diagnostics
  conv <- object@convergence
  cat("Raking fit summary\n")
  cat("  weight column: ", object@name, "\n", sep = "")
  cat("  units (n):     ", length(object@weights), "\n", sep = "")
  cat("  converged:     ", conv$converged, " (", conv$iterations,
      " iters)\n", sep = "")
  cat("  weight range:  ",
      paste(round(d$weight_range, 3), collapse = " to "), "\n", sep = "")
  cat("  Kish deff:     ", round(d$kish_deff, 3), "\n", sep = "")
  cat("  Kish ESS:      ", round(d$kish_ess, 2), "\n", sep = "")
  cat("  max |diff|:    ", signif(d$max_abs_diff, 4), "\n", sep = "")
  cat("\n  Balance (achieved vs target):\n")
  b <- object@balance
  b$target     <- round(b$target, 4)
  b$achieved   <- round(b$achieved, 4)
  b$difference <- round(b$difference, 4)
  print(b, row.names = FALSE)
  invisible(object)
}
