# --- The weighting model object (design doc 5) -----------------------------
# A data-free S7 object describing a weighting problem. Built up with pipe
# steps (target/calibrate/penalty/bounds) and then fit() to data.

raking_model <- S7::new_class(
  "raking_model",
  package = "raking",
  properties = list(
    loss = S7::new_property(S7::class_character, default = "KL"),
    name = S7::new_property(S7::class_character, default = "wt"),
    targets = S7::new_property(S7::class_list, default = list()),
    calibrations = S7::new_property(S7::class_list, default = list()),
    penalty = S7::new_property(S7::class_any, default = NULL),
    lambda = S7::new_property(S7::class_numeric, default = 1),
    bounds = S7::new_property(S7::class_any, default = NULL)
  ),
  validator = function(self) {
    valid_losses <- c("KL", "L1", "L2", "chi_square")
    if (length(self@loss) != 1 || !self@loss %in% valid_losses) {
      return(paste0("`loss` must be one of: ",
                    paste(valid_losses, collapse = ", ")))
    }
    if (length(self@name) != 1) {
      return("`name` must be a single string")
    }
    for (t in self@targets) {
      if (!is.null(t$tol) && t$tol < 0) {
        return("each target `tol` must be >= 0")
      }
    }
    b <- self@bounds
    if (!is.null(b)) {
      if (is.null(b$min) || is.null(b$max)) {
        return("`bounds` must have `min` and `max`")
      }
      if (b$min >= b$max) return("`bounds$min` must be < `bounds$max`")
      # Doc 4.5: min in [0, 1), max > 1. The mean weight is 1 by construction,
      # so min <= 1 <= max is required for the projection to be feasible.
      if (b$min < 0 || b$min >= 1) {
        return("`bounds$min` must be in [0, 1)")
      }
      if (b$max <= 1) return("`bounds$max` must be > 1")
    }
    NULL
  }
)

#' Create a raking model
#'
#' Starts a data-free weighting model. Extend it with pipe steps
#' (`penalty()`, `bounds()`, and later `target()`/`calibrate()`), then `fit()`
#' it to data.
#'
#' @param loss Discrepancy measure: one of "KL", "L1", "L2", "chi_square".
#'   Default "KL".
#' @param name Output weight-column name. Default "wt".
#' @return A `raking_model` object.
#' @examples
#' data <- data.frame(sex = c("M", "F", "M", "F"))
#' model <- raking() |>
#'   target(data.frame(sex = c("M", "F"), p = c(0.5, 0.5))) |>
#'   penalty(entropy())
#' fit <- model |> fit(data)
#' weights(fit)
#' summary(fit)
#' @export
raking <- function(loss = "KL", name = "wt") {
  raking_model(loss = loss, name = name)
}

# A short summary when a model is printed.
S7::method(print, raking_model) <- function(x, ...) {
  cat("<raking_model>\n")
  cat("  loss:        ", x@loss, "\n", sep = "")
  cat("  weight name: ", x@name, "\n", sep = "")
  cat("  targets:     ", length(x@targets), "\n", sep = "")
  cat("  calibrations:", length(x@calibrations), "\n", sep = "")
  cat("  penalty:     ", if (is.null(x@penalty)) "none" else "set", "\n", sep = "")
  b <- if (is.null(x@bounds)) list(min = 0, max = Inf) else x@bounds
  cat("  bounds:      [", b$min, ", ", b$max, "]\n", sep = "")
  invisible(x)
}
