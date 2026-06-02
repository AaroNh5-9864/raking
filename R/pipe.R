# --- Pipe steps (design doc 4.2-4.5) ---------------------------------------
# Each takes a model and returns an updated model. All validation here is
# data-free; data-dependent checks wait for fit().

#' Add a categorical margin target
#'
#' @param model A `raking_model`.
#' @param df A data.frame with one or more variable columns plus exactly one of
#'   `p` (proportions) or `n` (counts).
#' @param tol Tolerance band half-width (>= 0). `tol = 0` means exact matching;
#'   `tol > 0` allows each level within +/- tol of its target.
#' @return The updated `raking_model`.
#' @examples
#' raking() |> target(data.frame(age = c("Y", "O"), p = c(0.4, 0.6)))
#' # joint margin with counts, matched within +/- 0.02:
#' joint <- data.frame(sex = c("M", "M", "F", "F"), age = c("Y", "O", "Y", "O"),
#'                     n = c(110, 140, 130, 120))
#' raking() |> target(joint, tol = 0.02)
#' @export
target <- function(model, df, tol = 0) {
  if (!is.data.frame(df)) stop("`df` must be a data.frame", call. = FALSE)
  has_p <- "p" %in% names(df)
  has_n <- "n" %in% names(df)
  if (has_p == has_n) {
    stop("`df` must have exactly one of `p` or `n`", call. = FALSE)
  }
  if (tol < 0) stop("`tol` must be >= 0", call. = FALSE)

  # At least one variable column, and those columns must be categorical.
  vars <- setdiff(names(df), c("p", "n"))
  if (length(vars) < 1) {
    stop("`df` needs at least one variable column besides `p`/`n`",
         call. = FALSE)
  }
  for (v in vars) {
    if (!is.character(df[[v]]) && !is.factor(df[[v]])) {
      stop("target variable `", v, "` must be character or factor",
           call. = FALSE)
    }
  }

  # The p/n weights must be a clean non-negative numeric column; p must sum to 1.
  key <- if (has_p) "p" else "n"
  wts <- df[[key]]
  if (!is.numeric(wts) || anyNA(wts) || any(wts < 0)) {
    stop("`", key, "` must be numeric, non-negative, and free of NAs",
         call. = FALSE)
  }
  if (has_p && abs(sum(wts) - 1) > 1e-6) {
    stop("`p` must sum to 1 (got ", format(sum(wts)), ")", call. = FALSE)
  }
  if (has_n && sum(wts) <= 0) {
    stop("`n` counts must sum to a positive value", call. = FALSE)
  }

  model@targets <- c(model@targets, list(list(df = df, tol = tol)))
  model
}

#' Add a continuous mean calibration
#'
#' @param model A `raking_model`.
#' @param formula A formula `mean(var) ~ value`, captured unevaluated.
#' @param tol Tolerance band half-width (>= 0), in the variable's own units.
#'   `tol = 0` (default) matches the mean exactly; `tol > 0` allows the achieved
#'   mean within +/- tol of the target, pulled toward it by the model `loss`.
#' @return The updated `raking_model`.
#' @examples
#' raking() |> calibrate(mean(income) ~ 58000)
#' raking() |> calibrate(mean(age) ~ 40, tol = 0.5)
#' @export
calibrate <- function(model, formula, tol = 0) {
  if (!inherits(formula, "formula")) {
    stop("`calibrate()` needs a formula like mean(income) ~ 58000",
         call. = FALSE)
  }
  if (tol < 0) stop("`tol` must be >= 0", call. = FALSE)
  model@calibrations <- c(model@calibrations,
                          list(list(formula = formula, tol = tol)))
  model
}

#' Set the weight penalty
#'
#' @param model A `raking_model`.
#' @param reg A penalty from `entropy()`, `kl()`, or `sum_squares()`.
#' @param lambda Penalty multiplier (>= 0). Default 1.
#' @return The updated `raking_model`.
#' @examples
#' raking() |> penalty(entropy())
#' raking() |> penalty(sum_squares(), lambda = 2)
#' @export
penalty <- function(model, reg, lambda = 1) {
  if (!is.null(model@penalty)) message("replacing existing penalty")
  if (lambda < 0) stop("`lambda` must be >= 0", call. = FALSE)
  model@penalty <- reg
  model@lambda <- lambda
  model
}

#' Set element-wise weight bounds
#'
#' @param model A `raking_model`.
#' @param min Lower bound on weights (>= 0). Default 0.
#' @param max Upper bound on weights (> 0). Default Inf.
#' @return The updated `raking_model`.
#' @examples
#' raking() |> bounds(min = 0.5, max = 2)
#' @export
bounds <- function(model, min = 0, max = Inf) {
  if (!is.null(model@bounds)) message("replacing existing bounds")
  model@bounds <- list(min = min, max = max)
  model
}
