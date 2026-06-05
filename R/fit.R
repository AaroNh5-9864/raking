#' @include model.R
NULL

# --- fit(): attach data, assemble F, solve (design doc 4.6, 6.1) -----------

# Build F + losses + per-row scales from a model's targets/calibrations + data.
# Internal; driven by the model specs (not by a formula).
build_model_inputs <- function(model, data) {
  n <- nrow(data)
  if (n == 0) stop("`data` has no rows", call. = FALSE)
  blocks <- list()
  losses <- list()
  scales <- list()
  info <- list()   # per-row metadata: variable, level, tol (for the balance table)

  # Global loss key -> a soft loss object. Used as the inner loss of a
  # penalized band when tol > 0 (it pulls the achieved value toward the
  # target inside the band). Irrelevant when tol == 0 (exact matching).
  make_soft_loss <- function(target) {
    switch(model@loss,
      KL = KLLoss(target = target),
      L2 = LeastSquaresLoss(target = target),
      L1 = L1Loss(target = target),
      # Pearson chi-square == weighted least squares with weight 1/sqrt(target).
      # It divides by the target, so targets must be strictly positive.
      chi_square = {
        if (any(target <= 0)) {
          stop("chi_square loss requires strictly positive targets",
               call. = FALSE)
        }
        LeastSquaresLoss(target = target, diag_weight = 1 / sqrt(target))
      },
      stop("loss '", model@loss, "' is not implemented in fit()",
           call. = FALSE)
    )
  }

  # Categorical margins (one loss block per target() call).
  for (spec in model@targets) {
    df <- spec$df
    vars <- setdiff(names(df), c("p", "n"))
    if (length(vars) < 1) {
      stop("a target() df needs at least one variable column", call. = FALSE)
    }
    for (v in vars) {
      if (is.null(data[[v]])) {
        stop("variable `", v, "` not found in data", call. = FALSE)
      }
      if (anyNA(data[[v]])) {
        stop("target variable `", v, "` has NAs in the data", call. = FALSE)
      }
    }
    target <- if ("p" %in% names(df)) df$p else df$n / sum(df$n)
    # One indicator row per df row: a single level, or a cell (combination of
    # all variable columns) for a joint margin. 1 where the data matches that
    # row on every variable column.
    block <- t(vapply(seq_len(nrow(df)), function(i) {
      hit <- rep(TRUE, n)
      for (v in vars) {
        hit <- hit & (as.character(data[[v]]) == as.character(df[[v]][i]))
      }
      as.numeric(hit)
    }, numeric(n)))

    # Doc 8: every target level must occur in the data, and every data unit must
    # fall in exactly one target cell (no unmatched data levels, no overlap).
    row_counts <- rowSums(block)
    if (any(row_counts == 0)) {
      bad <- which(row_counts == 0)[1]
      lvl <- paste(vapply(vars, function(v) as.character(df[[v]][bad]), ""),
                   collapse = ":")
      stop("target level `", lvl, "` for `", paste(vars, collapse = ":"),
           "` does not occur in the data", call. = FALSE)
    }
    col_counts <- colSums(block)
    if (any(col_counts == 0)) {
      stop("some data rows have a level of `", paste(vars, collapse = ":"),
           "` with no matching target", call. = FALSE)
    }
    if (any(col_counts > 1)) {
      stop("target cells for `", paste(vars, collapse = ":"),
           "` overlap (a data row matches more than one)", call. = FALSE)
    }

    if (spec$tol == 0) {
      # tol == 0: hard exact matching (doc 4.2). A zero-width band is just the
      # target point, so the global loss is irrelevant here.
      loss <- EqualityLoss(target = target)
    } else {
      # tol > 0: penalized band (doc 6.2/6.5). The global loss pulls the
      # achieved value toward the target; the band [target +/- tol] is a hard
      # wall it cannot cross.
      loss <- BandedLoss(
        inner  = make_soft_loss(target),
        target = target,
        lower  = target - spec$tol,
        upper  = target + spec$tol
      )
    }
    var_label <- paste(vars, collapse = ":")
    level_labels <- vapply(seq_len(nrow(df)), function(i) {
      paste(vapply(vars, function(v) as.character(df[[v]][i]), ""),
            collapse = ":")
    }, "")
    blocks[[length(blocks) + 1]] <- block
    losses[[length(losses) + 1]] <- loss
    scales[[length(scales) + 1]] <- rep(1, nrow(df))
    info[[length(info) + 1]] <- data.frame(
      variable = rep(var_label, nrow(df)),
      level = level_labels,
      tol = rep(spec$tol, nrow(df)),
      stringsAsFactors = FALSE
    )
  }

  # Continuous mean calibrations (one normalized row each).
  for (cal in model@calibrations) {
    f <- cal$formula
    ctol <- cal$tol
    lhs <- f[[2]]
    if (!is.call(lhs) || as.character(lhs[[1]]) != "mean") {
      stop("calibrate() left-hand side must be mean(var)", call. = FALSE)
    }
    var <- as.character(lhs[[2]])
    value <- eval(f[[3]], environment(f))
    if (is.null(data[[var]])) {
      stop("variable `", var, "` not found in data", call. = FALSE)
    }
    values <- data[[var]]
    if (!is.numeric(values)) {
      stop("calibrate() requires a numeric variable: `", var, "`", call. = FALSE)
    }
    if (anyNA(values)) {
      stop("calibrate() variable `", var, "` has NAs in the data", call. = FALSE)
    }
    if (value == 0) stop("calibrate() target must be non-zero", call. = FALSE)
    blocks[[length(blocks) + 1]] <- matrix(values / value, nrow = 1)
    # The row is normalized so the target is 1. tol = 0 means match the mean
    # exactly; tol > 0 is a penalized band, with the tol expressed in the
    # variable's own units (so it scales by 1/value on the normalized row).
    if (ctol == 0) {
      losses[[length(losses) + 1]] <- EqualityLoss(target = 1)
    } else {
      tn <- abs(ctol / value)
      losses[[length(losses) + 1]] <- BandedLoss(
        inner = make_soft_loss(1), target = 1, lower = 1 - tn, upper = 1 + tn
      )
    }
    scales[[length(scales) + 1]] <- value
    info[[length(info) + 1]] <- data.frame(
      variable = var, level = "mean", tol = ctol, stringsAsFactors = FALSE
    )
  }

  if (length(blocks) == 0) {
    stop("model has no targets or calibrations to fit", call. = FALSE)
  }
  list(F = do.call(rbind, blocks), losses = losses,
       scales = unlist(scales), info = do.call(rbind, info))
}

# --- The fitted result -----------------------------------------------------

raking_fit <- S7::new_class(
  "raking_fit",
  package = "raking",
  properties = list(
    weights     = S7::class_numeric,                 # length n, sum = n
    balance     = S7::class_data.frame,              # doc 7 balance table
    convergence = S7::class_list,                    # converged/iters/residuals
    diagnostics = S7::class_list,                    # weight + margin diagnostics
    model       = raking_model,                      # the model that produced it
    call        = S7::class_any,                     # captured fit() call
    name        = S7::class_character                # weight column name
  )
)

S7::method(print, raking_fit) <- function(x, ...) {
  conv <- x@convergence
  d <- x@diagnostics
  cat("<raking_fit>\n")
  cat("  loss:           ", x@model@loss, "\n", sep = "")
  cat("  units (n):      ", length(x@weights), "\n", sep = "")
  cat("  converged:      ", conv$converged, " (", conv$iterations,
      " iters)\n", sep = "")
  cat("  weight range:   ",
      paste(round(d$weight_range, 3), collapse = " to "), "\n", sep = "")
  cat("  Kish ESS:       ", round(d$kish_ess, 2), "\n", sep = "")
  invisible(x)
}

# --- The fit() entry point -------------------------------------------------

# Resolve base/design weights to a probability distribution b (sums to 1).
resolve_base <- function(base, data, n) {
  if (is.null(base)) return(rep(1 / n, n))
  if (is.character(base) && length(base) == 1) {
    if (is.null(data[[base]])) {
      stop("base column `", base, "` not found in data", call. = FALSE)
    }
    base <- data[[base]]
  }
  if (!is.numeric(base) || length(base) != n) {
    stop("`base` must be a numeric vector of length nrow(data), or a column name",
         call. = FALSE)
  }
  if (any(base <= 0)) {
    stop("`base` weights must be strictly positive", call. = FALSE)
  }
  base / sum(base)
}

# Resolve the effective regularizer and lambda for the solver, given the
# model's penalty and the base distribution b (sums to 1). Returns a
# list(reg, lambda). Decisions (doc 4.4, 6.3):
#   * No penalty() step + uniform base  -> unpenalized (ZeroRegularizer, lam 0).
#   * No penalty() step + design base   -> pull toward the base (KL), lam 1, so
#                                          a supplied base is never silently
#                                          ignored.
#   * entropy() + design base           -> KL toward the base distribution.
#   * sum_squares()                     -> shrink toward the base weights.
#   * explicit kl()/other              -> kept as-is.
resolve_penalty <- function(penalty, lambda, b) {
  n <- length(b)
  # lambda == 0 means unpenalized: use the identity (Zero) regularizer. (An
  # entropy/KL prox with step lambda/rho = 0 would divide by zero.)
  if (!is.null(lambda) && lambda == 0) {
    return(list(reg = ZeroRegularizer(), lambda = 0))
  }
  uniform <- isTRUE(all.equal(b, rep(1 / n, n)))
  if (is.null(penalty)) {
    if (uniform) return(list(reg = ZeroRegularizer(), lambda = 0))
    return(list(reg = KLRegularizer(prior = b), lambda = 1))
  }
  if (S7::S7_inherits(penalty, EntropyRegularizer) && !uniform) {
    return(list(reg = KLRegularizer(prior = b, limit = penalty@limit),
                lambda = lambda))
  }
  if (S7::S7_inherits(penalty, SumSquaresRegularizer)) {
    return(list(reg = SumSquaresRegularizer(base = b), lambda = lambda))
  }
  list(reg = penalty, lambda = lambda)
}

# Re-export the fit() generic from the generics package (doc 4.6), so users can
# call fit() after library(raking) and it dispatches on raking_model.
#' @importFrom generics fit
#' @export
generics::fit

# fit() method for raking_model: attach data, assemble F from the model's
# targets and calibrations, run the solver, return a raking_fit. Arguments:
#   object  - a raking_model
#   data    - a data.frame to weight
#   base    - base/design weights (vector or column name); default uniform
#   control - a rake_control() (or list); extra settings via ...
S7::method(fit, raking_model) <- function(object, data, base = NULL,
                                          control = rake_control(),
                                          ...) {
  cl <- tryCatch(match.call(), error = function(e) NULL)
  inp <- build_model_inputs(object, data)
  n <- nrow(data)
  m <- nrow(inp$F)
  b <- resolve_base(base, data, n)
  pen <- resolve_penalty(object@penalty, object@lambda, b)
  # model bounds are on the sum=n scale; the solver works on sum=1, so /n.
  bnds <- if (is.null(object@bounds)) {
    NULL
  } else {
    list(lower = object@bounds$min / n, upper = object@bounds$max / n)
  }
  # Merge any settings passed via ...; explicit raw eps overrides margin_tol.
  dots <- list(...)
  control <- modifyList(as.list(control), dots)
  if (("eps_abs" %in% names(dots) || "eps_rel" %in% names(dots)) &&
      !("margin_tol" %in% names(dots))) {
    control$margin_tol <- NULL
  }
  ctrl <- resolve_control(control, m, n)
  sol <- admm(inp$F, inp$losses, pen$reg, pen$lambda, bounds = bnds,
              rho = ctrl$rho, maxiter = ctrl$maxiter,
              eps_abs = ctrl$eps_abs, eps_rel = ctrl$eps_rel,
              verbose = ctrl$verbose)

  if (!sol$converged) {
    warning(sprintf(
      "raking did not converge in %d iterations; returning best iterate (convergence$converged = FALSE)",
      sol$iterations), call. = FALSE)
  }

  w <- sol$w_best * n                                  # sum = n scale
  achieved <- as.numeric(inp$F %*% sol$w_best) * inp$scales
  target   <- unlist(lapply(inp$losses, function(l) l@target)) * inp$scales
  difference <- achieved - target
  # "satisfied": within the margin's tol band (plus a small numerical slack, so
  # a converged band sitting at its edge still reads as satisfied) or, for exact
  # rows, within a small tolerance scaled to the target magnitude.
  abstol <- inp$info$tol + 1e-4 * pmax(1, abs(target))
  balance <- data.frame(
    variable   = inp$info$variable,
    level      = inp$info$level,
    target     = target,
    achieved   = achieved,
    difference = difference,
    satisfied  = abs(difference) <= abstol + 1e-12,
    stringsAsFactors = FALSE
  )

  wmean <- mean(w)
  wvar  <- stats::var(w)
  deff  <- 1 + wvar / wmean^2
  denom <- abs(target)
  pct   <- ifelse(denom > 0, abs(difference) / denom, NA_real_)
  diagnostics <- list(
    weight_range = range(w),
    weight_mean  = wmean,
    weight_sd    = stats::sd(w),
    kish_deff    = deff,
    kish_ess     = n / deff,
    max_abs_diff = max(abs(difference)),
    max_pct_diff = 100 * max(pct, na.rm = TRUE)
  )
  convergence <- list(
    converged       = sol$converged,
    iterations      = sol$iterations,
    primal_residual = sol$primal_residual,
    dual_residual   = sol$dual_residual
  )

  raking_fit(
    weights     = w,
    balance     = balance,
    convergence = convergence,
    diagnostics = diagnostics,
    model       = object,
    call        = cl,
    name        = object@name
  )
}
