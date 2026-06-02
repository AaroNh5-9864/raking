# Project a vector onto the probability simplex {w : w >= 0, sum(w) = z}.
# Sort-based algorithm (Held et al. 1974; same as regrake and Python rsw).
# Internal, not exported.
projection_simplex <- function(v, z = 1) {
  n <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u > (cssv - z) / seq_len(n)))
  theta <- (cssv[rho] - z) / rho
  pmax(v - theta, 0)
}

# Project v onto {w : lower <= w <= upper, sum(w) = z}. Finds the threshold
# theta where sum(clip(v - theta, lower, upper)) = z by bisection over the
# breakpoints, then interpolates. Internal, not exported.
projection_bounded_simplex <- function(v, lower, upper, z = 1) {
  n <- length(v)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  # An infinite upper bound never binds on the simplex (no weight can exceed the
  # total z), so clamp it to z. This keeps the breakpoints finite -- otherwise a
  # `v - Inf = -Inf` breakpoint poisons the interpolation with NaN.
  upper[!is.finite(upper)] <- z

  clip <- function(theta) pmax(pmin(v - theta, upper), lower)

  # If plain clipping already sums to z, that's the projection.
  w <- pmin(pmax(v, lower), upper)
  if (abs(sum(w) - z) < 1e-12) return(w)

  # sum(clip(theta)) is decreasing in theta; bracket the solution.
  breakpoints <- sort(c(v - lower, v - upper))
  lo <- 1L
  hi <- 2L * n
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (sum(clip(breakpoints[mid])) > z) lo <- mid else hi <- mid
  }
  t_lo <- breakpoints[lo]
  t_hi <- breakpoints[hi]
  s_lo <- sum(clip(t_lo))
  s_hi <- sum(clip(t_hi))
  theta <- if (abs(s_lo - s_hi) < 1e-15) {
    t_lo
  } else {
    t_lo + (s_lo - z) * (t_hi - t_lo) / (s_lo - s_hi)
  }
  clip(theta)
}

# ADMM solver for regularized survey weighting.
#   F      : m x n design matrix (rows = constraints, cols = sample units).
#   losses : list of S7 loss objects; their @target blocks stack to F's m rows.
#   reg    : an S7 regularizer object.
#   lam    : regularization strength.
# Returns a list; w_best is the calibrated weight vector.
admm <- function(F, losses, reg, lam,
                 rho = 50, maxiter = 5000,
                 eps_abs = 1e-5, eps_rel = 1e-5,
                 bounds = NULL, verbose = FALSE) {
  F <- Matrix::Matrix(F, sparse = TRUE)
  m <- nrow(F)
  n <- ncol(F)

  # Each loss owns a contiguous block of rows of F. Work out the boundaries.
  sizes <- vapply(losses, function(l) length(l@target), integer(1))
  if (sum(sizes) != m) stop("loss target sizes must sum to nrow(F)")
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L

  # ADMM variables, warm-started at uniform weights.
  f <- as.numeric(Matrix::rowMeans(F))
  w <- rep(1 / n, n); w_bar <- w; w_tilde <- w
  y <- rep(0, m); z <- rep(0, n); u <- rep(0, n)

  # KKT matrix Q = [[2I, t(F)], [F, -I]], factorized ONCE and reused each
  # iteration. Tiny diagonal damping keeps the sparse LDL^T factor stable.
  Q <- rbind(
    cbind(2 * Matrix::Diagonal(n), Matrix::t(F)),
    cbind(F, -Matrix::Diagonal(m))
  )
  damp <- max(1e-8, max(Matrix::rowSums(abs(Q))) * 1e-6)
  Q <- Q + damp * Matrix::Diagonal(nrow(Q))
  Q_factor <- Matrix::Cholesky(Q, LDL = TRUE, perm = TRUE)

  rhs <- numeric(n + m)

  # Resolve the prox/evaluate methods once, up front, so the hot loop calls
  # them directly instead of paying S7 method dispatch on every iteration.
  loss_prox <- lapply(losses, function(l) S7::method(prox, S7::S7_class(l)))
  reg_prox <- S7::method(prox, S7::S7_class(reg))
  # The returned weights are the final feasible projection iterate w_bar,
  # consistent with rsw/regrake (best-iterate tracking in the paper applies to
  # the out-of-scope representative-selection path).

  converged <- FALSE
  r_norm <- s_norm <- NA_real_
  k <- 0L   # so `iterations` is defined even if maxiter < 1 (no loop body)

  for (k in seq_len(maxiter)) {
    Fw <- as.numeric(F %*% w)

    # f-update: every loss applies its own prox to its block. Because the S7
    # method already carries the loss's target/bounds, this call is uniform --
    # no branching on loss type. (This is the payoff of the S7 redesign.)
    for (i in seq_along(losses)) {
      idx <- starts[i]:ends[i]
      f[idx] <- loss_prox[[i]](losses[[i]], Fw[idx] - y[idx], 1 / rho)
    }

    # w_tilde: regularizer prox.  w_bar: projection onto the (bounded) simplex.
    w_tilde <- reg_prox(reg, w - z, lam / rho)
    w_bar <- if (is.null(bounds)) {
      projection_simplex(w - u)
    } else {
      projection_bounded_simplex(w - u, bounds$lower, bounds$upper)
    }

    # Solve the KKT system for the new w using the cached factorization.
    rhs[1:n] <- as.numeric(Matrix::crossprod(F, f + y)) + w_tilde + z + w_bar + u
    w_old <- w
    w <- as.numeric(Matrix::solve(Q_factor, rhs))[1:n]

    # Dual variable updates.
    Fw <- as.numeric(F %*% w)
    y <- y + f - Fw
    z <- z + w_tilde - w
    u <- u + w_bar - w

    # Primal/dual residuals and tolerances (Boyd et al. 2011).
    r_norm <- sqrt(sum((f - Fw)^2) + sum((w_tilde - w)^2) + sum((w_bar - w)^2))
    s_norm <- rho * sqrt(sum((Fw - f)^2) + 2 * sum((w - w_old)^2))
    p <- m + 2 * n
    Ax  <- sqrt(sum(f^2) + sum(w_tilde^2) + sum(w_bar^2))
    Bz  <- sqrt(3 * sum(w^2))
    ATy <- rho * sqrt(sum(y^2) + sum(z^2) + sum(u^2))
    eps_pri  <- sqrt(p) * eps_abs + eps_rel * max(Ax, Bz)
    eps_dual <- sqrt(p) * eps_abs + eps_rel * ATy

    if (verbose && k %% 50 == 0) {
      message(sprintf("it %d | r %.2e/%.2e | s %.2e/%.2e",
                      k, r_norm, eps_pri, s_norm, eps_dual))
    }
    # Converged on primal/dual residuals. A tol>0 margin is a penalized band
    # whose prox already clips the loss copy into the band each iteration, so a
    # converged fit satisfies the band up to the residual tolerance; no separate
    # band check is needed.
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      converged <- TRUE
      break
    }
  }

  w_best <- w_bar

  list(w = w, w_bar = w_bar, w_tilde = w_tilde, w_best = w_best,
       y = y, z = z, u = u, iterations = k, converged = converged,
       primal_residual = r_norm, dual_residual = s_norm)
}
