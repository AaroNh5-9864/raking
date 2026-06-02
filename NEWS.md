# raking 0.0.0.9000

First development version. An R-native implementation of Optimal Representative
Sample Weighting (Barratt, Angeris & Boyd, 2021), built on S7 with a data-free,
pipe-based model API and an ADMM solver that reproduces the `rsw` / `regrake`
results.

## Features

* `raking()` builds a data-free weighting model; `target()`, `calibrate()`,
  `penalty()`, and `bounds()` extend it; `fit()` attaches data and solves.
* Categorical margins and continuous mean calibrations, each either exact
  (`tol = 0`) or matched within a penalized tolerance band (`tol > 0`): the loss
  pulls the achieved value toward the target while a hard wall keeps it within
  `+/- tol`.
* Losses `"KL"`, `"L2"`, `"L1"`, `"chi_square"` (used inside bands); penalties
  `entropy()`, `kl()`, `sum_squares()`; element-wise weight bounds enforced by
  bounded-simplex projection.
* Base/design weights via `fit(base = ...)`.
* `fit()` returns a `raking_fit` with `weights`, a `balance` table,
  `convergence`, `diagnostics` (including Kish deff/ESS), the `model`, and the
  `call`. `weights()`, `summary()`, and `print()` methods are provided.

## Notes

* Returns the final feasible iterate, consistent with `rsw` / `regrake`.
* The bounded sparse solve has a ~2e-4 accuracy floor for soft (L2) problems;
  exact constraints match to the convergence tolerance.
* `calibrate(var())` and `calibrate(quantile())` are planned, not yet
  implemented.
