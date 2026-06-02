# raking

An R implementation of Optimal Representative Sample Weighting (Barratt,
Angeris & Boyd, 2021): find survey weights that make a sample match known
population targets, using the ADMM solver shared with the Python `rsw` and R
`regrake` packages. Built on S7 classes with a data-free, pipe-based model API.

## Installation

```r
# install.packages("remotes")
remotes::install_github("<your-username>/raking")
```

## Usage

You build a data-free model by piping steps together, then `fit()` it to data.

```r
library(raking)

model <- raking(name = "wt") |>
  target(data.frame(sex = c("M", "F"), p = c(0.48, 0.52))) |>  # exact margin
  calibrate(mean(income) ~ 58000) |>                            # exact mean
  penalty(entropy()) |>                                         # keep weights uniform
  bounds(min = 0.5, max = 2)                                    # weight bounds

fit <- model |> fit(survey_data)

survey_data$wt <- weights(fit)   # length nrow(data), sums to n
summary(fit)                     # balance table + Kish ESS + convergence
```

Categorical margins (`target`) and continuous means (`calibrate`) are matched
exactly by default; pass `tol > 0` to a margin to allow a tolerance band. The
entropy penalty keeps the weights as uniform as possible subject to the
constraints.

## Status

Research preview. See `RAKING-REVIEW-NOTES.md` (in the project) for the current
design notes and open decisions.

## References

Barratt, Angeris & Boyd (2021), *Optimal Representative Sample Weighting*,
SSRN 3604610.
