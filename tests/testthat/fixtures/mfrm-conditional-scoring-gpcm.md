# Conditional scoring regression fixture

Synthetic separate-owner MML GPCM, 120 persons and 720 ratings: Rater owns
steps, Criterion owns slopes, with an estimated intercept-only normal prior.
The retained synthetic input is from the 2026-09-26 joint-profile investigation
(`validation-results/portable-gpcm-development-probe.rds`, repository only).
It is not a portable calibration artifact. The fixture stores data and a native
fit; tests freshly evaluate likelihood, gradient, information and integration.

Regenerated on 2026-09-26 with:

```r
fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score", model = "GPCM",
  step_facet = "Rater", slope_facet = "Criterion", mml_integration = "adaptive",
  quad_points = 31, maxit = 400, reltol = 1e-10, optimizer = "BFGS")
```

Reusing this fit avoids repeatedly optimizing an unchanged positive-control
case in tests. The fixture retains global review status: local conditional
score eligibility must not rewrite global inference/boundary flags.
