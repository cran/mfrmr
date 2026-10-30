# Saved one-slope profile regression fixture

`gpcm-profile-separate.rds` is a serializable `mfrm_slope_intervals` result for
Criterion C1 from the native synthetic fit in
`mfrm-conditional-scoring-gpcm.rds`. It is a saved inference result, not a
portable calibration. Source responses/parameter metadata are retained by the
existing inference identity contract.

The 2026-09-27 development profile used the native order-31 likelihood and
order-61 reference, two starts, maxit 400, two BFGS tolerances and a conditional
single curvature polish. The original lower-root search failed at nuisance
gradient 1.56e-4 despite matching likelihoods. The completed validation reused
25 unchanged passing attempts and computed 13 failed/new attempts; the original
attempts remain under `validation-results/gpcm-profile-intervals-20260927/`.
The fixture's `validation_history` preserves this history. Its elapsed time is
the resumed call time, not a from-scratch performance benchmark.

Both 95% endpoints (approximately 0.6066673 and 0.9021057) passed an independent
continuous integral at the original numerical tolerances. Tests reuse the
result to check plots, tables, reports and export without repeating the whole
nuisance search. The quadratic constraint/search tests exercise the algorithm
independently; the repository runner records actual native-model endpoint
validation. This fixture does not establish repeated-sampling coverage.

See `inst/validation/gpcm-profile-intervals-20260927.md` for all results and limits.
