# Private segmented GMFRM result fixtures

`gmfrm-result.rds` contains the complete N20 binary 2-by-2 run from
`validation-results/gmfrm-run-20261011/live.rds`, together with that file's
digest and original manifest binding. The data have predetermined response
pattern counts, not new random draws. Both prescribed searches returned code 1;
their separate native finishes and independent references qualify locally.
The test never refits the data or recomputes a Hessian.

`gmfrm-result-held.rds` contains the unselected N20-b09 free result from the
retained procedure replay, with the original parent specification and native
start mapping reconstructed in `inst/validation/gmfrm-result-20261011.R`.
Comparison evidence receives the same data/model specification; its likelihoods,
constraints, boundary outcomes and selection decision remain unchanged. It
tests that a boundary reference in the attempt record does not become native
information in the candidate or permit Person scoring.

The current source-bound qualification is under
`validation-results/gmfrm-result-20261011/qualified/`. These are regression
fixtures for result handoff, fixed-calibration scoring and refusals. They are
not independent replications or evidence of sampling coverage.
