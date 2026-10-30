# Public adaptive fit from explicit native starts

`gmfrm-explicit-start.rds` retains an actual public `fit_mfrm()` call on the
existing N20 binary 2-by-2 dataset in `gmfrm-result.rds`. The supplied roster
contains the retained segmented local estimate and the original neutral vector.
The adaptive objective is evaluated and optimized independently for both;
neither the segmented likelihood nor its qualification is imported. Both
public optimizer/gradient checks pass. No response data are generated.

The fixture includes exact call arguments, the public fit and all its start
attempts, four component Wald intervals, one Judge-location contrast and
fixed-calibration scoring output. The latter uses the same response records
with different Person IDs; it is a numerical scoring witness, not independent
new-Person performance evidence. CI coverage remains unqualified.

Source/input-bound live fitting, public outputs, portable calibration and
fresh-process replay are recorded under
`validation-results/gmfrm-public-start-20261011/`. Fixture source hashes identify
the manifest and completed phase files. Unit tests replay the stored optimizer
outcomes while checking the current public dispatch, initial adaptive objectives,
metadata, refusals and downstream output contracts. They do not repeat the
optimizer, coverage simulation or full-information calculation.

On 2026-10-15, only the fixed-calibration scores were regenerated from the
retained fit and the same scoring rows to adopt the R-version-independent
source identifier. Estimates agreed within 1e-12; the fit, input arguments,
intervals and original provenance were unchanged. No optimizer was rerun.
