# Fixed-normal GPCM native-location fixture

`native-location-fixed-gpcm.rds` is the unchanged native fit of the prespecified
A/B full-panel N=60, SD=1, S-GPCM replicate 1. There are three Rater levels,
three Criterion levels, four categories and 540 ratings. Criterion owns both
steps and the geometric-mean-one relative slopes. This fixture was selected
before fitting; no data were generated or replaced for the check.

The single October 1 development fit uses direct fixed-grid MML, q=31, BFGS,
maxit=400, reltol=1e-9, preserved categories and the explicit
`gpcm_mml_identification="fixed_standard_normal"` restriction. Source/input
hashes, arguments, warnings and the initial interval error are retained under
`validation-results/native-location-branches-20261001/`. The fit legitimately
records `posterior_basis="legacy_mml"`; the location consumer initially rejected
it by expecting a different label. The repair changes the consumer's source
check, not any fitted parameter or readiness flag.

The repaired consumer passes the full-information and same-point q=61 review.
The scaled gradient is approximately 7.815e-6, grid score displacement .005920
and standardized covariance change .007559, within its existing thresholds.
Tests reuse this fit without optimization or simulation, verify covariance
transformation and reject inconsistent population coding. They also check
saved result/report/plot/export replay with numerical recomputation blocked.
This is a branch/regression witness, not a coverage study or a completed
broad-map calibration job.
