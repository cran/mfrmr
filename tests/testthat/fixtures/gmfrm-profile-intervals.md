# Two-family component profile reference

`gmfrm-profile-intervals.rds` reuses the synthetic 240-Person, three-task,
three-rater dataset documented in `gmfrm-joint-information.md`. It adds no
sample or coverage evidence. The fixture contains the fixed-N(0,1) fit and
actual public `confint()` results for Task t3 and Rater r3: the first target
is the dependent geometric-mean-one component; the second is a free component.
An earlier 61-node Rater profile with a genuinely unresolved upper endpoint is
also retained for output checks, with its original source identity.

The saved 31-node parameter vector was refined at 61 then 121 nodes using
`mfrm_gmfrm_em(..., maxit = 500, score_tol = 1e-7)` and assembled with the
native `mfrm_gmfrm_fit_result()` adapter. This is an explicit same-data
numerical refinement, not an automatic profiling fallback or a claim about
neutral-start recovery. The final public calls are:

```r
task <- confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3"))
rater <- confint(fit, parm = "slopes", method = "profile", slope = c(Rater = "r3"))
```

Both use unchanged profile controls: two starts, 400 BFGS iterations per stage,
up to eight bracketing steps, an initial log-slope step of .5, and 95% pointwise
chi-square(1) cutoffs. The 121-node objective is checked against 241 nodes at
every constrained solution. Complete optimizer stages and full free parameter
vectors remain in the profile attributes; hashes identify the computing source.

Task t3 gives [0.898828, 1.227248], compared with Wald [0.900647, 1.228492].
Rater r3 gives [1.069859, 1.585058], compared with Wald [1.074609, 1.590472].
The two searches took 133.683 and 127.729 seconds locally. These small
differences in one example do not establish a coverage advantage.

The earlier 61-node profiles each had an unresolved endpoint and remain under
`validation-results/gmfrm-profile-20260929/`. The Task search exposed a missed
curvature-refinement trigger: the fitting-grid nuisance gradient was just
inside 1e-4 while the higher-grid gradient was just outside. The shared helper
now attempts refinement for either score failure; no acceptance threshold
changed. The Rater search instead had a substantive quadrature discrepancy
at an outward trial. The subsequent explicit grid refinement returned both
bounds; that trial's failure alone does not establish that a finer grid is
necessary at the actual endpoint. Neither cause establishes mathematical
unboundedness.

The tests independently reconstruct the literal response equation and marginal
likelihood at all four endpoints, retain owner/level identity including renamed
overlapping labels, and exercise plots/reports/export/replay without rerunning
these searches. They also retain the low-grid failure rules. Sampling coverage,
weak/boundary cases and independent confirmation remain open.

## Retained numerical failures

`gmfrm-numerical-repair.rds` contains two original common-Person SD=.5
datasets from the fixed 2026-09-29 feasibility study (240 Persons, three tasks,
six raters, scores 0:2). `em` is replicate 2, which missed the declared 1e-7
score tolerance after 500 iterations. `profile` is replicate 4;
`log_slope` is its first failed Task t3 trial and `old_checks` retains both
constrained-start outcomes, including their NLL disagreement greater than 80.
`source_hash` identifies the original computing source. The observed data
are inside each fit's `gmfrm$specification`; these are development regressions,
not new independent sampling evidence.

The tests restart the first dataset from neutral parameters at the unchanged
121-node grid, 500-iteration ceiling and 1e-7 score tolerance. They verify
fixed-count M-step refinement, marginal/Q ascent, actual score convergence,
and preservation of failed refinement attempts. On the second dataset,
retained and neutral constrained starts must agree under the optimizer's
per-Person scale, while all likelihood-ratio and acceptance calculations
remain in total-likelihood units. Both tests compare the results with a
literal 241-node response-formula likelihood independent of the fitted kernel.
Analytic profile tests separately check contraction after an inaccurate
outward trial, genuine failures before an endpoint, and failed interior points.

`recovered` contains the revised public fit and Task t3 profile for rotating-pair
SD=.5 replicate 16 from the matched 2026-09-30 replay, plus its computing-source
hashes. Both endpoints are computed (0.9434804, 1.8242810), while two inaccurate
outward trials remain in the profile and attempt tables. Regression tests
verify that plots, reports and saved objects retain those failures without
reprofiling or recalculating covariance. The full public export/replay and
rendered plot were also inspected during validation. This is evidence about
numerical recovery and output preservation, not sampling coverage.

`weak_profile` contains the revised EM fit for common-Person SD=.5 replicate
11 and the original 400-iteration constrained attempts at its first lower
outward trial. Both attempts reach the optimizer limit as a nuisance Rater r2
slope becomes small. The regression retries from the original fit and neutral
parameters at the unchanged 400-iteration stage limit. It requires a recorded
curvature-coordinate restart, total-likelihood/gradient agreement between
starts, the exact component constraint, and the independent 241-node literal
likelihood. A deliberately inaccurate reference integral must still fail;
unavailable curvature scaling must retain the failed result and reason.
