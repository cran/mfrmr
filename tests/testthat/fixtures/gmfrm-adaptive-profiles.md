# Adaptive two-family component profile fixture

This fixture reuses the N20 fit in `gmfrm-rank-deficient-source.rds` and the
Task t3 profile saved by `inst/validation/gmfrm-adaptive-profile-20261010.R`.
It contains no newly generated data. Task t3 is the dependent first-family
slope: its log slope is minus the sum of the other two free log slopes.

Both endpoints remain unavailable because a constrained evaluation fails
the unchanged numerical checks. Failed candidates, both starts and their
optimizer records remain in the fixture. The fixture is evidence for routing,
failure reporting and replay, not for coverage or an unbounded interval.

`inst/validation/gmfrm-adaptive-profile-reference-20261010.R` independently
reconstructs the response formula, sum-zero constraints, posterior modes,
integration scales and Louis information. At the retained source, continuous
integration and independently constructed adaptive quadrature agree within
3.56e-13 in total negative log likelihood. The fixture retains that reference
and input/profile hashes. The test compares the native moving-node objective
and gradient with this reference, then uses finite differences at a changed
parameter vector to check that integration nodes are recalculated.

The test reuses the saved search instead of repeating minutes of optimization.
It also checks failed-endpoint tables, English plot labels, saved reports and
export/reload with refitting and re-profiling disabled. Real profile searches
for both slope families at N20 and N240, and independent checks of every
computed endpoint, are retained separately in
`validation-results/gmfrm-adaptive-profile-20261010/`. These four numerical
cases are not a sampling-performance study.
