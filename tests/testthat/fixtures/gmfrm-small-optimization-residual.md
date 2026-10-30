# Small local optimization residual

The RDS retains replicate 1 of the common-Person, ability-SD .5 condition in
`inst/validation/gmfrm-sparse-intervals-20260928.R` (seed 92810001): synthetic
data and the 31-point EM parameter vector at mean-score tolerance 1e-6.
Its curvature-scaled gradient is about 0.000190; the original 1e-4 interval
refusal threshold excluded it despite that small displacement in SE units.

A neutral-start 61-point fit at mean-score tolerance 1e-7 supplies the saved
reference parameters and intervals. Maximum component log-slope movement is
under 0.0001 of its original local SE; interval limits differ by less than
0.001 SE. This verifies the numerical meaning of the warning on a real fit,
not coverage. The regression reconstructs the original native object at its
retained optimum without refitting and verifies that the warning, estimates
and numerical checks reach print/plot/report output. The ordinary rank,
information, mean-score and quadrature requirements remain in force.

The `unstable` member retains rotating-pair SD .5 replicate 26 (seed 92810026).
Its mean score is below 1e-6, but its standardized displacement is about .0675,
with one rater slope near .00661. It verifies that a materially larger local
correction still withholds intervals despite a small raw score. It is not
used to claim that this finite point is a globally identified maximum.
