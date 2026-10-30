# Scoring-prior regression reference

This small extract reuses the completed fixed-calibration sensitivity study:
`inst/validation/gpcm-prior-sensitivity-20260926.R` and its adjacent result record.
It uses `mfrm-conditional-scoring-gpcm.rds`, not a newly fitted model. Persons
P1 (ordinary pattern) and P34 (all lowest scores) are retained at the declared
nested one/three/six-rating subsets, under retained, higher-mean and wider priors.

The CSV contains unrounded production EAP, posterior SD and interval endpoints
that were independently checked using scalar continuous integration. Across
the full study, maximum EAP/SD errors were 8.91e-8/1.55e-7 and endpoint CDF error
1.86e-12. Regression equality at a tighter tolerance checks API routing against
those saved calculations; it does not tighten the oracle's error guarantee.
The source calibration is fixed and this is not a coverage experiment.
