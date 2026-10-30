# Two-family slope interval numerical reference

`gmfrm-joint-information.rds` retains the 240-Person, three-task, three-rater,
three-category example from `inst/validation/gmfrm-mml-em-20260927.R` (seed
20260927). The generator uses the literal adjacent-category product-slope
equation and fixed standard normal ability distribution. This is synthetic data.

The fixture contains data, the retained 31-node EM parameter vector and the
observed marginal information obtained by independently differentiating the
analytic marginal score in that validation. Parameters are two centered task
locations, three free rater locations, three centered rater-step coordinates,
two centered task-log-slope coordinates and three free rater log slopes.

The interval test uses the current native adapter at the retained optimum and
an independently specified 6-by-13 component Jacobian. It compares the covariance
and log-Wald limits against the inverse of the stored full information, including
cross-family covariance. The numerical qualification also reevaluates the
current likelihood and 31-versus-61-node information and scores without an EM
iteration in the reference calculation. Separate public-entry tests refit this
same dataset at five and 31 nodes to distinguish numerical convergence from
integration stability. No new sample is generated. This is a calculation and output regression
reference, not a sampling-coverage or global-identification study.
