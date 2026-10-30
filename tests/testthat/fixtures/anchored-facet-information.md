# Conditional anchored facet information

The three retained fits use all 768 ratings in packaged `example_core`, with
48 Persons, four raters, four criteria and four categories. Criterion is
uncentered, Rater is whole-facet centered, and the population is fixed N(0,1).
RSM fixes R01 at .7; PCM fixes R01 at -.4 and has Criterion-owned steps.
The group RSM fixes R01 at .7 and the mean of R01--R03 at .1. Thus R04 is
also fixed at -.3, while R02 and R03 remain separately estimated.

The ordinary RSM/PCM fits are reused unchanged. Only the missing group case
was fitted for this check. All use q81 and pass their ordinary readiness
checks. The independent reference uses the retained TAM 4.3.25 response-design
arrays and manually constructed affine coefficient maps. For each Person,
the marginal score is the posterior mean of the complete-data score; observed
information is the posterior mean complete-data information minus the
posterior score covariance. The full joint inverse includes criterion and
step nuisance coordinates. Independent finite differences check this identity.

Reference q161 calculations at the same retained vectors check integration
sensitivity. Maximum full-covariance standardized change is 1.686e-5 and
maximum Rater SE change is 3.147e-9. Tests compare the production covariance,
scores, model/sandwich contrasts and fixed targets to these stored references.
They also verify larger clusters, rank refusal and saved output reuse.

The existing execution evidence and recipe are in
`tmp/anchor-alignment-20261014/a4-conditional-information.*`; source identity
and numerical discrepancies are retained in this fixture's provenance and
checks. This is a constraint/implementation check, not an empirical coverage
study, a source-anchor uncertainty calculation or a claim of equivalence to
TAM's printed standard errors.
