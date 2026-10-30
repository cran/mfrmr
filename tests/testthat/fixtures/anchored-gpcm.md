# Anchored one-family GPCM regression fixture

`anchored-gpcm.rds` contains three native MML fits to the existing synthetic
`example_core` ratings (48 Persons, four Raters, four Criteria, scores 1--4).
These are implementation checks, not a coverage study or real rater data.

All fits fix R01 at 0.7, use whole-facet centering, Criterion-owned positive
relative slopes with geometric mean one, BFGS, `maxit = 400`, `reltol = 1e-10`
and `category_policy = "preserve"`.

- `shared`: Criterion owns steps and slopes; default estimated normal mean
  and variance; 61 fixed quadrature points.
- `separate_group`: Rater owns steps, Criterion owns slopes; default estimated
  normal population; mean R01:R03 fixed at 0.1, implying R04 = -0.3.
  A 61-point fit was retained for sensitivity assessment; the fixture uses
  its 121-point refit.
- `fixed_population`: Criterion owns steps and slopes; fixed N(0,1) with
  `noncenter_facet = "Criterion"`; 61 fixed quadrature points.

Each entry retains independent affine coefficient equations, adjacent-category
probabilities, normal-mixture log likelihood, numerical information inverse,
Person score sandwich and curve Jacobian. Reference integration uses statmod;
parameter expansion does not use the package constraint helpers. A derivative
check near a zero coordinate required absolute-scaled differences to avoid
cancellation in the independent numerical reference. Production gradients and
numerical thresholds were unchanged.

The `scoring` element holds one new-Person response pattern and independent
continuous-integration EAP/SD. Shared-owner and fixed-population source fits do
not pass the stricter fitted-object scoring integration checks: default scoring
must refuse them; explicit review scores remain review-only. The refined
separate-owner source passes conditional scoring. All three retain global
inference readiness FALSE; target-specific slope/curve checks do not change it.

Regeneration and original results are retained in the repository under
`tmp/anchor-alignment-20261014/a4-gpcm-probe.R`, `a4-gpcm-refine.R`,
`a4-gpcm-reference-finish.R` and `a4-gpcm-scoring.R`. Failed initial checks
remain in their original logs. This fixture adds no claims about global
optimality, source-anchor uncertainty or repeated-sample interval coverage.
