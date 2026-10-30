# FACETS group-origin regression reference

`facets-group-anchor.csv` contains output measures and Rasch-Andrich steps from
four FACETS 4.5.0 64-bit executions on 2026-10-24. Each completed with exit 0.
Executable SHA256: `dfb0afb0faa18f026d1b3b4175f22e42cc3764430eb83cbd368c7a572b3593a1`.

Responses are the retained RSM/PCM rows of `tam-jml-0.2.4-response.csv`
(SHA256 `2c7a39e77c4fe38e25030210d679e3a0e3b978a1b427756449a0da365886cbde`),
324 rows in each model. No responses were simulated or changed. Controls:
Person positive, Person noncentered, Task centered, RSM shared R3 scale or
PCM R3 scale per Task; `Convergence=.0000001,.0000001,0,0`, `Iterations=500`,
`Udecimals=6`, `TotalScore=Yes`, `Xtreme=.3,.5`.

Rater uses `G`, with contributions .2 and 1 for levels 1 and 2 in group 1.
Their sum is 1.2, hence the mfrmr target mean is .6. `origin` leaves Rater 3
unconstrained. `restricted` puts Rater 3 alone in group 2 with value -.2.
These are illustrative constraints, not empirically validated calibrations.

The origin cases agree within 1e-4. The restricted cases are deliberately
**not** equality targets: both native group means hold, but the likelihood
derivative for increasing Rater 1 and decreasing Rater 2 is .215081 (RSM) or
.168563 (PCM). A feasible .001 shift increases log likelihood. Tests retain
the independently evaluated constrained likelihood and its stationary mfrmr
solution, rather than requiring the native discrepancy to be reproduced.
No general SE, coverage, extreme-group or iterative-update equivalence follows.

Specification references: [group elements](https://www.winsteps.com/facetman64/elements.htm)
and [centering](https://www.winsteps.com/facetman64/centerfacet.htm).

The unchanged Rater Scorefiles from the two RSM executions are also retained
as `facets-group-origin-score.txt` (SHA256
`ebe5d9b5d58973c607e6b4857150871d5bc749bc5dcb761f8e90d807e87377f2`) and
`facets-group-restricted-score.txt` (SHA256
`ca578b99598cd4ea14245b41c67e929897c1ec559dde63f4316cd639dad4f57e`).
They preserve fixed-width fields and raw status/group codes for import tests.
Levels 1 and 2 are estimated group members, yet both have Status=1, as does
the singleton in the restricted case. Group codes can also label reporting
groups under direct anchoring; importing Status/Group does not recover the
original Labels= mode. Tests require that ambiguity to remain explicit,
without changing points, reported SEs or numerical availability.

## Extreme group members: import contract

The four `facets-group-extreme-{low,high}-{g,x}-score.txt` files are unchanged
Person Scorefiles from FACETS 4.5.0 executions on 2026-10-24, using the same
executable identity above. All returned exit 0. The retained RSM responses
have P03's six responses set to 0 (`low`) or 3 (`high`); the other responses
are unchanged. P24 has score 8 in seven observations. Person uses G or X,
with P03 and P24 each contributing -0.4 to group 1. Rater is noncentered,
Task centered. Other numerical controls are as above. This is a reporting
probe, not evidence of finite ordinary JML or calibrated uncertainty.

| Boundary/mode | Scorefile SHA256 |
| --- | --- |
| low G | `b64854633bea4ca2e77a92c0c7512c8f862079c4262d0a88fa0eed42f5640d6e` |
| low X | `442a730f15c2e0e426086527bed11f4116c3aba3949bae5ad8369267806cec9b` |
| high G | `d2d018df5646476cabefef24c404394f0629409f576dd915ec42808aa7454e7b` |
| high X | `e042ddce182e46ad0c6900148f8840a7c9901651549eeeec6ad5c8d9d6f1f7eb` |

P03 retains Status -2/-3, Group 1, and numeric Infit/Outfit placeholders
1.00; Table 7 explicitly prints Minimum/Maximum instead of those fit values.
P24 retains Status 1, Group 1 and interior fit statistics. Under G, the two
reported measures average -0.4; under X, both report -0.4. The reader must
preserve these source values without asserting a common tail-adjustment
algorithm or replacing the extreme placeholders with interpretable fit.
These four cases do not establish the behavior of all group configurations,
models or FACETS versions.
Table 7 marks the interior P24 with G under both input modes; it is not a
substitute for the original Labels= specification when distinguishing G/X.
