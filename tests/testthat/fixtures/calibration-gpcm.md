# Portable GPCM regression fixtures

`calibration-gpcm-shared.rds` and `calibration-gpcm-separate.rds` are test
containers, each with `calibration`, `rows`, and `native` elements. Only
`calibration` is a portable artifact: it contains no native fit or training
responses. `rows` contains four synthetic scoring patterns and `native` retains
the corresponding fitted-object predictions for regression comparison.

Both artifacts were extracted, validated and frozen through the public API on
2026-09-26, with passing fresh conditional source checks. They retain global
boundary/identification states `not_evaluated` and inference readiness FALSE.
Neither artifact is evidence that the source has a finite global maximum.

- Shared Rater slope/step owner: the synthetic fit in
  `validation-results/gpcm-probability-refit-20260925/refit-788.rds` (`fit`).
  Source integration order 61, comparison 121. Artifact scoring order 121.
- Separate Criterion slope/Rater step owners: native fixture
  `mfrm-conditional-scoring-gpcm.rds` (`fit`), documented alongside that file.
  Source integration order 31, comparison 61. Artifact scoring order 31.

Scoring rows comprise the first two observed synthetic persons plus copies of
one person's contexts with all scores at the minimum or maximum. Native scores
use the same scoring order as each artifact. The initial shared-owner attempt
at scoring order 31 failed the per-person numerical check on these patterns;
the refined fixture retains the same calibration. No optimizer was rerun.

`test-portable-gpcm.R` checks the artifact against native results and a scalar
adjacent-logit implementation integrated independently over the real line.
It also tests altered priors, incomplete/malformed records, coarse scoring,
new-process replay, summary/plot/CSV preservation and existing lifecycle rules.
This fixed-calibration test does not measure parameter recovery or coverage.
See `inst/validation/portable-gpcm-20260927.md` for the implementation record.
