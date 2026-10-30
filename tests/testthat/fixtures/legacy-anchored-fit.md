# Legacy anchored fit

`legacy-anchored-fit.rds` is the actual native fit saved before the whole-facet
centering change in `tmp/anchor-compatibility-20261014/tam-probe.rds` (`$native`).
It has not been relabeled or re-estimated. It uses the packaged synthetic
`example_core` ratings: 48 Persons, four Raters, four Criteria, scores 1--4.

The RSM MML fit has R01 fixed at 0.7, Person uncentered, legacy remainder
centering, fixed N(0,1), 31 integration points, maxit 400 and reltol 1e-9.
Its stored top-level, replay and facet specifications have no `centering_scope`
metadata. The other raters sum to zero, so all raters together sum to 0.7.
Use it to check preservation of old constraints and reporting, not as a
coverage study or evidence that these settings suit all datasets.

Copied without alteration by:

```r
fit <- readRDS("tmp/anchor-compatibility-20261014/tam-probe.rds")$native
saveRDS(fit, "tests/testthat/fixtures/legacy-anchored-fit.rds", compress = "xz")
```
