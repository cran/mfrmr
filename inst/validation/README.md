# mfrmr validation artifacts

This directory contains non-exported release-review helpers and evidence
artifacts. They are included with the package so that release decisions can be
reconstructed from source files, check logs, and documented validation
criteria.

These files are not user-facing analysis functions. They support release
review, CRAN submission preparation, and future maintenance. Public release
notes stay in `NEWS.md`; implementation-level evidence, source-grounding, and
long-run validation details belong here.

## Primary files

- `release-readiness.R`: release-readiness review. Source this file and run
  `mfrmr_release_readiness_review(pkg_dir = ".")` from the package root. The
  review checks version labels, the local check log, the CI check workflow,
  public terminology, and the release evidence files.
- `release-evidence-map-0.2.0.md`: narrative review map linking release
  claims to mathematical, statistical, UX, documentation, and engineering
  evidence.
- `release-evidence-map-0.2.1.md`: source-grounded evidence map for the
  0.2.1 bounded-`GPCM` recovery-review refinements, including the boundary
  between cited model literature and package-specific validation labels.
- `release-evidence-checklist-0.2.1.csv`: structured checklist used by the
  readiness helper and by manual release review for the current release. The
  0.2.0 checklist is retained as historical release evidence.
- `gpcm-post-0.2.1-roadmap.md`: maintenance roadmap for bounded-`GPCM`
  surfaces that remain caveated, `blocked`, or `deferred` after 0.2.1,
  including score-side review, report/QC bundles, design and screening
  operating characteristics, linking synthesis, posterior predictive checks,
  and heavy-backend extensions.
- `external-parameter-recovery-simulation-0.2.0.md`: compact review of the
  separate common-data parameter-recovery simulation workflow. The large
  generated datasets and engine outputs are not bundled with the package; this
  file records the release-relevant evidence and its limits.
- `external-recovery-audit.R`: optional audit helper that reads a local
  `Parameter_Recovery_Simulation/` output directory, checks expected CSV
  schemas, records file fingerprints, and regenerates the compact evidence
  summary tables used for release review.

## Recommended local sequence

Run these commands from the package root after any source, roxygen, vignette, or
compiled-code change:

```sh
R CMD build .
R CMD check --no-manual --as-cran mfrmr_0.2.1.tar.gz
```

Then run:

```r
source("inst/validation/release-readiness.R")
readiness <- mfrmr_release_readiness_review(pkg_dir = ".")
summary(readiness)
```

The release candidate should have `Status: OK` in the local check log and no
`concern` rows in `readiness$gate_summary`. The check log must also report the
same package version as `DESCRIPTION`; stale logs from an earlier release are
reported as a package-check concern. If the local environment cannot verify
external clock time, record that environment-only NOTE in `cran-comments.md`
and rerun the package check with the clock check disabled to confirm that
package checks are otherwise clean.

CRAN-time tests are intentionally lightweight because CRAN check hosts have
strict timing constraints. Run the full non-CRAN regression surface separately
when release evidence is needed:

```sh
NOT_CRAN=true Rscript -e 'testthat::test_local(".")'
```

If the external common-data simulation workflow has been refreshed, audit it
from the package side before updating the evidence summary:

```r
source("inst/validation/external-recovery-audit.R")
external_review <- mfrmr_review_external_recovery_simulation(
  "../Parameter_Recovery_Simulation"
)
summary(external_review)

source("inst/validation/release-readiness.R")
readiness <- mfrmr_release_readiness_review(
  pkg_dir = ".",
  external_recovery_dir = "../Parameter_Recovery_Simulation"
)
summary(readiness)$external_recovery_status
```

## Cross-platform evidence

GitHub Actions runs the package on macOS, Windows, and Linux across release,
oldrel, and devel R. Warnings are treated as check failures. The workflow also
uploads the check directory as an artifact for each matrix job so that release
review can compare local and CI evidence instead of relying only on the final
job status.

The readiness helper checks the workflow contract from source. It does not
replace reading the uploaded CI artifacts before release submission. The
external parameter-recovery summary is an additional source-grounded review
artifact, not a substitute for rerunning the package tests or the optional
long-running validation scripts.
