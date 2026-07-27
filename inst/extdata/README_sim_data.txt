Synthetic example datasets
==========================

All bundled datasets are synthetic long-format rating data. The compact
examples and the Eckes and Jin (2021)-inspired designs are available as
lazy-loaded R objects and through `load_mfrmr_data()`.

Access from R:

    # Canonical loader (returns a long-format data.frame):
    mfrmr::load_mfrmr_data("example_operational")
    mfrmr::load_mfrmr_data("example_core")
    mfrmr::load_mfrmr_data("example_bias")
    mfrmr::load_mfrmr_data("study1")
    mfrmr::load_mfrmr_data("study2")
    mfrmr::load_mfrmr_data("combined")
    mfrmr::load_mfrmr_data("study1_itercal")
    mfrmr::load_mfrmr_data("study2_itercal")
    mfrmr::load_mfrmr_data("combined_itercal")

    # Compare intended teaching roles before selecting a dataset:
    mfrmr::list_mfrmr_data(details = TRUE)

    # Equivalent base-R form (the lazy-loaded objects are named ej2021_*):
    data("ej2021_study1", package = "mfrmr")

    # Export one of the datasets to a CSV file for external tooling:
    write.csv(mfrmr::load_mfrmr_data("study1"),
              "ej2021_study1.csv",
              row.names = FALSE)

Source
------

Synthetic example datasets:
- `example_operational`: 48 persons, 6 raters, 3 criteria, connected
  two-rater assignment, and six planned omissions. Scores are sampled from
  explicit RSM category probabilities with a fixed seed; omissions are absent
  long-format rows rather than NA or sentinel scores.
- `example_core`: idealized complete crossing for fast examples.
- `example_bias`: balanced partial assignment with deliberately planted
  group-by-criterion and rater-by-criterion effects, plus group latent means of
  -0.1 and +0.1 logits documented in `?mfrmr_example_data`.

Synthetic designs inspired by Eckes & Jin (2021):
- Study 1 design target: 307 examinees, 18 raters, 3 criteria, 4-category scale
- Study 2 design target: 206 examinees, 12 raters, 9 criteria, 4-category scale

Reference
---------

Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in
writing assessment: A Bayesian facets modeling approach.
Psychological Test and Assessment Modeling, 63(1), 65-94.

Notes
-----

- The exact response-generation code and random seed for the legacy
  Eckes-and-Jin-inspired objects are not available; do not use them as
  parameter-recovery evidence.
- The `*_itercal` objects are legacy synthetic variants, not empirical
  calibration standards or external validation evidence.
- `combined` objects reuse person and rater labels across studies. Relabeling
  prevents accidental collisions but creates disconnected components; it does
  not link the study scales. Analyze studies separately unless an explicit
  identity and anchor/linking design establishes the common scale.
- These are synthetic datasets (not original TestDaF operational
  records).
