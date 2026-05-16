Packaged simulation datasets
============================

Starting in mfrmr 0.1.6, the Eckes and Jin (2021)-inspired simulation
datasets are shipped as lazy-loaded R objects under `data/` rather
than as duplicate CSV files. This keeps the installed package size
below 5 MB without removing any data.

Access from R:

    # Canonical loader (returns a long-format data.frame):
    mfrmr::load_mfrmr_data("study1")
    mfrmr::load_mfrmr_data("study2")
    mfrmr::load_mfrmr_data("combined")
    mfrmr::load_mfrmr_data("study1_itercal")
    mfrmr::load_mfrmr_data("study2_itercal")
    mfrmr::load_mfrmr_data("combined_itercal")

    # Equivalent base-R form (the lazy-loaded objects are named ej2021_*):
    data("ej2021_study1", package = "mfrmr")

    # Export one of the datasets to a CSV file for external tooling:
    write.csv(mfrmr::load_mfrmr_data("study1"),
              "ej2021_study1.csv",
              row.names = FALSE)

Source
------

Simulated for package validation, inspired by Eckes & Jin (2021):
- Study 1 design target: 307 examinees, 18 raters, 3 criteria, 4-category scale
- Study 2 design target: 206 examinees, 12 raters, 9 criteria, 4-category scale

Reference
---------

Eckes, T., & Jin, K.-Y. (2021). Measuring rater centrality effects in
writing assessment: A Bayesian facets modeling approach.
Psychological Test and Assessment Modeling, 63(1), 65-94.

Notes
-----

- The `*_itercal` variants are iterative-calibrated synthetic datasets
  tuned to match target score distribution / fit / reliability
  profiles; the other variants are the baseline synthetic datasets.
- These are synthetic datasets (not original TestDaF operational
  records).
