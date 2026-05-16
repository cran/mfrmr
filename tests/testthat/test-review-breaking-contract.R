test_that("review helpers no longer expose audit compatibility classes", {
  toy <- load_mfrmr_data("example_core")

  anchor_review <- review_mfrm_anchors(toy, "Person", c("Rater", "Criterion"), "Score")
  expect_s3_class(anchor_review, "mfrm_anchor_review")
  expect_false(inherits(anchor_review, "mfrm_anchor_audit"))

  fit <- suppressWarnings(
    fit_mfrm(toy, "Person", c("Rater", "Criterion"), "Score",
             method = "JML", maxit = 10)
  )
  facet_review <- facet_small_sample_review(fit)
  expect_s3_class(facet_review, "mfrm_facet_sample_review")
  expect_false(inherits(facet_review, "mfrm_facet_sample_audit"))
  expect_false("anchor_audit" %in% names(fit$config))
})

test_that("old audit spellings are not exported", {
  old_names <- c(
    "audit_mfrm_anchors",
    "precision_audit_report",
    "audit_conquest_overlap",
    "facet_small_sample_audit",
    "build_weighting_audit",
    "reference_case_audit"
  )
  expect_false(any(old_names %in% getNamespaceExports("mfrmr")))
})
