test_that("compatibility_alias_table exposes retained aliases and preferred names", {
  alias_tbl <- compatibility_alias_table()

  expect_true(is.data.frame(alias_tbl))
  expect_true(all(c("Alias", "PreferredName", "Surface", "Lifecycle", "RetainedFor", "RemovalPlan", "Notes") %in% names(alias_tbl)))
  expect_identical(anyDuplicated(alias_tbl$Alias), 0L)

  expect_identical(
    alias_tbl$PreferredName[match("mfrmRFacets", alias_tbl$Alias)],
    "run_mfrm_facets"
  )
  expect_identical(
    alias_tbl$PreferredName[match("analyze_dif", alias_tbl$Alias)],
    "analyze_dff"
  )
  expect_identical(
    alias_tbl$PreferredName[match("JMLE", alias_tbl$Alias)],
    "JML"
  )
  expect_identical(
    alias_tbl$PreferredName[match("ReadyForAPA", alias_tbl$Alias)],
    "DraftReady"
  )

  removed_review_names <- c(
    "audit_mfrm_anchors",
    "precision_audit_report",
    "audit_conquest_overlap",
    "facet_small_sample_audit",
    "build_weighting_audit",
    "reference_case_audit"
  )
  expect_false(any(removed_review_names %in% alias_tbl$Alias))

  removed_field_names <- c(
    "config$anchor_audit",
    "diagnostics$precision_audit",
    "summary(conquest_overlap)$audit_scope"
  )
  expect_false(any(removed_field_names %in% alias_tbl$Alias))
})

test_that("compatibility_alias_table supports surface filters", {
  function_aliases <- compatibility_alias_table("functions")
  expect_true(nrow(function_aliases) >= 2L)
  expect_true(all(function_aliases$Surface == "function"))
  expect_true(all(c(
    "mfrmRFacets",
    "analyze_dif"
  ) %in% function_aliases$Alias))

  column_aliases <- compatibility_alias_table("columns")
  expect_true(all(column_aliases$Surface == "column"))
  expect_true(all(c("ReadyForAPA", "SE", "Fair(M) Average", "Fair(Z) Average") %in% column_aliases$Alias))

  argument_aliases <- compatibility_alias_table("arguments")
  expect_true(all(argument_aliases$Surface == "argument"))
  expect_true("JMLE" %in% argument_aliases$Alias)
  expect_false("build_model_choice_review(run_weighting_audit)" %in% argument_aliases$Alias)

  field_aliases <- compatibility_alias_table("fields")
  expect_true(all(field_aliases$Surface == "field"))
  expect_identical(nrow(field_aliases), 0L)

  plot_metric_aliases <- compatibility_alias_table("plot_metrics")
  expect_true(all(plot_metric_aliases$Surface == "plot_metric"))
  expect_true(all(c("FairM", "FairZ") %in% plot_metric_aliases$Alias))
})
