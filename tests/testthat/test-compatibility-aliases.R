test_that("compatibility_alias_table exposes retained aliases and preferred names", {
  alias_tbl <- compatibility_alias_table()

  expect_true(is.data.frame(alias_tbl))
  expect_true(all(c("Alias", "PreferredName", "Surface", "Lifecycle", "RetainedFor", "Notes") %in% names(alias_tbl)))
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
})

test_that("compatibility_alias_table supports surface filters", {
  function_aliases <- compatibility_alias_table("functions")
  expect_true(nrow(function_aliases) >= 2L)
  expect_true(all(function_aliases$Surface == "function"))
  expect_true(all(c("mfrmRFacets", "analyze_dif") %in% function_aliases$Alias))

  column_aliases <- compatibility_alias_table("columns")
  expect_true(all(column_aliases$Surface == "column"))
  expect_true(all(c("ReadyForAPA", "SE", "Fair(M) Average", "Fair(Z) Average") %in% column_aliases$Alias))

  plot_metric_aliases <- compatibility_alias_table("plot_metrics")
  expect_true(all(plot_metric_aliases$Surface == "plot_metric"))
  expect_true(all(c("FairM", "FairZ") %in% plot_metric_aliases$Alias))
})
