test_that("gpcm_capability_matrix exposes the bounded GPCM support contract", {
  tbl <- gpcm_capability_matrix()

  expect_s3_class(tbl, "data.frame")
  expect_true(all(c(
    "Area", "Helpers", "Status", "PrimaryUse", "Boundary", "Evidence"
  ) %in% names(tbl)))
  expect_true(all(tbl$Status %in% c(
    "supported", "supported_with_caveat", "blocked", "deferred"
  )))

  expect_true(any(
    tbl$Area == "Core fitting and summaries" &
      tbl$Status == "supported"
  ))
  expect_true(any(
    tbl$Area == "Exploratory diagnostics and residual follow-up" &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_apa_outputs\\(\\)", tbl$Helpers) &
      tbl$Status == "blocked"
  ))
  expect_true(any(
    tbl$Area == "Design planning and forecasting" &
      tbl$Status == "deferred"
  ))
  expect_true(any(
    grepl("build_misfit_casebook\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_weighting_audit\\(\\)", tbl$Helpers) &
      tbl$Status == "supported_with_caveat"
  ))
  expect_true(any(
    grepl("build_linking_review\\(\\)", tbl$Helpers) &
      tbl$Status == "deferred"
  ))
})

test_that("gpcm_capability_matrix filters by status", {
  full_tbl <- gpcm_capability_matrix()
  blocked_tbl <- gpcm_capability_matrix("blocked")
  supported_tbl <- gpcm_capability_matrix("supported")

  expect_true(nrow(blocked_tbl) > 0)
  expect_true(nrow(supported_tbl) > 0)
  expect_true(all(blocked_tbl$Status == "blocked"))
  expect_true(all(supported_tbl$Status == "supported"))
  expect_equal(nrow(blocked_tbl), sum(full_tbl$Status == "blocked"))
  expect_equal(nrow(supported_tbl), sum(full_tbl$Status == "supported"))
})
