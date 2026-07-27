test_that("declared operational roster separates planned omissions from unassigned cells", {
  dat <- load_mfrmr_data("example_operational")
  env <- new.env(parent = emptyenv())
  utils::data(
    list = "mfrmr_example_operational_design",
    package = "mfrmr",
    envir = env
  )

  review <- describe_mfrm_data(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    expected_design = env$mfrmr_example_operational_design
  )
  structural <- review$structural_missingness$summary

  expect_equal(structural$Status, "declared")
  expect_equal(structural$ExpectedCells, 288L)
  expect_equal(structural$ObservedCells, 282L)
  expect_equal(structural$MatchedCells, 282L)
  expect_equal(structural$MissingExpectedCells, 6L)
  expect_equal(structural$UnexpectedObservedCells, 0L)
  expect_equal(structural$CoverageRate, 282 / 288)
  expect_equal(nrow(review$structural_missingness$missing_expected_cells), 6L)
  expect_equal(nrow(review$structural_missingness$unexpected_observed_cells), 0L)
  expect_true(all(review$design_connectivity$Connected))
  expect_setequal(review$design_connectivity$Basis, c("observed", "declared_expected"))
  expect_false(any(nzchar(review$design_components$PersonLabels)))

  compact <- summary(review)
  expect_equal(compact$structural_missingness$MissingExpectedCells, 6L)
  expect_true("planned_cells_not_observed" %in% compact$caveats$Condition)
})

test_that("structural missingness is not inferred without a declared roster", {
  dat <- load_mfrmr_data("example_operational")
  review <- describe_mfrm_data(
    dat,
    "Person",
    c("Rater", "Criterion"),
    "Score"
  )

  structural <- review$structural_missingness$summary
  expect_equal(structural$Status, "not_declared")
  expect_true(is.na(structural$ExpectedCells))
  expect_true(is.na(structural$MissingExpectedCells))
  expect_match(review$structural_missingness$settings$note, "not assessed")
  expect_false("planned_cells_not_observed" %in% summary(review)$caveats$Condition)
})

test_that("disconnected Person-facet components are explicit and privacy-safe", {
  dat <- data.frame(
    Person = sprintf("P%02d", 1:12),
    Rater = rep(c("R1", "R2"), each = 6L),
    Score = rep(c(1L, 2L), 6L),
    stringsAsFactors = FALSE
  )
  review <- describe_mfrm_data(
    dat,
    person = "Person",
    facets = "Rater",
    score = "Score",
    include_agreement = FALSE
  )

  observed <- review$design_connectivity[
    review$design_connectivity$Basis == "observed", , drop = FALSE
  ]
  expect_false(observed$Connected)
  expect_equal(observed$Components, 2L)
  expect_equal(nrow(review$design_components), 2L)
  expect_setequal(review$design_components$LevelLabels, c("R1", "R2"))
  expect_true(all(review$design_components$PersonLabels == ""))

  compact <- summary(review)
  expect_true(any(grepl("disconnected_observed", compact$caveats$Condition)))
  printed <- paste(capture.output(print(compact)), collapse = "\n")
  expect_false(grepl("P01|P02|P07|P08", printed))
})

test_that("declared roster validation and duplicate-cell audit are strict", {
  dat <- data.frame(
    Person = c("P01", "P01", sprintf("P%02d", 2:11)),
    Rater = rep(c("R1", "R2"), each = 6L),
    Score = rep(c(1L, 2L), 6L),
    stringsAsFactors = FALSE
  )
  valid_design <- unique(dat[c("Person", "Rater")])

  expect_error(
    describe_mfrm_data(dat, "Person", "Rater", "Score", expected_design = valid_design["Person"]),
    "missing required column"
  )
  blank_design <- valid_design
  blank_design$Rater[1L] <- " "
  expect_error(
    describe_mfrm_data(dat, "Person", "Rater", "Score", expected_design = blank_design),
    "missing or blank IDs"
  )
  duplicate_design <- rbind(valid_design, valid_design[1L, , drop = FALSE])
  expect_error(
    describe_mfrm_data(dat, "Person", "Rater", "Score", expected_design = duplicate_design),
    "duplicate planned cells"
  )
  expect_error(
    describe_mfrm_data(dat, "Person", "Rater", "Score", min_linking_persons = 0),
    "positive integer"
  )

  review <- suppressWarnings(describe_mfrm_data(
    dat,
    "Person",
    "Rater",
    "Score",
    expected_design = valid_design,
    include_agreement = FALSE,
    min_linking_persons = 7L
  ))
  expect_equal(review$duplicate_cell_summary$DuplicateCells, 1L)
  expect_equal(review$duplicate_cell_summary$RowsInDuplicateCells, 2L)
  expect_equal(review$duplicate_cell_summary$ExtraRows, 1L)
  expect_true(review$duplicate_cell_summary$HasDuplicates)
  expect_equal(nrow(review$duplicate_cell_detail), 1L)
  expect_gt(sum(review$linkage_summary$SparseLevels), 0L)
  expect_true("duplicate_person_facet_cells" %in% summary(review)$caveats$Condition)
})
