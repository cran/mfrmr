make_resampling_validation_data <- function() {
  people <- sprintf("P%02d", 1:10)
  regions <- c(rep("A", 6), rep("B", 3), "C")
  raters <- c("R1", "R2", "R3", "R4", "R1", "R2", "R3", "R4", "R1", "R2")
  rows <- lapply(seq_along(people), function(i) {
    data.frame(
      Study = "ObservedToy",
      Person = people[i],
      Region = regions[i],
      Rater = raters[i],
      Criterion = c("C1", "C2"),
      Score = c((i %% 4) + 1L, ((i + 1L) %% 4) + 1L),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

test_that("observed-data resampling specs preserve strata and rater coverage", {
  dat <- make_resampling_validation_data()

  spec <- build_mfrm_resampling_spec(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    strata = "Region",
    preserve_facets = "Rater",
    reps = 3,
    sample_n = 1,
    seed = 20260525
  )

  expect_s3_class(spec, "mfrm_resampling_spec")
  expect_equal(spec$design, "stratified_subsample")
  expect_true(all(c("A", "B", "C") %in% spec$stratum_overview$Stratum))
  expect_equal(spec$target_plan$TargetPersons, c(1L, 1L, 1L))

  spec_summary <- summary(spec)
  expect_s3_class(spec_summary, "summary.mfrm_resampling_spec")
  expect_true(all(c("overview", "stratum_overview", "target_plan",
                    "preserve_overview", "terminology") %in%
                    names(spec_summary)))
  expect_match(paste(spec_summary$terminology, collapse = " "), "not known true parameters")
  expect_output(print(spec_summary), "Observed-Data Resampling Specification")
  expect_output(print(spec), "Stratum target plan")

  draws <- draw_mfrm_resamples(spec)
  expect_s3_class(draws, "mfrm_resamples")
  expect_length(draws$samples, 3)
  expect_equal(nrow(draws$manifest), 3)
  expect_true(all(draws$manifest$StrataRepresented == 3L))
  expect_true(all(draws$manifest$PreserveCoverageComplete))
  expect_true(any(draws$manifest$TopupPersonClusters > 0L))
  expect_true(all(c(".mfrm_original_person", ".mfrm_draw_unit",
                    ".mfrm_stratum", ".mfrm_topup") %in%
                    names(draws$samples[[1]])))
  expect_true(setequal(unique(draws$samples[[1]]$Region), c("A", "B", "C")))
  expect_true(setequal(unique(draws$samples[[1]]$Rater), c("R1", "R2", "R3", "R4")))

  draws_summary <- summary(draws)
  expect_s3_class(draws_summary, "summary.mfrm_resamples")
  expect_equal(draws_summary$overview$Reps[1], 3L)
  expect_equal(draws_summary$overview$PreserveCoverageCompleteRate[1], 1)
  expect_true(nrow(draws_summary$stratum_summary) == 3L)
  expect_true(nrow(draws_summary$preserve_summary) == 1L)
  expect_output(print(draws_summary), "Observed-Data Resampling Draws Summary")
  expect_output(print(draws), "Preserve-facet summary")
})

test_that("observed-data resampling draws are reproducible from the spec seed", {
  dat <- make_resampling_validation_data()
  spec <- build_mfrm_resampling_spec(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    strata = "Region",
    preserve_facets = "Rater",
    reps = 2,
    sample_n = 1,
    seed = 77
  )

  first <- draw_mfrm_resamples(spec)
  second <- draw_mfrm_resamples(spec)

  expect_equal(first$manifest, second$manifest)
  expect_equal(first$stratum_manifest, second$stratum_manifest)
  expect_equal(first$preserve_manifest, second$preserve_manifest)
  expect_equal(first$samples, second$samples)
})

test_that("bootstrap resampling re-keys person clusters and allows manifest-only draws", {
  dat <- make_resampling_validation_data()
  spec <- build_mfrm_resampling_spec(
    dat,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    strata = "Region",
    preserve_facets = NULL,
    design = "stratified_bootstrap",
    reps = 1,
    sample_fraction = 1,
    seed = 123
  )

  draws <- draw_mfrm_resamples(spec)
  expect_true(draws$manifest$Replace[1])
  expect_equal(nrow(draws$preserve_manifest), 0L)
  expect_true(all(grepl("__b[0-9]{3}$", unique(draws$samples[[1]]$Person))))
  expect_true(".mfrm_original_person" %in% names(draws$samples[[1]]))

  manifest_only <- draw_mfrm_resamples(spec, keep_data = FALSE)
  expect_length(manifest_only$samples, 0L)
  expect_equal(nrow(manifest_only$manifest), 1L)
  expect_equal(nrow(manifest_only$preserve_manifest), 0L)
})

test_that("observed-data resampling rejects conflicting person strata", {
  dat <- make_resampling_validation_data()
  dat$Region[dat$Person == "P01" & dat$Criterion == "C2"] <- "B"

  expect_error(
    build_mfrm_resampling_spec(
      dat,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      strata = "Region"
    ),
    "at most one stratum"
  )
})
