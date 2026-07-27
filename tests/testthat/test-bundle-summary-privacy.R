test_that("bundle summaries suppress unexpected-response person identifiers by default", {
  private_person_id <- "PRIVATE_PERSON_ID_7XZ"
  bundle <- structure(
    list(
      table = data.frame(
        Row = 12L,
        Person = private_person_id,
        Facet = "Rater",
        Level = "R1",
        StdResidual = 3.2,
        stringsAsFactors = FALSE
      ),
      summary = data.frame(
        TotalObservations = 20L,
        UnexpectedN = 1L,
        stringsAsFactors = FALSE
      ),
      thresholds = list(abs_z_min = 2, prob_max = 0.3, rule = "either")
    ),
    class = c("mfrm_unexpected", "mfrm_bundle", "list")
  )

  default_summary <- summary(bundle)
  expect_false("Person" %in% names(default_summary$preview))
  expect_false(any(grepl(private_person_id, unlist(default_summary), fixed = TRUE)))
  expect_false(grepl(
    private_person_id,
    paste(capture.output(print(default_summary)), collapse = "\n"),
    fixed = TRUE
  ))

  person_summary <- summary(bundle, include_person = TRUE)
  expect_identical(person_summary$preview$Person[1], private_person_id)
  expect_match(
    paste(capture.output(print(person_summary)), collapse = "\n"),
    private_person_id,
    fixed = TRUE
  )
  expect_identical(bundle$table$Person[1], private_person_id)
})

test_that("bundle summaries redact person-facet levels without removing aggregates", {
  private_person_id <- "PRIVATE_FACET_LEVEL_4QR"
  bundle <- structure(
    list(
      summary = data.frame(Rows = 2L, MeanInfit = 1.05),
      table = data.frame(
        Facet = c("Person", "Rater"),
        Level = c(private_person_id, "R1"),
        Infit = c(1.2, 0.9),
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bundle", "list")
  )

  default_summary <- summary(bundle)
  expect_equal(default_summary$summary$Rows[1], 2L)
  expect_equal(default_summary$summary$MeanInfit[1], 1.05)
  expect_identical(default_summary$preview$Level[1], "<suppressed>")
  expect_identical(default_summary$preview$Level[2], "R1")
  expect_false(any(grepl(private_person_id, unlist(default_summary), fixed = TRUE)))

  person_summary <- summary(bundle, include_person = TRUE)
  expect_identical(person_summary$preview$Level[1], private_person_id)
})

test_that("bundle summaries suppress local directory details but retain source provenance", {
  private_root <- "/Users/private-user/restricted-study"
  private_file <- file.path(private_root, "conquest_parameters.csv")
  bundle <- structure(
    list(
      table = data.frame(
        Artifact = "parameters",
        Path = private_file,
        stringsAsFactors = FALSE
      ),
      settings = data.frame(
        Setting = c("parameter_file", "delimiter"),
        Value = c(private_file, ","),
        stringsAsFactors = FALSE
      )
    ),
    class = c("mfrm_bundle", "list")
  )

  summarized <- summary(bundle)
  printed <- paste(capture.output(print(summarized)), collapse = "\n")
  expect_false(grepl(private_root, paste(unlist(summarized), collapse = "\n"), fixed = TRUE))
  expect_false(grepl(private_root, printed, fixed = TRUE))
  expect_true(any(grepl("conquest_parameters.csv", unlist(summarized), fixed = TRUE)))
  expect_identical(bundle$table$Path[1], private_file)
  expect_identical(bundle$settings$Value[1], private_file)
})
