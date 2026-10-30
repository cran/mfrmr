hierarchical_identity_fixture <- local({
  fits <- readRDS(test_path("fixtures", "mml-comparison-intervals.rds"))
  fit <- fits$RSM$fit
  list(fit = fit, diagnostics = diagnose_mfrm(fit, residual_pca = "none"),
       review = analyze_hierarchical_structure(fit, facets = c("Rater", "Criterion"),
         compute_icc = FALSE, igraph_layout = FALSE))
})

test_that("the reporting checklist verifies reviewed data and names its scope", {
  x <- hierarchical_identity_fixture
  local_mocked_bindings(analyze_hierarchical_structure = function(...) stop("Unexpected review"),
    diagnose_mfrm = function(...) stop("Unexpected diagnostics"), .package = "mfrmr")
  consume <- function(review) reporting_checklist(x$fit, diagnostics = x$diagnostics,
                                                  hierarchical_structure = review)
  out <- consume(x$review)$checklist
  row <- out[out$Item == "Hierarchical structure review", ]
  expect_true(row$Available)
  expect_true(row$DraftReady)
  expect_match(row$Detail, "Rater, Criterion (person column: Person)", fixed = TRUE)
  expect_match(row$Detail, "384 retained rating rows", fixed = TRUE)
  expect_match(row$Detail, "ICC not available", fixed = TRUE)
  old <- x$review
  old$data_identity <- NULL
  expect_error(consume(old), "no MFRM refit is needed")
  expect_output(summary(old), "Summary")
  expect_error(consume(list()), "must come from analyze_hierarchical_structure")
  altered <- x$review
  altered$facets <- c("Rater", "Group")
  expect_error(consume(altered), "fit/hierarchical_structure mismatch")
})

test_that("structural review follows selected observations rather than fitted estimates", {
  x <- hierarchical_identity_fixture
  data <- x$fit$prep$data
  reordered <- data[rev(seq_len(nrow(data))), ]
  reordered$Rater <- factor(reordered$Rater, levels = rev(levels(reordered$Rater)))
  reordered$Unrelated <- seq_len(nrow(data))
  reordered$Score <- reordered$Score + 1
  reordered$Weight <- 10
  review <- analyze_hierarchical_structure(reordered, c("Rater", "Criterion"),
    compute_icc = FALSE, igraph_layout = FALSE)
  expect_identical(review$data_identity, x$review$data_identity)
  expect_no_error(reporting_checklist(x$fit, diagnostics = x$diagnostics,
                                     hierarchical_structure = review))
  altered <- x$fit
  altered$facets$others$Estimate <- altered$facets$others$Estimate + .1
  expect_no_error(.ha_validate_data_identity(altered, review))
  for (change in list(
    function(d) { d$Person <- as.character(d$Person); d$Person[1] <- "OTHER"; d },
    function(d) { d$Rater[1] <- d$Rater[2]; d },
    function(d) d[-1, ],
    function(d) rbind(d[-1, ], d[2, ])
  )) {
    other <- analyze_hierarchical_structure(change(data), c("Rater", "Criterion"),
      compute_icc = FALSE, igraph_layout = FALSE)
    expect_error(reporting_checklist(x$fit, diagnostics = x$diagnostics,
      hierarchical_structure = other), "fit/hierarchical_structure mismatch")
  }
  # A review without Person only describes the selected facet assignments.
  no_person <- analyze_hierarchical_structure(data, c("Rater", "Criterion"),
    person = NULL, compute_icc = FALSE, igraph_layout = FALSE)
  expect_no_error(.ha_validate_data_identity(x$fit, no_person))
  expect_null(no_person$data_identity$person)
})

test_that("raw-data column mappings and saved structural reviews remain reusable", {
  x <- hierarchical_identity_fixture
  data <- x$fit$prep$data
  names(data)[names(data) == "Person"] <- "Examinee"
  data$Rater <- as.character(data$Rater)
  review <- analyze_hierarchical_structure(data, c("Rater", "Criterion"),
    person = "Examinee", compute_icc = FALSE, igraph_layout = FALSE)
  mapped <- x$fit
  mapped$config$source_columns$person <- "Examinee"
  expect_no_error(.ha_validate_data_identity(mapped, review))
  expect_error(.ha_validate_data_identity(x$fit, review), "reviewed columns or rating rows")
  data$School <- data$Group
  outside <- analyze_hierarchical_structure(data, c("Rater", "School"),
    person = "Examinee", compute_icc = FALSE, igraph_layout = FALSE)
  expect_error(.ha_validate_data_identity(mapped, outside), "reviewed columns or rating rows")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(list(fit = mapped, review = review), path)
  saved <- readRDS(path)
  expect_no_error(.ha_validate_data_identity(saved$fit, saved$review))
})

test_that("ICC-containing reviews also verify the analyzed scores", {
  skip_if_not_installed("lme4")
  x <- hierarchical_identity_fixture
  data <- x$fit$prep$data
  names(data)[names(data) == "Score"] <- "Rating"
  data$Rating <- factor(data$Rating, levels = rev(sort(unique(data$Rating))))
  review <- analyze_hierarchical_structure(data, c("Rater", "Criterion"),
    score = "Rating", igraph_layout = FALSE)
  mapped <- x$fit
  mapped$config$source_columns$score <- "Rating"
  expect_no_error(.ha_validate_data_identity(mapped, review))
  expect_identical(review$data_identity$score, "Rating")
  precise <- mapped$prep$data
  precise$Score[1] <- precise$Score[1] + 1e-15
  expect_false(identical(
    .ha_data_identity(precise, c("Rater", "Criterion"), "Person", "Score"),
    .ha_data_identity(mapped$prep$data, c("Rater", "Criterion"), "Person", "Score")))
  changed <- mapped
  changed$prep$data$Score[1] <- changed$prep$data$Score[1] + 1
  expect_error(.ha_validate_data_identity(changed, review), "fit/hierarchical_structure mismatch")
  # Missing-value omission in ICC does not erase rows from the design review.
  data$Rating[1] <- NA
  omitted <- analyze_hierarchical_structure(data, c("Rater", "Criterion"),
    score = "Rating", missing = "omit", igraph_layout = FALSE)
  expect_identical(omitted$data_identity$rows, nrow(data))
  expect_error(.ha_validate_data_identity(mapped, omitted), "fit/hierarchical_structure mismatch")
})
