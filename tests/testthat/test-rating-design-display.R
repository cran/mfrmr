test_that("observed rater panels retain shared-person linkage separately", {
  ratings <- expand.grid(Person = paste0("P", 1:8), Rater = paste0("R", 1:4),
                         Criterion = c("C1", "C2"), stringsAsFactors = FALSE)
  ratings$Score <- (match(ratings$Person, unique(ratings$Person)) +
                    match(ratings$Rater, unique(ratings$Rater))) %% 4L
  same_panel <- (ratings$Person %in% paste0("P", 1:4)) == (ratings$Rater %in% c("R1", "R2"))
  audit <- function(d) check_mfrm_connectivity(d, "Person", c("Rater", "Criterion"), "Score")
  separated <- audit(ratings[same_panel, ])
  expect_true(separated$connected)
  expect_equal(separated$rater_linkage$n_components, 2L)
  expect_setequal(separated$rater_linkage$components$LevelLabels, c("R1, R2", "R3, R4"))
  bridge <- audit(ratings[same_panel | ratings$Person == "P1", ])
  expect_equal(bridge$rater_linkage$n_components, 1L)
  expect_equal(bridge$rater_linkage$bridge_persons, "P1")
  redundant <- audit(ratings[same_panel | ratings$Person %in% c("P1", "P2"), ])
  expect_length(redundant$rater_linkage$bridge_persons, 0L)
})

test_that("saved-calibration information validates and counts the declared rating plan", {
  data <- load_mfrmr_data("example_core")
  fit <- read_example_bundle()$fit
  context <- data[1L, c("Rater", "Criterion"), drop = FALSE]
  compute <- function(...) compute_information(fit, c(-1, 1), 3L, ...)
  one <- compute(rating_plan = transform(context, Count = 1L))
  two <- compute(rating_plan = transform(context, Count = 2L))
  expect_equal(two$test_information$Information, 2 * one$test_information$Information)
  expect_equal(two$test_information$SE, one$test_information$SE / sqrt(2))
  expect_error(compute(rating_plan = transform(context, Rater = "unknown", Count = 1L)),
               "Unknown or missing rating-plan level")
  p <- predict(fit, transform(context, theta = 0), type = "probs")
  scores <- as.numeric(colnames(p))
  expect_equal(one$test_information$Information[2L], sum(p * (scores - sum(p * scores))^2))
  average <- compute(exposure_basis = "per_person_average")
  expect_equal(sum(average$exposures$Count), nrow(data) / length(unique(data$Person)))
  expect_identical(compute()$exposure_basis, "unique_contexts")
  expect_equal(sum(compute(exposure_basis = "observed_total")$exposures$Count), nrow(data))
})
