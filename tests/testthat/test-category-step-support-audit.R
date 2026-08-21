.category_support_data <- function(c2_categories = 0:3, categories = 0:3) {
  data <- expand.grid(
    Person = paste0("P", sprintf("%02d", 1:24)),
    Rater = c("R1", "R2"),
    Criterion = c("C1", "C2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  index <- as.integer(factor(data$Person)) + as.integer(factor(data$Rater))
  data$Score <- ifelse(
    data$Criterion == "C1",
    categories[(index %% length(categories)) + 1L],
    c2_categories[(index %% length(c2_categories)) + 1L]
  )
  data
}

.category_support_audit <- function(data, model = "PCM", rating_min = 0L,
                                    rating_max = 3L) {
  prep <- mfrmr:::prepare_mfrm_data(
    data,
    person_col = "Person",
    facet_cols = c("Rater", "Criterion"),
    score_col = "Score",
    rating_min = rating_min,
    rating_max = rating_max,
    keep_original = TRUE
  )
  step_facet <- if (identical(model, "RSM")) NULL else "Criterion"
  slope_facet <- if (identical(model, "GPCM")) "Criterion" else NULL
  signs <- mfrmr:::build_facet_signs(prep$facet_names, character(0))
  config <- mfrmr:::build_estimation_config(
    prep = prep,
    model = model,
    method = "JML",
    step_facet = step_facet,
    slope_facet = slope_facet,
    interaction_specs = NULL,
    weight_col = NULL,
    facet_signs = signs$signs,
    positive_facets = signs$positive_facets,
    noncenter_facet = "Person",
    dummy_facets = character(0),
    anchor_df = NULL,
    group_anchor_df = NULL,
    population = NULL
  )
  mfrmr:::audit_mfrm_category_support(
    prep = prep,
    config = config$config,
    sizes = config$sizes
  )
}

test_that("balanced PCM records the complete category and step contract", {
  audit <- .category_support_audit(.category_support_data())
  required <- c(
    "ScaleScope", "StepScope", "DeclaredCategories", "ObservedGlobal",
    "ObservedWithinScope", "RetainedForFit", "FreeStepCount",
    "FixedStepCount", "DerivedStepCount", "UnsupportedFreeStepCount",
    "UnsupportedCategory", "UnsupportedStep", "ZeroType",
    "MinimumObservedCategoryCount", "MaximumCategoryFraction",
    "NormalizedCategoryEntropy", "InformationState", "ReasonCode"
  )

  expect_identical(
    audit$contract_version,
    "mfrmr-readiness-0.2.3-v3"
  )
  expect_identical(audit$readiness$CategoryState, "adequate")
  expect_true(all(required %in% names(audit$support_table)))
  expect_identical(nrow(audit$support_table), 2L)
  expect_identical(nrow(audit$category_table), 8L)
  expect_identical(nrow(audit$step_status), 6L)
  expect_true(all(audit$support_table$FreeStepCount == 2L))
  expect_true(all(audit$support_table$FixedStepCount == 0L))
  expect_true(all(audit$support_table$DerivedStepCount == 1L))
  expect_true(all(audit$support_table$UnsupportedStep == ""))
  expect_true(all(audit$step_status$ParameterStatus == "estimable"))
  expect_false(any(audit$step_status$Fixed))
  expect_true(all(audit$category_table$RetainedForFit))

  fit <- suppressWarnings(fit_mfrm(
    .category_support_data(),
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    rating_min = 0,
    rating_max = 3,
    model = "PCM",
    method = "JML",
    step_facet = "Criterion",
    maxit = 5
  ))
  expect_identical(
    fit$data_review$category_support,
    fit$config$category_support_audit
  )
  expect_identical(
    fit$data_review$category_support$readiness$CategoryState,
    "adequate"
  )
})

test_that("PCM unsupported internal-category contrast stops before optimization", {
  data <- .category_support_data(c2_categories = c(0L, 2L, 3L))
  audit <- .category_support_audit(data)
  c2 <- audit$support_table[audit$support_table$StepScope == "Criterion=C2", ]
  affected <- audit$step_status[
    audit$step_status$StepScope == "Criterion=C2" &
      audit$step_status$ParameterStatus == "unsupported", , drop = FALSE
  ]

  expect_identical(audit$readiness$CategoryState, "unsupported_coordinate")
  expect_identical(c2$ObservedGlobal, "0;1;2;3")
  expect_identical(c2$ObservedWithinScope, "0;2;3")
  expect_identical(c2$RetainedForFit, "0;1;2;3")
  expect_identical(c2$UnsupportedCategory, "1")
  expect_identical(c2$UnsupportedStep, "Step1;Step2")
  expect_identical(c2$UnsupportedFreeStepCount, 2L)
  expect_identical(nrow(audit$unsupported_contrasts), 1L)
  expect_identical(audit$unsupported_contrasts$Direction, "Step1 - Step2")
  expect_identical(audit$unsupported_contrasts$FreeCoordinates,
                   "FreeStep1;FreeStep2")
  expect_identical(affected$Step, c("Step1", "Step2"))
  expect_true(all(affected$CrossingSupport))
  expect_true(all(affected$FreeCoordinateAffected))

  condition <- tryCatch(
    suppressWarnings(fit_mfrm(
      data,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      rating_min = 0,
      rating_max = 3,
      model = "PCM",
      method = "JML",
      step_facet = "Criterion",
      maxit = 5
    )),
    error = identity
  )
  expect_s3_class(condition, "mfrmr_category_readiness_error")
  expect_s3_class(condition, "mfrmr_readiness_condition")
  expect_identical(condition$readiness$CategoryState,
                   "unsupported_coordinate")
  expect_match(condition$readiness$ReasonCodes,
               "step_scope_category_unobserved", fixed = TRUE)
  expect_match(condition$readiness$ReasonCodes,
               "unsupported_step_coordinate", fixed = TRUE)
  expect_false(grepl("P01", conditionMessage(condition), fixed = TRUE))
})

test_that("RSM local absence is review evidence, not a common-step blocker", {
  audit <- .category_support_audit(
    .category_support_data(c2_categories = 0:2),
    model = "RSM"
  )
  c2 <- audit$local_support[
    audit$local_support$Facet == "Criterion" &
      audit$local_support$Level == "C2", , drop = FALSE
  ]

  expect_identical(audit$readiness$CategoryState, "weak_information")
  expect_identical(audit$support_table$InformationState, "adequate")
  expect_identical(c2$MissingCategories, "3")
  expect_identical(c2$InformationState, "weak_information")
  expect_length(audit$unsupported_step_coordinates, 0L)
  expect_warning(
    mfrmr:::mfrmr_warn_category_support(audit),
    class = "mfrmr_category_readiness_warning"
  )
})

test_that("empty intermediate categories identify exact step recession directions", {
  data <- .category_support_data(
    c2_categories = c(0L, 2L, 3L),
    categories = c(0L, 2L, 3L)
  )
  audit <- .category_support_audit(data)

  expect_identical(audit$readiness$CategoryState, "unsupported_coordinate")
  expect_true(any(audit$step_status$ParameterStatus == "unsupported"))
  expect_true(all(audit$support_table$UnsupportedCategory == "1"))
  expect_true(all(audit$support_table$UnsupportedStep == "Step1;Step2"))
  expect_true(any(audit$category_table$ZeroType ==
                    "sampling_zero_global_internal"))
  expect_identical(nrow(audit$unsupported_contrasts), 2L)
})

test_that("unsupported contrast suppresses only its missing internal category", {
  eta <- 0.35
  steps <- c(-0.6, -0.1, 0.2, 0.5)
  direction <- c(0, 1, -1, 0)
  log_weights <- function(step_values) {
    0:4 * eta - c(0, cumsum(step_values))
  }

  shift <- log_weights(steps + 7 * direction) - log_weights(steps)
  expect_equal(shift, c(0, 0, -7, 0, 0), tolerance = 1e-12)

  observed <- c(0L, 1L, 3L, 4L)
  log_likelihood <- function(multiplier) {
    values <- log_weights(steps + multiplier * direction)
    normalizer <- max(values) + log(sum(exp(values - max(values))))
    sum(values[observed + 1L] - normalizer)
  }
  expect_gt(log_likelihood(7), log_likelihood(0))
  expect_gt(log_likelihood(20), log_likelihood(7))
})

test_that("bounded GPCM inherits the PCM unsupported-step gate", {
  audit <- .category_support_audit(
    .category_support_data(c2_categories = c(0L, 2L, 3L)),
    model = "GPCM"
  )

  expect_identical(audit$model, "GPCM")
  expect_identical(audit$readiness$CategoryState, "unsupported_coordinate")
  expect_identical(audit$support_table$UnsupportedCategory[2], "1")
  expect_identical(audit$support_table$UnsupportedStep[2], "Step1;Step2")
  expect_false(any(audit$step_status$ParameterClass != "step"))
})

test_that("unobserved boundary categories remain review evidence for WP3", {
  audit <- .category_support_audit(
    .category_support_data(c2_categories = 0:2)
  )
  c2 <- audit$support_table[audit$support_table$StepScope == "Criterion=C2", ]

  expect_identical(audit$readiness$CategoryState, "weak_information")
  expect_identical(c2$UnsupportedCategory, "")
  expect_identical(c2$UnsupportedStep, "")
  expect_identical(c2$UnsupportedFreeStepCount, 0L)
  expect_length(audit$unsupported_step_coordinates, 0L)
  expect_identical(nrow(audit$unsupported_contrasts), 0L)
  expect_match(c2$ZeroType, "sampling_zero_step_scope_boundary", fixed = TRUE)
})

test_that("binary ladders do not invent unsupported free step coordinates", {
  data <- .category_support_data(
    c2_categories = 0L,
    categories = 0:1
  )
  audit <- .category_support_audit(
    data,
    model = "PCM",
    rating_min = 0L,
    rating_max = 1L
  )
  c2 <- audit$support_table[audit$support_table$StepScope == "Criterion=C2", ]
  c2_step <- audit$step_status[
    audit$step_status$StepScope == "Criterion=C2", , drop = FALSE
  ]

  expect_identical(audit$readiness$CategoryState, "weak_information")
  expect_identical(c2$FreeStepCount, 0L)
  expect_identical(c2$UnsupportedFreeStepCount, 0L)
  expect_identical(c2$UnsupportedStep, "")
  expect_identical(c2_step$ParameterStatus, "not_estimated")
  expect_false(c2_step$FreeCoordinateAffected)
})

test_that("category support decisions are invariant to row order and labels", {
  data <- .category_support_data(c2_categories = c(0L, 2L, 3L))
  reversed <- data[rev(seq_len(nrow(data))), , drop = FALSE]
  relabelled <- data
  relabelled$Criterion <- ifelse(
    relabelled$Criterion == "C1", "Scale_B", "Scale_A"
  )
  relabelled$Rater <- ifelse(relabelled$Rater == "R1", "Judge_Y", "Judge_X")

  original_audit <- .category_support_audit(data)
  reversed_audit <- .category_support_audit(reversed)
  relabelled_audit <- .category_support_audit(relabelled)

  summarize <- function(audit) {
    table(
      audit$support_table$ObservedWithinScope,
      audit$support_table$InformationState,
      audit$support_table$UnsupportedStep
    )
  }
  expect_identical(reversed_audit$readiness$CategoryState,
                   original_audit$readiness$CategoryState)
  expect_identical(relabelled_audit$readiness$CategoryState,
                   original_audit$readiness$CategoryState)
  expect_identical(summarize(reversed_audit), summarize(original_audit))
  expect_identical(summarize(relabelled_audit), summarize(original_audit))
})
