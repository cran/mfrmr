fixed_calibration_fit_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(model = c("RSM", "PCM")) {
    model <- match.arg(model)
    if (exists(model, envir = cache, inherits = FALSE)) {
      return(get(model, envir = cache, inherits = FALSE))
    }
    toy <- load_mfrmr_data("example_core")
    keep <- unique(toy$Person)[seq_len(if (identical(model, "RSM")) 18L else 24L)]
    toy <- toy[toy$Person %in% keep, , drop = FALSE]
    fit <- suppressWarnings(fit_mfrm(
      toy,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = model,
      step_facet = if (identical(model, "PCM")) "Criterion" else NULL,
      quad_points = 5,
      maxit = if (identical(model, "RSM")) 15 else 40
    ))
    assign(model, list(fit = fit, data = toy), envir = cache)
    get(model, envir = cache, inherits = FALSE)
  }
})

fixed_calibration_draft_fixture <- function(model = c("RSM", "PCM")) {
  model <- match.arg(model)
  fixture <- fixed_calibration_fit_fixture(model)
  draft <- mfrmr:::mfrmr_extract_calibration_draft(
    fixture$fit,
    calibration_id = paste0("calibration-", tolower(model)),
    source_fit_id = paste0("fit-", tolower(model)),
    created_at_utc = "2026-08-22T00:00:00Z"
  )
  list(fit = fixture$fit, data = fixture$data, draft = draft)
}

fixed_calibration_frozen_fixture <- function(model = c("RSM", "PCM")) {
  fixture <- fixed_calibration_draft_fixture(match.arg(model))
  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    fixture$draft, validated_at_utc = "2026-08-22T00:01:00Z"
  )
  fixture$frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-22T00:02:00Z"
  )
  fixture
}

fixed_calibration_scoring_rows <- function(data, n_person = 3L) {
  source_persons <- unique(as.character(data$Person))[seq_len(n_person)]
  rows <- data[data$Person %in% source_persons,
               c("Person", "Rater", "Criterion", "Score"), drop = FALSE]
  rows$Person <- paste0("NEW", match(as.character(rows$Person), source_persons))
  rownames(rows) <- NULL
  rows
}

fixed_calibration_anchored_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    fixture <- fixed_calibration_fit_fixture("PCM")
    anchors <- data.frame(
      Facet = "Rater", Level = "R01", Anchor = 0,
      stringsAsFactors = FALSE
    )
    groups <- data.frame(
      Facet = c("Criterion", "Criterion"),
      Level = c("Accuracy", "Content"),
      Group = c("G1", "G1"),
      GroupValue = c(0.1, 0.1),
      stringsAsFactors = FALSE
    )
    fit <- suppressWarnings(fit_mfrm(
      fixture$data,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = "RSM",
      anchors = anchors,
      group_anchors = groups,
      anchor_policy = "error",
      quad_points = 5,
      maxit = 40
    ))
    draft <- mfrmr:::mfrmr_extract_calibration_draft(
      fit,
      calibration_id = "calibration-anchored-rsm",
      source_fit_id = "fit-anchored-rsm",
      created_at_utc = "2026-08-22T00:00:00Z"
    )
    cache <<- list(fit = fit, draft = draft)
    cache
  }
})

fixed_calibration_interaction_fixture <- local({
  cache <- list()
  function(model = c("RSM", "PCM")) {
    model <- match.arg(model)
    if (!is.null(cache[[model]])) return(cache[[model]])
    data <- mfrmr:::sample_mfrm_data(seed = 42)
    fit <- suppressWarnings(suppressMessages(fit_mfrm(
      data,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "MML",
      model = model,
      step_facet = if (identical(model, "PCM")) "Criterion" else NULL,
      facet_interactions = "Rater:Criterion",
      min_obs_per_interaction = 0,
      quad_points = 5,
      maxit = 80,
      mml_engine = "direct"
    )))
    draft <- mfrmr:::mfrmr_extract_calibration_draft(
      fit,
      calibration_id = "calibration-interaction-rsm",
      source_fit_id = "fit-interaction-rsm",
      created_at_utc = "2026-08-22T00:00:00Z"
    )
    cache[[model]] <<- list(fit = fit, draft = draft, data = data)
    cache[[model]]
  }
})

capture_calibration_error <- function(expr) {
  tryCatch(
    force(expr),
    mfrm_calibration_error = function(e) e
  )
}

test_that("RSM draft contains expanded calibration state and no fit payload", {
  fixture <- fixed_calibration_draft_fixture("RSM")
  draft <- fixture$draft

  expect_s3_class(draft, "mfrm_calibration")
  expect_identical(draft$lifecycle$state, "draft")
  expect_identical(draft$model$family, "RSM")
  expect_identical(draft$model$estimator, "MML")
  expect_identical(draft$model$step_owner, "shared")
  expect_identical(draft$scoring_basis$type, "fixed_standard_normal")
  expect_identical(draft$scoring_basis$prior_mean, 0)
  expect_identical(draft$scoring_basis$prior_sd, 1)
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(draft)), 0L)
  expect_true(all(c("facet", "shared_step") %in%
    draft$parameters$coordinates$ParameterClass))
  expect_false(any(draft$parameters$coordinates$ParameterClass == "owned_step"))
  expect_identical(
    mfrmr:::mfrmr_calibration_find_prohibited(draft), character(0)
  )
  expect_false(any(fixture$data$Person %in% unlist(draft, use.names = FALSE)))

  facet_coordinates <- draft$parameters$coordinates[
    draft$parameters$coordinates$ParameterClass == "facet",
    c("OwnerFacet", "Level", "Value")
  ]
  expected <- as.data.frame(fixture$fit$facets$others)[
    , c("Facet", "Level", "Estimate")
  ]
  expect_identical(facet_coordinates$OwnerFacet, as.character(expected$Facet))
  expect_identical(facet_coordinates$Level, as.character(expected$Level))
  expect_equal(facet_coordinates$Value, expected$Estimate, tolerance = 0)
  expect_equal(sum(draft$scoring_basis$weights), 1, tolerance = 1e-15)
  expect_true(all(diff(draft$scoring_basis$nodes) > 0))
})

test_that("PCM draft stores every owner-specific step at full precision", {
  fixture <- fixed_calibration_draft_fixture("PCM")
  draft <- fixture$draft
  owned <- draft$parameters$coordinates[
    draft$parameters$coordinates$ParameterClass == "owned_step", , drop = FALSE
  ]

  expect_identical(draft$model$family, "PCM")
  expect_identical(draft$model$step_owner, "Criterion")
  expect_identical(
    draft$eligibility$support_profile_id,
    "pcm_mml_fixed_standard_normal_v1"
  )
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(draft)), 0L)
  expect_identical(
    nrow(owned),
    as.integer(length(fixture$fit$config$facet_levels$Criterion) *
      (fixture$fit$config$n_cat - 1L))
  )
  expect_true(all(owned$OwnerFacet == "Criterion"))
  step_matrix <- matrix(
    owned$Value,
    nrow = length(fixture$fit$config$facet_levels$Criterion),
    byrow = TRUE
  )
  expect_equal(rowSums(step_matrix), rep(0, nrow(step_matrix)), tolerance = 1e-12)
})

test_that("normalized direct and group anchors become typed declarations", {
  fixture <- fixed_calibration_anchored_fixture()
  anchors <- fixture$draft$constraints$anchors

  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(fixture$draft)), 0L)
  expect_identical(
    anchors$AnchorType,
    c("direct", "group", "group")
  )
  expect_identical(anchors$OwnerFacet, c("Rater", "Criterion", "Criterion"))
  expect_identical(anchors$Level, c("R01", "Accuracy", "Content"))
  expect_equal(anchors$Value, c(0, 0.1, 0.1), tolerance = 0)
  expect_identical(anchors$CoordinateSystem, rep("expanded_logit", 3))
  expect_identical(anchors$DeclarationOrder, 1:3)

  review <- suppressWarnings(mml_quadrature_sensitivity(
    fixture$fit, fixed_calibration_fit_fixture("PCM")$data,
    quad_points = c(5L, 7L), theta_points = 41L
  ))
  reviewed <- extract_mfrm_calibration(review$fits$q7,
                                      quadrature_review = review)
  expect_identical(reviewed$constraints$anchors, anchors)
})

test_that("RSM interaction coordinates retain the complete identified matrix", {
  fixture <- fixed_calibration_interaction_fixture()
  draft <- fixture$draft
  interactions <- draft$model$interactions
  coordinates <- draft$parameters$coordinates[
    draft$parameters$coordinates$ParameterClass == "interaction", , drop = FALSE
  ]
  expected <- fixture$fit$interactions$effects

  expect_true(mfrmr:::mfrm_inference_ready(fixture$fit))
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(draft)), 0L)
  expect_identical(nrow(interactions), 1L)
  expect_identical(interactions$InteractionId, "Rater:Criterion")
  expect_identical(interactions$Identification, "two_way_sum_to_zero_margins")
  expect_identical(nrow(coordinates), nrow(expected))
  expect_equal(coordinates$Value, expected$Estimate, tolerance = 0)
  mat <- matrix(
    coordinates$Value,
    nrow = interactions$LevelCountA,
    ncol = interactions$LevelCountB
  )
  expect_equal(rowSums(mat), rep(0, nrow(mat)), tolerance = 1e-12)
  expect_equal(colSums(mat), rep(0, ncol(mat)), tolerance = 1e-12)
})

test_that("validation and freezing are immutable registered transitions", {
  draft <- fixed_calibration_draft_fixture("RSM")$draft
  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    draft, validated_at_utc = "2026-08-22T00:01:00Z"
  )
  frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-22T00:02:00Z"
  )

  expect_identical(draft$lifecycle$state, "draft")
  expect_identical(validated$lifecycle$state, "validated")
  expect_identical(frozen$lifecycle$state, "frozen")
  expect_identical(draft$lifecycle$revision, 1L)
  expect_identical(validated$lifecycle$revision, 2L)
  expect_identical(frozen$lifecycle$revision, 3L)
  expect_identical(
    frozen$lifecycle$events$Operation,
    c("extract_draft", "validate", "freeze")
  )
  expect_false(draft$validation$schema_valid)
  expect_true(validated$validation$schema_valid)
  expect_true(frozen$validation$semantic_valid)
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(frozen)), 0L)

  err <- capture_calibration_error(mfrmr:::mfrmr_freeze_calibration(draft))
  expect_s3_class(err, "mfrm_calibration_error")
  expect_identical(err$code, "LIFECYCLE_TRANSITION_INVALID")
})

test_that("canonical RDS persistence is lossless and refuses overwrite", {
  draft <- fixed_calibration_draft_fixture("RSM")$draft
  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    draft, validated_at_utc = "2026-08-22T00:01:00Z"
  )
  frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-22T00:02:00Z"
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  saved <- mfrmr:::mfrmr_save_calibration(frozen, path)
  expect_true(file.exists(saved))
  restored <- mfrmr:::mfrmr_load_calibration(path)
  expect_identical(restored, frozen)
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(restored)), 0L)

  err <- capture_calibration_error(
    mfrmr:::mfrmr_save_calibration(frozen, path)
  )
  expect_s3_class(err, "mfrm_calibration_error")
  expect_identical(err$code, "PERSISTENCE_TARGET_EXISTS")
})

test_that("mutations fail closed with stable refusal codes", {
  draft <- fixed_calibration_draft_fixture("RSM")$draft

  newer <- draft
  newer$header$schema_version <- 7L
  review <- mfrmr:::mfrmr_review_calibration(newer)
  expect_true("SCHEMA_VERSION_UNSUPPORTED" %in% review$Code)

  partial <- draft
  partial$scoring_basis$nodes <- NULL
  review <- mfrmr:::mfrmr_review_calibration(partial)
  expect_true("SCHEMA_FIELD_MISSING" %in% review$Code)

  duplicate <- draft
  duplicate$parameters$coordinates$CoordinateKey[2] <-
    duplicate$parameters$coordinates$CoordinateKey[1]
  review <- mfrmr:::mfrmr_review_calibration(duplicate)
  expect_true("PARAMETER_COORDINATE_DUPLICATE" %in% review$Code)

  nonfinite <- draft
  nonfinite$parameters$coordinates$Value[1] <- Inf
  review <- mfrmr:::mfrmr_review_calibration(nonfinite)
  expect_true("PARAMETER_NONFINITE" %in% review$Code)

  remapped <- draft
  remapped$response$score_map$OriginalScore[1] <- 999
  review <- mfrmr:::mfrmr_review_calibration(remapped)
  expect_true("IDENTITY_COMPONENT_MISMATCH" %in% review$Code)

  edge <- draft
  edge$scoring_basis$weights[1] <- 0
  review <- mfrmr:::mfrmr_review_calibration(edge)
  expect_true("QUADRATURE_WEIGHTS_INVALID" %in% review$Code)

  injected <- draft
  injected$provenance$training_design_matrix <- matrix(1, 1, 1)
  review <- mfrmr:::mfrmr_review_calibration(injected)
  expect_true("PROHIBITED_TRAINING_STATE" %in% review$Code)
  expect_true("SCHEMA_FIELD_UNEXPECTED" %in% review$Code)

  err <- capture_calibration_error(
    mfrmr:::mfrmr_validate_calibration_draft(nonfinite)
  )
  expect_s3_class(err, "mfrm_calibration_error")
  expect_identical(err$code, "PARAMETER_NONFINITE")
  expect_true(nrow(err$refusals) >= 1L)
})

test_that("load rejects corrupt, partial, altered, and newer artifacts", {
  draft <- fixed_calibration_draft_fixture("RSM")$draft
  paths <- vapply(seq_len(4), function(i) tempfile(fileext = ".rds"), character(1))
  on.exit(unlink(paths), add = TRUE)

  writeBin(charToRaw("not an RDS artifact"), paths[1])
  err <- capture_calibration_error(mfrmr:::mfrmr_load_calibration(paths[1]))
  expect_s3_class(err, "mfrm_calibration_error")
  expect_identical(err$code, "PERSISTENCE_READ_FAILED")

  partial <- draft
  partial$parameters$coordinates <- NULL
  saveRDS(partial, paths[2], version = 3)
  err <- capture_calibration_error(mfrmr:::mfrmr_load_calibration(paths[2]))
  expect_identical(err$code, "SCHEMA_FIELD_MISSING")

  altered <- draft
  altered$response$score_map$OriginalScore[1] <- -999
  saveRDS(altered, paths[3], version = 3)
  err <- capture_calibration_error(mfrmr:::mfrmr_load_calibration(paths[3]))
  expect_identical(err$code, "IDENTITY_COMPONENT_MISMATCH")

  newer <- draft
  newer$header$schema_version <- 7L
  saveRDS(newer, paths[4], version = 3)
  err <- capture_calibration_error(mfrmr:::mfrmr_load_calibration(paths[4]))
  expect_identical(err$code, "SCHEMA_VERSION_UNSUPPORTED")
})

test_that("source readiness and optional lanes do not inherit core eligibility", {
  fit <- fixed_calibration_fit_fixture("RSM")$fit

  ineligible <- fit
  ineligible$readiness$fit$InferenceReady[1] <- FALSE
  draft <- mfrmr:::mfrmr_extract_calibration_draft(
    ineligible,
    calibration_id = "ineligible",
    source_fit_id = "fit-ineligible",
    created_at_utc = "2026-08-22T00:00:00Z"
  )
  expect_identical(draft$eligibility$source_readiness_status, "ineligible")
  expect_true("SOURCE_READINESS_INELIGIBLE" %in%
    mfrmr:::mfrmr_review_calibration(draft)$Code)
  err <- capture_calibration_error(
    mfrmr:::mfrmr_validate_calibration_draft(draft)
  )
  expect_identical(err$code, "SOURCE_READINESS_INELIGIBLE")

  jml <- fit
  jml$config$method <- "JML"
  err <- capture_calibration_error(
    mfrmr:::mfrmr_extract_calibration_draft(jml)
  )
  expect_identical(err$code, "SOURCE_READINESS_INELIGIBLE")

  for (par in list(head(fit$opt$par, -1L), c(fit$opt$par, 0))) {
    malformed <- fit; malformed$opt$par <- par
    err <- capture_calibration_error(extract_mfrm_calibration(malformed))
    expect_identical(err$code, "SOURCE_READINESS_INELIGIBLE")
    expect_match(conditionMessage(err), "Parameter vector length")
  }

  gpcm <- fit
  gpcm$config$model <- "unsupported"
  err <- capture_calibration_error(
    mfrmr:::mfrmr_extract_calibration_draft(gpcm)
  )
  expect_identical(err$code, "MODEL_FAMILY_UNSUPPORTED")

  latent <- fit
  latent$config$population_spec$active <- TRUE
  err <- capture_calibration_error(
    mfrmr:::mfrmr_extract_calibration_draft(latent)
  )
  expect_identical(err$code, "SCORING_BASIS_UNSUPPORTED")
})

test_that("portable scoring grid is independent of a one-point source fit grid", {
  toy <- load_mfrmr_data("example_core")
  fit <- suppressWarnings(fit_mfrm(
    toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 1,
    maxit = 100
  ))
  expect_true(mfrmr:::mfrm_inference_ready(fit))

  draft <- mfrmr:::mfrmr_extract_calibration_draft(
    fit,
    calibration_id = "source-one-scoring-thirty-one",
    source_fit_id = "fit-source-one",
    created_at_utc = "2026-08-24T00:00:00Z"
  )
  expect_identical(fit$config$estimation_control$quad_points, 1L)
  expect_identical(draft$scoring_basis$quadrature_order, 31L)
  expect_identical(
    draft$scoring_basis$scoring_algorithm, "quadrature_eap_v2"
  )

  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    draft, validated_at_utc = "2026-08-24T00:01:00Z"
  )
  frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-24T00:02:00Z"
  )
  rows <- data.frame(
    Person = c("LOW", "LOW", "HIGH", "HIGH"),
    Rater = c("R01", "R02", "R01", "R02"),
    Criterion = c("Accuracy", "Content", "Accuracy", "Content"),
    Score = c(1, 1, 4, 4),
    stringsAsFactors = FALSE
  )
  scored <- mfrmr:::mfrmr_score_calibration(frozen, rows)
  expect_true(all(scored$estimates$SD > 0))
  expect_gt(scored$estimates$Estimate[scored$estimates$Person == "HIGH"],
            scored$estimates$Estimate[scored$estimates$Person == "LOW"])

  err <- capture_calibration_error(
    mfrmr:::mfrmr_extract_calibration_draft(
      fit, scoring_quad_points = 1
    )
  )
  expect_identical(err$code, "QUADRATURE_ORDER_INVALID")
})

test_that("high scoring orders survive the full calibration lifecycle", {
  fixture <- fixed_calibration_fit_fixture("RSM")
  rows <- fixed_calibration_scoring_rows(fixture$data)
  expect_error(predict_mfrm_units(fixture$fit, rows, scoring_quad_points = 61L),
               "Posterior scoring integration did not pass")
  for (n in c(121L, 181L, 241L)) {
    draft <- mfrmr:::mfrmr_extract_calibration_draft(
      fixture$fit, scoring_quad_points = n
    )
    expect_identical(draft$scoring_basis$quadrature_order, n)
    expect_identical(draft$scoring_basis$quadrature_rule,
                     "gauss_hermite_standard_normal_recurrence_v2")
    expect_true(all(draft$scoring_basis$weights > 0))
    expect_equal(nrow(mfrmr:::mfrmr_review_calibration(draft)), 0L)
    frozen <- mfrmr:::mfrmr_freeze_calibration(
      mfrmr:::mfrmr_validate_calibration_draft(draft)
    )
    portable <- mfrmr:::mfrmr_score_calibration(frozen, rows)
    fitted <- predict_mfrm_units(fixture$fit, rows, scoring_quad_points = n)
    expect_true(all(fitted$estimates$ScoreIntegrationReady))
    expect_equal(portable$estimates$Person, fitted$estimates$Person)
    expect_equal(portable$estimates$Estimate, unname(fitted$estimates$Estimate),
                 tolerance = 1e-12)
    expect_equal(portable$estimates$SD, unname(fitted$estimates$SD), tolerance = 1e-12)
  }
  legacy <- draft
  legacy$scoring_basis$quadrature_rule <- "gauss_hermite_standard_normal_golub_welsch_v1"
  expect_false("QUADRATURE_RULE_INVALID" %in%
                 mfrmr:::mfrmr_review_calibration(legacy)$Code)
  legacy$scoring_basis$weights[1] <- 0
  expect_true("QUADRATURE_WEIGHTS_INVALID" %in%
                mfrmr:::mfrmr_review_calibration(legacy)$Code)
})

test_that("calibration print and summary expose scope without source Persons", {
  draft <- fixed_calibration_draft_fixture("RSM")$draft
  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    draft, validated_at_utc = "2026-08-22T00:01:00Z"
  )
  frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-22T00:02:00Z"
  )
  s <- summary(frozen)
  printed <- capture.output(print(frozen))
  summary_printed <- capture.output(print(s))

  expect_s3_class(s, "summary.mfrm_calibration")
  expect_identical(s$state, "frozen")
  expect_identical(
    s$support_profile, "rsm_mml_fixed_standard_normal_v1"
  )
  expect_true(any(grepl("<mfrm_calibration>", printed, fixed = TRUE)))
  expect_true(any(grepl("mfrmr Calibration Summary", summary_printed, fixed = TRUE)))
  expect_true(any(grepl("Scoring prior: standard normal", printed, fixed = TRUE)))
  expect_true(any(grepl("Scoring prior: standard normal", summary_printed, fixed = TRUE)))
  frozen_visibility <- NULL
  summary_visibility <- NULL
  capture.output(frozen_visibility <- withVisible(print(frozen)))
  capture.output(summary_visibility <- withVisible(print(s)))
  expect_false(frozen_visibility$visible)
  expect_false(summary_visibility$visible)
  expect_true(withVisible(summary(frozen))$visible)
  expect_error(
    mfrmr:::summary.mfrm_calibration(list()),
    "must be a portable calibration artifact",
    fixed = TRUE
  )
  expect_error(
    summary(structure(list(), class = "mfrm_calibration")),
    "not a complete portable calibration artifact",
    fixed = TRUE
  )
  malformed_artifact <- frozen
  malformed_artifact$parameters$coordinates <- NULL
  expect_error(
    summary(malformed_artifact),
    "not a complete portable calibration artifact",
    fixed = TRUE
  )
  expect_error(
    print(structure(list(), class = "summary.mfrm_calibration")),
    "not a complete portable calibration summary",
    fixed = TRUE
  )
  public_print <- paste(c(printed, summary_printed), collapse = "\n")
  expect_false(grepl("\\bLane\\b|core_|G1|OPT-[0-9]+|Support profile|fixed_standard_normal_v1", public_print,
                     perl = TRUE))
  expect_false(any(grepl(fixed_calibration_fit_fixture("RSM")$data$Person[1], printed, fixed = TRUE)))
})

test_that("portable artifacts use a user-facing support profile contract", {
  draft <- fixed_calibration_draft_fixture("PCM")$draft

  expect_identical(
    names(draft$eligibility),
    c(
      "support_profile_id", "source_readiness_contract",
      "source_readiness_status", "parameter_class_status"
    )
  )
  expect_false("lane_id" %in% names(draft$eligibility))
  expect_identical(
    draft$integrity$semantic_components$support_profile,
    "pcm_mml_fixed_standard_normal_v1"
  )
  expect_false(
    "eligibility_lane" %in% names(draft$integrity$semantic_components)
  )
  expect_identical(
    draft$provenance$created_by, "mfrmr::extract_mfrm_calibration"
  )

  malformed <- draft
  malformed$eligibility$support_profile_id <-
    "rsm_mml_fixed_standard_normal_v1"
  refusals <- mfrmr:::mfrmr_review_calibration(malformed)
  expect_true("SUPPORT_PROFILE_INVALID" %in% refusals$Code)
  expect_true(
    "eligibility.support_profile_id" %in% refusals$FieldPath
  )
})

test_that("fixed-calibration implementation helpers remain unexported", {
  exports <- getNamespaceExports("mfrmr")
  expect_false(any(c(
    "mfrmr_extract_calibration_draft", "mfrmr_validate_calibration_draft",
    "mfrmr_freeze_calibration", "mfrmr_save_calibration",
    "mfrmr_load_calibration"
  ) %in% exports))
})

test_that("artifact coordinates independently reproduce RSM and PCM fit scoring", {
  for (model in c("RSM", "PCM")) {
    fixture <- fixed_calibration_fit_fixture(model)
    fixture$frozen <- freeze_mfrm_calibration(validate_mfrm_calibration(
      mfrmr:::mfrmr_extract_calibration_draft(fixture$fit, scoring_quad_points = 121L)
    ))
    rows <- fixed_calibration_scoring_rows(fixture$data)
    rows$Weight <- rep(c(0.5, 1, 2), length.out = nrow(rows))

    artifact_score <- mfrmr:::mfrmr_score_calibration(
      fixture$frozen, rows, weight = "Weight", interval_level = 0.90
    )
    fit_score <- predict_mfrm_units(
      fixture$fit, rows, weight = "Weight", interval_level = 0.90,
      scoring_quad_points = 121L
    )
    expect_true(all(fit_score$estimates$ScoreIntegrationReady))
    artifact_estimates <- artifact_score$estimates[
      order(artifact_score$estimates$Person), , drop = FALSE
    ]
    fit_estimates <- as.data.frame(fit_score$estimates)[
      order(fit_score$estimates$Person), , drop = FALSE
    ]
    rownames(artifact_estimates) <- NULL
    rownames(fit_estimates) <- NULL

    expect_s3_class(artifact_score, "mfrm_calibration_score")
    expect_equal(
      artifact_estimates[c("Estimate", "SD", "Lower", "Upper")],
      fit_estimates[c("Estimate", "SD", "Lower", "Upper")],
      tolerance = 1e-14
    )
    expect_identical(
      artifact_estimates[c("Person", "Observations", "WeightedN")],
      fit_estimates[c("Person", "Observations", "WeightedN")]
    )
    expect_identical(artifact_score$settings$engine_identity,
                     "artifact_coordinates_v1")
    expect_identical(artifact_score$settings$lifecycle_state, "frozen")
  }

  scorer_source <- paste(deparse(body(mfrmr:::mfrmr_score_calibration)),
                         collapse = "\n")
  expect_false(grepl(
    "predict_mfrm_units|compute_person_posterior_summary|build_indices|expand_params|fit_mfrm",
    scorer_source
  ))
})

test_that("an independent direct oracle reproduces one-row RSM scoring", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  artifact <- fixture$frozen
  row <- fixed_calibration_scoring_rows(fixture$data, n_person = 1L)[1, , drop = FALSE]
  scored <- mfrmr:::mfrmr_score_calibration(
    artifact, row, interval_level = 0.80
  )

  coordinates <- artifact$parameters$coordinates
  base_eta <- 0
  for (facet in artifact$model$facet_names) {
    coordinate <- coordinates[
      coordinates$ParameterClass == "facet" &
        coordinates$OwnerFacet == facet &
        coordinates$Level == as.character(row[[facet]]),
      , drop = FALSE
    ]
    expect_identical(nrow(coordinate), 1L)
    base_eta <- base_eta + artifact$model$facet_signs[[facet]] * coordinate$Value
  }
  step_rows <- coordinates[coordinates$ParameterClass == "shared_step", , drop = FALSE]
  step_rows <- step_rows[order(as.integer(step_rows$Step)), , drop = FALSE]
  cumulative_step <- c(0, cumsum(step_rows$Value))
  score_internal <- artifact$response$score_map$InternalScore[
    match(row$Score, artifact$response$score_map$OriginalScore)
  ]
  score_k <- score_internal - artifact$response$rating_min
  category <- 0:(artifact$response$n_categories - 1L)
  log_likelihood <- vapply(artifact$scoring_basis$nodes, function(theta) {
    logit <- category * (theta + base_eta) - cumulative_step
    logit[score_k + 1L] - log(sum(exp(logit)))
  }, numeric(1))
  unnormalized <- artifact$scoring_basis$weights * exp(log_likelihood)
  posterior <- unnormalized / sum(unnormalized)
  oracle_estimate <- sum(artifact$scoring_basis$nodes * posterior)
  oracle_sd <- sqrt(sum(
    posterior * (artifact$scoring_basis$nodes - oracle_estimate)^2
  ))
  oracle_quantile <- function(probability) {
    artifact$scoring_basis$nodes[which(cumsum(posterior) >= probability)[1]]
  }

  expect_equal(scored$estimates$Estimate, oracle_estimate, tolerance = 1e-15)
  expect_equal(scored$estimates$SD, oracle_sd, tolerance = 1e-15)
  density <- function(theta) vapply(theta, function(value) {
    logits <- category * (value + base_eta) - cumulative_step
    p <- exp(logits - max(logits))
    p[score_k+1L] / sum(p) * dnorm(value)
  }, 0)
  mass <- integrate(density, -12, 12, rel.tol=1e-11)$value
  cdf <- function(value) integrate(density, -12, value, rel.tol=1e-11)$value / mass
  expect_equal(cdf(scored$estimates$Lower), 0.10, tolerance=1e-8)
  expect_equal(cdf(scored$estimates$Upper), 0.90, tolerance=1e-8)

  # Historical algorithm identities still reproduce their original grid bounds.
  legacy <- fixture$draft
  legacy$scoring_basis$scoring_algorithm <- "quadrature_eap_v1"
  legacy$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(legacy)
  legacy <- mfrmr:::mfrmr_freeze_calibration(mfrmr:::mfrmr_validate_calibration_draft(legacy))
  old <- mfrmr:::mfrmr_score_calibration(legacy,row,interval_level=0.80)
  expect_equal(old$estimates$Lower,oracle_quantile(0.10),tolerance=0)
  expect_equal(old$estimates$Upper,oracle_quantile(0.90),tolerance=0)
  expect_match(paste(old$notes,collapse=" "),"legacy v1")
})

test_that("interaction artifacts score from their complete stored cell matrix", {
  fixture <- fixed_calibration_interaction_fixture()
  validated <- mfrmr:::mfrmr_validate_calibration_draft(
    mfrmr:::mfrmr_extract_calibration_draft(fixture$fit, scoring_quad_points = 181L),
    validated_at_utc = "2026-08-22T00:01:00Z"
  )
  frozen <- mfrmr:::mfrmr_freeze_calibration(
    validated, frozen_at_utc = "2026-08-22T00:02:00Z"
  )
  persons <- unique(as.character(fixture$fit$prep$data$Person))[1:2]
  rows <- fixture$fit$prep$data[
    as.character(fixture$fit$prep$data$Person) %in% persons,
    c("Person", "Rater", "Task", "Criterion", "Score"), drop = FALSE
  ]
  rows$Person <- paste0("INTERACTION", match(as.character(rows$Person), persons))
  rows$Score <- frozen$response$score_map$OriginalScore[
    match(as.integer(as.character(rows$Score)), frozen$response$score_map$InternalScore)
  ]
  rownames(rows) <- NULL

  artifact_score <- mfrmr:::mfrmr_score_calibration(frozen, rows)
  fit_score <- predict_mfrm_units(fixture$fit, rows, scoring_quad_points = 181L)
  expect_true(all(fit_score$estimates$ScoreIntegrationReady))
  artifact_estimates <- artifact_score$estimates[
    order(artifact_score$estimates$Person), , drop = FALSE
  ]
  fit_estimates <- as.data.frame(fit_score$estimates)[
    order(fit_score$estimates$Person), , drop = FALSE
  ]
  expect_equal(
    artifact_estimates[c("Estimate", "SD", "Lower", "Upper")],
    fit_estimates[c("Estimate", "SD", "Lower", "Upper")],
    tolerance = 1e-14
  )
})

test_that("public RSM and PCM interaction calibrations retain reviewed fit scoring", {
  for (model in c("RSM", "PCM")) {
    fixture <- fixed_calibration_interaction_fixture(model)
    review <- suppressWarnings(mml_quadrature_sensitivity(
      fixture$fit, fixture$data, quad_points = c(5L, 7L), theta_points = 41L
    ))
    fit <- review$fits$q7
    draft <- extract_mfrm_calibration(fit, quadrature_review = review,
                                     scoring_quad_points = 241L)
    artifact <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    path <- tempfile(fileext = ".rds")
    save_mfrm_calibration(artifact, path)
    restored <- load_mfrm_calibration(path)
    unlink(path)
    expect_identical(restored, artifact)
    persons <- unique(fixture$data$Person)[1:2]
    rows <- fixture$data[fixture$data$Person %in% persons, , drop = FALSE]
    rows$Person <- paste0("NEW", match(rows$Person, persons))
    score <- score_mfrm_calibration(restored, rows)
    fitted <- predict_mfrm_units(fit, rows, scoring_quad_points = 241L)
    expect_true(all(fitted$estimates$ScoreIntegrationReady))
    fields <- c("Person", "Estimate", "SD", "Lower", "Upper")
    expect_equal(score$estimates[fields], as.data.frame(fitted$estimates)[fields],
                 tolerance = 1e-12)
    expect_identical(nrow(restored$model$interactions), 1L)
    expect_identical(restored$scoring_basis$quadrature_order, 241L)
    expect_identical(as.integer(fit$config$estimation_control$quad_points), 7L)
  }
})

test_that("artifact scoring is row-order, chunk-order, score-map, RNG, and option invariant", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  artifact <- fixture$frozen
  rows <- fixed_calibration_scoring_rows(fixture$data)
  baseline <- mfrmr:::mfrmr_score_calibration(artifact, rows)

  reversed <- mfrmr:::mfrmr_score_calibration(
    artifact, rows[rev(seq_len(nrow(rows))), , drop = FALSE]
  )
  baseline_sorted <- baseline$estimates[order(baseline$estimates$Person), , drop = FALSE]
  reversed_sorted <- reversed$estimates[order(reversed$estimates$Person), , drop = FALSE]
  rownames(baseline_sorted) <- NULL
  rownames(reversed_sorted) <- NULL
  expect_equal(reversed_sorted, baseline_sorted, tolerance = 1e-15)

  chunks <- split(rows, rows$Person)
  chunk_estimates <- do.call(rbind, lapply(chunks, function(chunk) {
    mfrmr:::mfrmr_score_calibration(artifact, chunk)$estimates
  }))
  chunk_estimates <- chunk_estimates[order(chunk_estimates$Person), , drop = FALSE]
  rownames(chunk_estimates) <- NULL
  expect_equal(chunk_estimates, baseline_sorted, tolerance = 1e-15)

  remapped <- artifact
  replacement <- seq(10, by = 10, length.out = nrow(remapped$response$score_map))
  remapped$response$score_map$OriginalScore <- replacement
  remapped$integrity$semantic_components <-
    mfrmr:::mfrmr_calibration_semantic_components(remapped)
  remapped_rows <- rows
  remapped_rows$Score <- replacement[
    match(rows$Score, artifact$response$score_map$OriginalScore)
  ]
  remapped_score <- mfrmr:::mfrmr_score_calibration(remapped, remapped_rows)
  remapped_sorted <- remapped_score$estimates[
    order(remapped_score$estimates$Person), , drop = FALSE
  ]
  rownames(remapped_sorted) <- NULL
  expect_equal(remapped_sorted, baseline_sorted, tolerance = 1e-15)

  renamed <- rows
  names(renamed) <- c("Subject", "Judge", "Domain", "Rating")
  renamed$Judge <- factor(
    renamed$Judge, levels = rev(unique(as.character(renamed$Judge)))
  )
  renamed$Domain <- factor(
    renamed$Domain, levels = rev(unique(as.character(renamed$Domain)))
  )
  mapped_score <- mfrmr:::mfrmr_score_calibration(
    artifact,
    renamed,
    person = "Subject",
    facets = c(Rater = "Judge", Criterion = "Domain"),
    score = "Rating"
  )
  mapped_sorted <- mapped_score$estimates[
    order(mapped_score$estimates$Person), , drop = FALSE
  ]
  rownames(mapped_sorted) <- NULL
  expect_equal(mapped_sorted, baseline_sorted, tolerance = 1e-15)

  set.seed(240824)
  rng_before <- .Random.seed
  old_options <- options(contrasts = c("contr.sum", "contr.poly"), digits = 4)
  on.exit(options(old_options), add = TRUE)
  option_score <- mfrmr:::mfrmr_score_calibration(artifact, rows)
  expect_identical(.Random.seed, rng_before)
  option_sorted <- option_score$estimates[
    order(option_score$estimates$Person), , drop = FALSE
  ]
  rownames(option_sorted) <- NULL
  expect_equal(option_sorted, baseline_sorted, tolerance = 1e-15)
})

test_that("artifact scoring fails closed with stable input and operational codes", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  artifact <- fixture$frozen
  rows <- fixed_calibration_scoring_rows(fixture$data)
  artifact_before <- artifact
  rows_before <- rows

  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(fixture$draft, rows)
    )$code,
    "LIFECYCLE_NOT_FROZEN"
  )
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, rows[, -1, drop = FALSE])
    )$code,
    "SCORING_COLUMN_MISSING"
  )
  unknown_level <- rows
  unknown_level$Rater[1] <- "UNSEEN"
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, unknown_level)
    )$code,
    "SCORING_FACET_LEVEL_UNKNOWN"
  )
  unknown_score <- rows
  unknown_score$Score[1] <- 999
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, unknown_score)
    )$code,
    "SCORING_SCORE_UNKNOWN"
  )
  missing_value <- rows
  missing_value$Criterion[1] <- NA_character_
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, missing_value)
    )$code,
    "SCORING_VALUE_MISSING"
  )
  bad_weight <- rows
  bad_weight$Weight <- 1
  bad_weight$Weight[1] <- 0
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, bad_weight, weight = "Weight")
    )$code,
    "SCORING_WEIGHT_INVALID"
  )
  bad_weight$Weight[1] <- Inf
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, bad_weight, weight = "Weight")
    )$code,
    "SCORING_WEIGHT_INVALID"
  )
  duplicate <- rbind(rows, rows[1, , drop = FALSE])
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, duplicate)
    )$code,
    "SCORING_EVENT_DUPLICATE"
  )
  malformed <- artifact
  facet_row <- which(malformed$parameters$coordinates$ParameterClass == "facet")[1]
  malformed$parameters$coordinates$Level[facet_row] <-
    malformed$parameters$coordinates$Level[facet_row + 1L]
  malformed$integrity$semantic_components <-
    mfrmr:::mfrmr_calibration_semantic_components(malformed)
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(malformed, rows)
    )$code,
    "PARAMETER_COORDINATE_INVALID"
  )
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(artifact, rows, interval_level = 1)
    )$code,
    "SCORING_INTERVAL_INVALID"
  )
  expect_identical(artifact, artifact_before)
  expect_identical(rows, rows_before)
})

test_that("operational scoring applies an explicit missing-response policy", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  rows <- fixed_calibration_scoring_rows(fixture$data, n_person = 2L)
  rows$Score[rows$Person == "NEW1"] <- NA_real_
  rows$Score[which(rows$Person == "NEW2")[1L]] <- NA_real_

  default_error <- capture_calibration_error(
    mfrmr:::mfrmr_score_calibration(fixture$frozen, rows)
  )
  expect_identical(default_error$code, "SCORING_SCORE_INVALID")
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(
        fixture$frozen, rows, missing_response = "guess"
      )
    )$code,
    "SCORING_MISSING_POLICY_INVALID"
  )

  scored <- mfrmr:::mfrmr_score_calibration(
    fixture$frozen, rows, missing_response = "omit"
  )
  expect_identical(scored$settings$missing_response_policy, "omit")
  expect_identical(nrow(scored$row_dispositions), nrow(rows))
  expect_identical(
    scored$row_dispositions$Disposition[is.na(rows$Score)],
    rep("omitted", sum(is.na(rows$Score)))
  )
  expect_identical(
    scored$row_dispositions$ReasonCode[is.na(rows$Score)],
    rep("RESPONSE_MISSING_OMITTED", sum(is.na(rows$Score)))
  )
  expect_identical(scored$row_review$InputRows, as.integer(nrow(rows)))
  expect_identical(scored$row_review$ScoredRows, as.integer(sum(!is.na(rows$Score))))
  expect_identical(scored$row_review$OmittedRows, as.integer(sum(is.na(rows$Score))))

  new1 <- scored$person_dispositions[
    scored$person_dispositions$Person == "NEW1", , drop = FALSE
  ]
  new2 <- scored$person_dispositions[
    scored$person_dispositions$Person == "NEW2", , drop = FALSE
  ]
  expect_identical(new1$Disposition, "not_scored")
  expect_identical(new1$ReasonCodes, "ZERO_VALID_RESPONSES")
  expect_identical(new1$ValidResponses, 0L)
  expect_true(new1$OmittedResponses > 0L)
  expect_false("NEW1" %in% scored$estimates$Person)
  expect_true("NEW2" %in% scored$estimates$Person)
  expect_identical(new2$OmittedResponses, 1L)
  expect_match(new2$ReasonCodes, "RESPONSES_MISSING_OMITTED", fixed = TRUE)
  expect_identical(
    new2$AdministrationCompleteness,
    "partial_supplied_missing_response"
  )
})

test_that("event identifiers distinguish intentional repeated observations", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  rows <- fixed_calibration_scoring_rows(fixture$data, n_person = 1L)[1, , drop = FALSE]
  repeated <- rbind(rows, rows)
  repeated$EventId <- c("occasion-1", "occasion-2")

  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(fixture$frozen, repeated)
    )$code,
    "SCORING_EVENT_DUPLICATE"
  )
  scored <- mfrmr:::mfrmr_score_calibration(
    fixture$frozen, repeated, event_id = "EventId"
  )
  expect_identical(scored$estimates$Observations, 2L)
  expect_identical(scored$settings$event_id_column, "EventId")
  expect_identical(scored$row_dispositions$EventId, repeated$EventId)

  repeated$EventId[2L] <- repeated$EventId[1L]
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(
        fixture$frozen, repeated, event_id = "EventId"
      )
    )$code,
    "SCORING_EVENT_DUPLICATE"
  )
})

test_that("endpoint and sparse EAP scores carry explicit conditional dispositions", {
  fixture <- fixed_calibration_frozen_fixture("RSM")
  artifact <- fixture$frozen
  rows <- fixed_calibration_scoring_rows(fixture$data, n_person = 2L)
  rows <- rows[!duplicated(rows$Person), , drop = FALSE]
  rows$Score <- c(
    min(artifact$response$score_map$OriginalScore),
    max(artifact$response$score_map$OriginalScore)
  )

  scored <- mfrmr:::mfrmr_score_calibration(artifact, rows)
  dispositions <- scored$person_dispositions[
    match(c("NEW1", "NEW2"), scored$person_dispositions$Person), , drop = FALSE
  ]
  expect_identical(
    dispositions$EndpointStatus,
    c("all_lower_endpoint", "all_upper_endpoint")
  )
  expect_true(all(dispositions$VerySparsePattern))
  expect_true(all(dispositions$Disposition == "scored_review"))
  expect_true(all(grepl("VERY_SPARSE_RESPONSE_PATTERN", dispositions$ReasonCodes)))
  expect_true(all(grepl("ALL_RESPONSES_.*_ENDPOINT", dispositions$ReasonCodes)))
  expect_true(all(is.finite(dispositions$QuadratureEdgeMass)))
  expect_true(all(dispositions$QuadratureEdgeMass >= 0 &
                    dispositions$QuadratureEdgeMass <= 1))
  expect_identical(
    grepl("QUADRATURE_EDGE_MASS_REVIEW", dispositions$ReasonCodes,
          fixed = TRUE),
    dispositions$QuadratureEdgeMass >= dispositions$QuadratureEdgeThreshold
  )
  expect_identical(
    dispositions$PriorSensitivityStatus,
    rep("not_evaluated_fixed_basis", 2L)
  )
  expect_identical(
    dispositions$EstimateBasis,
    rep("posterior_eap_fixed_calibration", 2L)
  )
  expect_identical(
    dispositions$UncertaintyBasis,
    rep("conditional_on_frozen_point_calibration", 2L)
  )
  expect_true(all(is.finite(scored$estimates$Estimate)))
  expect_true(all(is.finite(scored$estimates$SD)))
  expect_false(any(grepl("JML", scored$estimates$EstimateBasis, fixed = TRUE)))
})

test_that("every operational score carries calibration and basis identities", {
  fixture <- fixed_calibration_frozen_fixture("PCM")
  rows <- fixed_calibration_scoring_rows(fixture$data, n_person = 2L)
  scored <- mfrmr:::mfrmr_score_calibration(fixture$frozen, rows)

  expect_true(all(c(
    "row_dispositions", "person_dispositions", "settings"
  ) %in% names(scored)))
  expect_identical(
    unique(scored$estimates$CalibrationId),
    fixture$frozen$header$calibration_id
  )
  expect_identical(
    unique(scored$estimates$SchemaVersion),
    fixture$frozen$header$schema_version
  )
  expect_identical(
    scored$settings$semantic_components,
    fixture$frozen$integrity$semantic_components
  )
  expect_identical(
    unique(scored$estimates$ScoringBasis),
    fixture$frozen$scoring_basis$type
  )
  expect_identical(
    scored$settings$scoring_contract_id,
    "mfrmr_operational_scoring_v1"
  )
  expect_identical(
    scored$settings$prior_identity,
    list(type = "fixed_standard_normal", mean = 0, sd = 1)
  )
  expect_identical(scored$settings$score_map, fixture$frozen$response$score_map)
  expect_identical(
    scored$settings$source_readiness_contract,
    fixture$frozen$eligibility$source_readiness_contract
  )
  expect_true(any(grepl(
    "conditional on the frozen point calibration", scored$notes, fixed = TRUE
  )))
})

test_that("supersession and retirement create immutable terminal lineage records", {
  artifact <- fixed_calibration_frozen_fixture("RSM")$frozen
  artifact_before <- artifact
  superseded <- mfrmr:::mfrmr_supersede_calibration(
    artifact,
    record_id = "calibration-rsm-superseded-record",
    superseded_at_utc = "2026-08-22T00:03:00Z"
  )
  retired <- mfrmr:::mfrmr_retire_calibration(
    artifact,
    record_id = "calibration-rsm-retired-record",
    retired_at_utc = "2026-08-22T00:04:00Z"
  )

  expect_identical(artifact, artifact_before)
  expect_identical(superseded$lifecycle$state, "superseded")
  expect_identical(retired$lifecycle$state, "retired")
  expect_identical(superseded$lifecycle$revision, 4L)
  expect_identical(retired$lifecycle$revision, 4L)
  expect_identical(
    superseded$lifecycle$events$Operation,
    c("extract_draft", "validate", "freeze", "supersede")
  )
  expect_identical(
    retired$lifecycle$events$Operation,
    c("extract_draft", "validate", "freeze", "retire")
  )
  expect_identical(
    superseded$provenance$parent_calibration_id,
    artifact$header$calibration_id
  )
  expect_identical(
    retired$provenance$parent_calibration_id,
    artifact$header$calibration_id
  )
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(superseded)), 0L)
  expect_identical(nrow(mfrmr:::mfrmr_review_calibration(retired)), 0L)

  paths <- vapply(seq_len(2), function(i) tempfile(fileext = ".rds"), character(1))
  on.exit(unlink(paths), add = TRUE)
  mfrmr:::mfrmr_save_calibration(superseded, paths[1])
  mfrmr:::mfrmr_save_calibration(retired, paths[2])
  expect_identical(mfrmr:::mfrmr_load_calibration(paths[1]), superseded)
  expect_identical(mfrmr:::mfrmr_load_calibration(paths[2]), retired)

  rows <- fixed_calibration_scoring_rows(
    fixed_calibration_fit_fixture("RSM")$data, n_person = 1L
  )
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_score_calibration(superseded, rows)
    )$code,
    "LIFECYCLE_NOT_FROZEN"
  )
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_supersede_calibration(
        artifact, artifact$header$calibration_id
      )
    )$code,
    "PROVENANCE_PARENT_INVALID"
  )
  expect_identical(
    capture_calibration_error(
      mfrmr:::mfrmr_retire_calibration(
        superseded, "second-terminal-record"
      )
    )$code,
    "LIFECYCLE_TRANSITION_INVALID"
  )
})

test_that("provenance and event-chain mutations fail semantic review", {
  artifact <- fixed_calibration_frozen_fixture("RSM")$frozen

  bad_validation_time <- artifact
  bad_validation_time$validation$validated_at_utc <- "not-a-time"
  expect_true("VALIDATION_TIME_INVALID" %in%
    mfrmr:::mfrmr_review_calibration(bad_validation_time)$Code)

  bad_package_version <- artifact
  bad_package_version$provenance$source_package_version <- "not a version"
  expect_true("PROVENANCE_PACKAGE_VERSION_INVALID" %in%
    mfrmr:::mfrmr_review_calibration(bad_package_version)$Code)

  bad_creator <- artifact
  bad_creator$provenance$created_by <- "unregistered-constructor"
  expect_true("PROVENANCE_CREATOR_INVALID" %in%
    mfrmr:::mfrmr_review_calibration(bad_creator)$Code)

  bad_parent <- artifact
  bad_parent$provenance$parent_calibration_id <- "invented-parent"
  expect_true("PROVENANCE_PARENT_INVALID" %in%
    mfrmr:::mfrmr_review_calibration(bad_parent)$Code)

  bad_event <- artifact
  bad_event$lifecycle$events$Operation[3] <- "retire"
  expect_true("LIFECYCLE_EVENT_CHAIN_INVALID" %in%
    mfrmr:::mfrmr_review_calibration(bad_event)$Code)
})

test_that("artifact scorer and lifecycle functions remain unexported", {
  exports <- getNamespaceExports("mfrmr")
  expect_false(any(c(
    "mfrmr_score_calibration", "mfrmr_supersede_calibration",
    "mfrmr_retire_calibration"
  ) %in% exports))
})
