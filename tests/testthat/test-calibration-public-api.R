calibration_public_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    data <- load_mfrmr_data("example_core")
    persons <- unique(as.character(data$Person))[seq_len(18L)]
    training <- data[as.character(data$Person) %in% persons, , drop = FALSE]
    fit <- suppressWarnings(fit_mfrm(
      training,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      method = "MML",
      model = "RSM",
      quad_points = 5,
      maxit = 20
    ))
    quadrature_review <- suppressWarnings(mml_quadrature_sensitivity(
      fit, training, quad_points = c(5L, 7L), theta_points = 41L
    ))
    fit <- quadrature_review$fits$q7
    draft <- extract_mfrm_calibration(
      fit,
      calibration_id = "public-api-rsm",
      source_fit_id = "public-api-source-fit",
      created_at_utc = "2026-08-26T00:00:00Z",
      quadrature_review = quadrature_review
    )
    validated <- validate_mfrm_calibration(
      draft, validated_at_utc = "2026-08-26T00:01:00Z"
    )
    frozen <- freeze_mfrm_calibration(
      validated, frozen_at_utc = "2026-08-26T00:02:00Z"
    )
    rows <- training[
      as.character(training$Person) %in% persons[seq_len(2L)],
      c("Person", "Rater", "Criterion", "Score"),
      drop = FALSE
    ]
    person_map <- stats::setNames(c("NEW_LOW", "NEW_HIGH"), persons[seq_len(2L)])
    rows$Person <- unname(person_map[as.character(rows$Person)])
    rows$Score[rows$Person == "NEW_LOW"] <- min(training$Score)
    rows$Score[rows$Person == "NEW_HIGH"] <- max(training$Score)
    rownames(rows) <- NULL
    cache <<- list(
      data = training, fit = fit, quadrature_review = quadrature_review,
      draft = draft, validated = validated, frozen = frozen, rows = rows
    )
    cache
  }
})

test_that("portable calibration exports only the reviewed public workflow", {
  expected <- c(
    "mfrm_calibration_capabilities", "extract_mfrm_calibration",
    "review_mfrm_calibration", "validate_mfrm_calibration",
    "freeze_mfrm_calibration", "supersede_mfrm_calibration",
    "retire_mfrm_calibration", "save_mfrm_calibration",
    "load_mfrm_calibration", "score_mfrm_calibration"
  )
  exports <- getNamespaceExports("mfrmr")
  expect_true(all(expected %in% exports))
  expect_false(any(c(
    "mfrmr_extract_calibration_draft", "mfrmr_review_calibration",
    "mfrmr_validate_calibration_draft", "mfrmr_freeze_calibration",
    "mfrmr_score_calibration"
  ) %in% exports))
})

test_that("public capability matrix states the bounded artifact envelope", {
  capabilities <- mfrm_calibration_capabilities()
  expect_identical(nrow(capabilities), 8L)
  expect_identical(
    names(capabilities),
    c(
      "Model", "Estimator", "ScoringBasis", "PortableCalibration",
      "AnchorSupport", "InteractionSupport", "ExistingAlternative",
      "Limitation"
    )
  )
  available <- capabilities$PortableCalibration == "available"
  expect_identical(which(available), c(1L, 2L, 4L, 5L, 6L, 7L, 8L))
  expect_identical(capabilities$Model[available], c("RSM", "PCM", "GPCM", "RSM/PCM", "GPCM", "GPCM", "GPCM (two families)"))
  expect_identical(capabilities$Estimator[available], c("MML", "MML", "MML", "JML", "JML", "Corrected JML", "MML"))
  expect_true(all(
    capabilities$ScoringBasis[1:2] == "fixed standard normal"
  ))
  expect_true(all(nzchar(capabilities$ExistingAlternative[!available])))
  expect_true(all(nzchar(capabilities$Limitation)))

  guide <- mfrmr_output_guide("calibration")
  expect_identical(nrow(guide), 7L)
  decision_route <- guide[guide$ObjectRole == "conditional ability decisions",]
  expect_equal(nrow(decision_route),1L)
  expect_match(decision_route$NextStep,"apa_table",fixed=TRUE)
  expect_match(decision_route$DecisionBoundary,"no automatic decision rule",fixed=TRUE)
  expect_true(all(guide$Scope == "calibration"))
  expect_true(all(guide$RecommendedEntry))
  expect_true(any(grepl(
    "mfrm_calibration_capabilities()", guide$MainFunction, fixed = TRUE
  )))
  expect_true(any(grepl(
    "score_mfrm_calibration()", guide$MainFunction, fixed = TRUE
  )))
  expect_true(any(grepl(
    "plot(scores, type = \"interval\")", guide$MainFunction, fixed = TRUE
  )))
  expect_true(any(grepl(
    "GPCM currently requires unit weights and no anchors or interactions", guide$Notes, fixed = TRUE
  )))
  expect_true(any(grepl(
    "excludes calibration-parameter uncertainty",
    guide$DecisionBoundary,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "does not establish source-fit quality",
    guide$DecisionBoundary,
    fixed = TRUE
  )))
})

test_that("public lifecycle is explicit and artifact scoring performs no refit", {
  fixture <- calibration_public_fixture()
  expect_s3_class(fixture$draft, "mfrm_calibration")
  expect_identical(fixture$draft$lifecycle$state, "draft")
  expect_identical(nrow(review_mfrm_calibration(fixture$draft)), 0L)
  expect_identical(fixture$validated$lifecycle$state, "validated")
  expect_identical(fixture$frozen$lifecycle$state, "frozen")

  scored <- score_mfrm_calibration(
    fixture$frozen, fixture$rows, interval_level = 0.90
  )
  expect_s3_class(scored, "mfrm_calibration_score")
  expect_identical(scored$settings$engine_identity, "artifact_coordinates_v1")
  expect_identical(scored$settings$lifecycle_state, "frozen")
  expect_true(all(scored$estimates$SD > 0))
  expect_gt(
    scored$estimates$Estimate[scored$estimates$Person == "NEW_HIGH"],
    scored$estimates$Estimate[scored$estimates$Person == "NEW_LOW"]
  )

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saved <- save_mfrm_calibration(fixture$frozen, path)
  expect_identical(normalizePath(path), saved)
  restored <- load_mfrm_calibration(path)
  expect_identical(restored, fixture$frozen)
  expect_identical(
    score_mfrm_calibration(restored, fixture$rows),
    score_mfrm_calibration(fixture$frozen, fixture$rows)
  )

  retired <- retire_mfrm_calibration(
    fixture$frozen,
    record_id = "public-api-rsm-retired",
    retired_at_utc = "2026-08-26T00:03:00Z"
  )
  expect_identical(retired$lifecycle$state, "retired")
  error <- tryCatch(
    score_mfrm_calibration(retired, fixture$rows),
    mfrm_calibration_error = identity
  )
  expect_s3_class(error, "mfrm_calibration_error")
  expect_identical(error$code, "LIFECYCLE_NOT_FROZEN")
})

test_that("portable score methods foreground review and conditional uncertainty", {
  fixture <- calibration_public_fixture()
  scored <- score_mfrm_calibration(fixture$frozen, fixture$rows)

  printed <- capture.output(print(scored))
  summarized <- summary(scored, digits = 2)
  summary_printed <- capture.output(print(summarized))

  expect_s3_class(summarized, "summary.mfrm_calibration_score")
  expect_identical(
    names(summarized),
    c("overview", "estimates", "review", "row_review", "settings", "notes", "digits")
  )
  expect_identical(summarized$overview$Persons, 2L)
  expect_identical(summarized$overview$Scored, 2L)
  expect_identical(summarized$overview$Review, 2L)
  expect_identical(summarized$overview$NotScored, 0L)
  expect_true(all(c(
    "Person", "Estimate", "SD", "Lower", "Upper", "Observations",
    "Disposition", "ReasonCodes"
  ) %in% names(summarized$estimates)))
  semantic_columns <- c(
    "Person", "Disposition", "EstimateBasis", "UncertaintyBasis",
    "CalibrationId", "SchemaVersion", "ScoringBasis", "ScoringAlgorithm",
    "IntervalLevel"
  )
  expect_identical(
    summarized$estimates[semantic_columns], scored$estimates[semantic_columns]
  )
  score_columns <- c("Estimate", "SD", "Lower", "Upper")
  expect_identical(
    summarized$estimates[score_columns],
    round(scored$estimates[score_columns], digits = 2)
  )
  csv_path <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_path), add = TRUE)
  utils::write.csv(summarized$estimates, csv_path, row.names = FALSE)
  exported <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  expect_identical(
    exported[semantic_columns], scored$estimates[semantic_columns]
  )
  expect_true(all(summarized$review$Disposition == "scored_review"))
  expect_true(any(grepl("<mfrm_calibration_score>", printed, fixed = TRUE)))
  expect_true(any(grepl(
    "2 scored (2 requiring review); 0 not scored", printed, fixed = TRUE
  )))
  expect_true(any(grepl(
    "2 scored (2 requiring review); 0 not scored", summary_printed,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "calibration-parameter uncertainty excluded", printed, fixed = TRUE
  )))
  expect_true(any(grepl(
    "Interval interpretation", summary_printed, fixed = TRUE
  )))
  public_output <- paste(c(printed, summary_printed), collapse = "\n")
  expect_false(grepl("CORE-[0-9]|G[0-6] exit|internal", public_output, perl = TRUE))
  expect_false(grepl("scored_review|ALL_RESPONSES_", public_output))
  expect_match(public_output, "all responses at the minimum score", fixed = TRUE)

  old_width <- options(width = 80L)
  on.exit(options(old_width), add = TRUE)
  direct_visibility <- NULL
  direct_lines <- capture.output(
    direct_visibility <- withVisible(print(scored))
  )
  summary_visibility <- NULL
  summary_lines <- capture.output(
    summary_visibility <- withVisible(print(summarized))
  )
  expect_false(direct_visibility$visible)
  expect_identical(direct_visibility$value, scored)
  expect_false(summary_visibility$visible)
  expect_identical(summary_visibility$value, summarized)
  expect_lte(max(nchar(c(direct_lines, summary_lines))), 80L)
  expect_true(withVisible(summary(scored))$visible)

  long_id <- scored
  long_id$settings$calibration_id <- paste(rep("LONG-CALIBRATION-ID", 12L),
                                           collapse = "-")
  long_lines <- capture.output({
    print(long_id)
    print(summary(long_id))
  })
  expect_lte(max(nchar(long_lines)), 80L)
})

test_that("portable score plots preserve selection and non-scored dispositions", {
  fixture <- calibration_public_fixture()
  rows <- fixture$rows
  rows$Score[rows$Person == "NEW_LOW"] <- NA_real_
  scored <- score_mfrm_calibration(
    fixture$frozen, rows, missing_response = "omit"
  )

  summarized <- summary(scored)
  expect_identical(summarized$overview$NotScored, 1L)
  expect_true("NEW_LOW" %in% summarized$review$Person)

  interval <- plot(scored, type = "interval", draw = FALSE)
  precision <- plot(scored, type = "precision", draw = FALSE)
  edge_mass <- plot(scored, type = "edge_mass", draw = FALSE)
  expect_s3_class(interval, "mfrm_plot_data")
  expect_identical(interval$name, "calibration_score_interval")
  expect_identical(precision$name, "calibration_score_precision")
  expect_identical(edge_mass$name, "calibration_score_edge_mass")
  expect_true(all(c(
    "data", "selection_summary", "unplotted_dispositions",
    "interpretation", "display_encoding", "settings", "legend",
    "reference_lines"
  ) %in% names(interval$data)))
  expect_identical(
    interval$data$display_encoding$Status, c("Scored", "Review")
  )
  expect_identical(interval$data$display_encoding$Shape, c(16L, 17L))
  expect_true(all(
    interval$data$legend$aesthetic == "colour and shape"
  ))
  expect_true("NEW_LOW" %in% interval$data$unplotted_dispositions$Person)
  expect_identical(interval$data$selection_summary$NotScoredPersons, 1L)
  expect_identical(
    interval$data$subtitle,
    "1 of 1 scored Persons shown; 1 review disposition shown; 1 not scored"
  )
  expect_true(all(c(
    "ReviewFlag", "QuadratureEdgeMass", "QuadratureEdgeThreshold",
    "DisplayOrder", "LabelPosition", "LabelX", "LabelY"
  ) %in% names(interval$data$data)))
  expect_match(
    interval$data$interpretation$Guidance[3L],
    "calibration-parameter uncertainty excluded",
    fixed = TRUE
  )
  expect_match(
    interval$data$caption,
    "calibration-parameter uncertainty excluded",
    fixed = TRUE
  )

  mixed_rows <- fixture$rows[fixture$rows$Person == "NEW_LOW", , drop = FALSE]
  mixed_rows$Person <- "NEW_MIXED"
  mixed_rows$Score <- rep(c(2, 3), length.out = nrow(mixed_rows))
  prioritized <- score_mfrm_calibration(
    fixture$frozen, rbind(fixture$rows, mixed_rows)
  )
  precision_prioritized <- plot(
    prioritized, type = "precision", draw = FALSE
  )
  edge_prioritized <- plot(
    prioritized, type = "edge_mass", draw = FALSE
  )
  expect_true(all(
    precision_prioritized$data$data$LabelPosition[
      precision_prioritized$data$data$ReviewFlag
    ] == 3L
  ))
  precision_review <- precision_prioritized$data$data[
    precision_prioritized$data$data$ReviewFlag, , drop = FALSE
  ]
  expect_true(all(is.finite(precision_review$LabelX)))
  expect_true(all(is.finite(precision_review$LabelY)))
  expect_true(all(precision_review$LabelY > precision_review$SD))
  expect_identical(length(unique(precision_review$LabelY)),
                   nrow(precision_review))
  edge_review <- edge_prioritized$data$data[
    edge_prioritized$data$data$ReviewFlag, , drop = FALSE
  ]
  expect_true(all(
    edge_review$LabelPosition[edge_review$Estimate < 0] == 4L
  ))
  expect_true(all(
    edge_review$LabelPosition[edge_review$Estimate > 0] == 2L
  ))
  truncated <- plot(prioritized, top_n = 1L, draw = FALSE)
  expect_true(truncated$data$data$ReviewFlag)
  expect_identical(truncated$data$selection_summary$PlottedPersons, 1L)
  expect_identical(truncated$data$selection_summary$OmittedFromPlot, 2L)
  expect_true(truncated$data$selection_summary$ReviewRowsPrioritized)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(scored, type = "interval", draw = TRUE))
  expect_no_error(plot(scored, type = "precision", draw = TRUE))
  expect_no_error(plot(scored, type = "edge_mass", draw = TRUE))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplots <- lapply(
      list(interval, precision, edge_mass), as_ggplot
    )
    for (gg in ggplots) {
      expect_no_error(ggplot2::ggplot_build(gg))
      expect_identical(
        levels(gg$data$ReviewStatus), c("Scored", "Review")
      )
      shape_scales <- vapply(
        gg$scales$scales,
        function(scale) "shape" %in% scale$aesthetics,
        logical(1)
      )
      expect_true(any(shape_scales))
    }
  }
  expect_error(plot(scored, top_n = 0, draw = FALSE), "top_n")
  expect_error(plot(scored, main = character(), draw = FALSE), "main")
  expect_error(plot(scored, main = NA_character_, draw = FALSE), "main")
})

test_that("portable score print handles a batch with no score coordinates", {
  fixture <- calibration_public_fixture()
  rows <- fixture$rows
  rows$Score <- NA_real_
  scored <- score_mfrm_calibration(
    fixture$frozen, rows, missing_response = "omit"
  )

  printed <- capture.output(print(scored))
  expect_true(any(grepl(
    "No score coordinate is available", printed, fixed = TRUE
  )))
  expect_false(any(grepl("plot(x)", printed, fixed = TRUE)))
  expect_error(
    print(structure(list(), class = "summary.mfrm_calibration_score")),
    "not a complete portable calibration score summary",
    fixed = TRUE
  )
  malformed_summary <- summary(scored)
  expect_identical(malformed_summary$estimates, scored$estimates)
  malformed_summary$estimates$Person <- NULL
  expect_error(
    print(malformed_summary),
    "not a complete portable calibration score summary",
    fixed = TRUE
  )
})

test_that("public loader explicitly refuses incompatible schema fixtures", {
  fixture <- calibration_public_fixture()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  newer <- fixture$frozen
  newer$header$schema_version <- 7L
  saveRDS(newer, path, version = 3)
  newer_error <- tryCatch(
    load_mfrm_calibration(path),
    mfrm_calibration_error = identity
  )
  expect_s3_class(newer_error, "mfrm_calibration_error")
  expect_identical(newer_error$code, "SCHEMA_VERSION_UNSUPPORTED")

  unknown <- fixture$frozen
  unknown$future_section <- list(value = 1L)
  saveRDS(unknown, path, version = 3)
  unknown_error <- tryCatch(
    load_mfrm_calibration(path),
    mfrm_calibration_error = identity
  )
  expect_s3_class(unknown_error, "mfrm_calibration_error")
  expect_identical(unknown_error$code, "SCHEMA_FIELD_UNEXPECTED")

  partial <- fixture$frozen
  partial$response$score_map <- NULL
  saveRDS(partial, path, version = 3)
  partial_error <- tryCatch(
    load_mfrm_calibration(path),
    mfrm_calibration_error = identity
  )
  expect_s3_class(partial_error, "mfrm_calibration_error")
  expect_true(partial_error$code %in% c(
    "SCHEMA_FIELD_MISSING", "SCORE_MAP_INVALID"
  ))
})

test_that("public extraction refusals give actionable non-internal alternatives", {
  fixture <- calibration_public_fixture()
  cases <- list(
    family = list(
      code = "MODEL_FAMILY_UNSUPPORTED",
      mutate = function(fit) {
        fit$config$model <- "unsupported"
        fit
      }
    ),
    estimator = list(
      code = "MODEL_ESTIMATOR_UNSUPPORTED",
      mutate = function(fit) {
        fit$config$method <- "unsupported"
        fit
      }
    ),
    population = list(
      code = "SCORING_BASIS_UNSUPPORTED",
      mutate = function(fit) {
        fit$config$population_spec <- list(active = TRUE)
        fit
      }
    )
  )
  for (case in cases) {
    error <- tryCatch(
      extract_mfrm_calibration(case$mutate(fixture$fit)),
      mfrm_calibration_error = identity
    )
    expect_s3_class(error, "mfrm_calibration_error")
    expect_identical(error$code, case$code)
    expect_match(conditionMessage(error), "fitted-object scoring", fixed = TRUE)
    expect_false(grepl("OPT-[0-9]|CORE-[0-9]|G[0-6] exit|internal",
                       conditionMessage(error), perl = TRUE))
  }
})

test_that("public extraction requires user-reviewed highest-grid fit", {
  fixture <- calibration_public_fixture()
  review <- fixture$quadrature_review

  missing <- tryCatch(
    extract_mfrm_calibration(fixture$fit),
    mfrm_calibration_error = identity
  )
  expect_s3_class(missing, "mfrm_calibration_error")
  expect_identical(missing$code, "QUADRATURE_REVIEW_REQUIRED")

  wrong_type <- tryCatch(
    extract_mfrm_calibration(fixture$fit, quadrature_review = list()),
    mfrm_calibration_error = identity
  )
  expect_identical(wrong_type$code, "QUADRATURE_REVIEW_INVALID")

  lower <- tryCatch(
    extract_mfrm_calibration(
      review$fits$q5, quadrature_review = review
    ),
    mfrm_calibration_error = identity
  )
  expect_identical(lower$code, "QUADRATURE_SOURCE_NOT_HIGHEST")

  incomplete <- review
  incomplete$runs$EstimationConverged[1L] <- FALSE
  unfinished <- tryCatch(
    extract_mfrm_calibration(
      fixture$fit, quadrature_review = incomplete
    ),
    mfrm_calibration_error = identity
  )
  expect_identical(unfinished$code, "QUADRATURE_REVIEW_INCOMPLETE")

  expect_s3_class(
    extract_mfrm_calibration(
      fixture$fit, quadrature_review = review,
      calibration_id = "reviewed-highest-grid"
    ),
    "mfrm_calibration"
  )
})

test_that("portable extraction checks the retained fits behind a grid review", {
  fixture <- calibration_public_fixture()
  review <- fixture$quadrature_review
  changes <- list(
    duplicate_grid = function(x) { x$fits$q5 <- x$fits$q7; x },
    anchors = function(x) {
      x$fits$q5$config$replay_inputs$anchors <- data.frame(
        Facet = "Rater", Level = "R01", Value = 0.1
      )
      x
    },
    resolved_sign = function(x) {
      x$fits$q5$config$facet_signs[1L] <- 1
      x
    },
    score_map = function(x) {
      x$fits$q5$prep$score_map$OriginalScore <-
        x$fits$q5$prep$score_map$OriginalScore + 10
      x
    },
    integration = function(x) {
      x$fits$q5$config$estimation_control$mml_integration <- "adaptive"
      x
    },
    optimizer = function(x) {
      x$fits$q5$config$replay_inputs$reltol <- 0.1
      x
    },
    convergence = function(x) { x$fits$q5$opt$convergence <- 1L; x }
  )
  for (name in names(changes)) {
    error <- tryCatch(
      extract_mfrm_calibration(fixture$fit,
                               quadrature_review = changes[[name]](review)),
      mfrm_calibration_error = identity
    )
    expect_s3_class(error, "mfrm_calibration_error")
    expect_identical(error$code, if (name == "convergence") {
      "QUADRATURE_REVIEW_INCOMPLETE"
    } else "QUADRATURE_REVIEW_INVALID", info = name)
  }
  reordered <- review
  reordered$fits$q5$prep$data <-
    reordered$fits$q5$prep$data[nrow(reordered$fits$q5$prep$data):1L, ]
  expect_s3_class(
    extract_mfrm_calibration(fixture$fit, quadrature_review = reordered),
    "mfrm_calibration"
  )
})

test_that("portable interval identity survives saved results, summaries and CSV", {
  fixture <- calibration_public_fixture()
  for (algorithm in c("quadrature_eap_v1", "quadrature_eap_v2",
                      "adaptive_quadrature_eap_v1", "adaptive_quadrature_eap_v2")) {
    draft <- fixture$draft
    draft$scoring_basis$scoring_algorithm <- algorithm
    draft$integrity$semantic_components <-
      mfrmr:::mfrmr_calibration_semantic_components(draft)
    artifact <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    scored <- score_mfrm_calibration(artifact, fixture$rows, interval_level = 0.955)
    summarized <- summary(scored, digits = 1)
    fields <- c("ScoringAlgorithm", "IntervalLevel")
    expect_identical(unique(scored$estimates$ScoringAlgorithm), algorithm)
    expect_identical(unique(scored$estimates$IntervalLevel), 0.955)
    expect_identical(summarized$estimates[fields], scored$estimates[fields])

    path <- tempfile(fileext = ".csv")
    utils::write.csv(summarized$estimates, path, row.names = FALSE)
    exported <- utils::read.csv(path, stringsAsFactors = FALSE)
    unlink(path)
    expect_identical(exported$Person, scored$estimates$Person)
    expect_identical(as.list(exported[fields]), as.list(scored$estimates[fields]))

    older <- scored
    older$estimates[fields] <- NULL
    expect_identical(summary(older, digits = 1)$estimates, summarized$estimates)
    expect_identical(scored$estimates[c("Estimate", "SD", "Lower", "Upper")],
                     older$estimates[c("Estimate", "SD", "Lower", "Upper")])
    interval <- plot(older, draw = FALSE)
    expected <- if (endsWith(algorithm, "_v1")) {
      "posterior mass may differ from the requested level"
    } else "continuous posterior quantiles"
    displayed <- gsub("[[:space:]]+", " ", paste(capture.output({
      print(older)
      print(summary(older))
    }), collapse = " "))
    expect_match(displayed, expected, fixed = TRUE)
    expect_match(displayed, "95.5% intervals", fixed = TRUE)
    expect_match(interval$data$caption, expected, fixed = TRUE)
    expect_match(interval$data$caption, "95.5% intervals", fixed = TRUE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      gg <- as_ggplot(interval)
      expect_match(gsub("[[:space:]]+", " ", paste(plot_data(gg, "notes")$Text, collapse = "\n")), expected, fixed = TRUE)
    }
  }
  rows <- fixture$rows
  rows$Score <- NA_real_
  empty <- score_mfrm_calibration(fixture$frozen, rows, missing_response = "omit")
  expect_identical(empty$estimates$IntervalLevel, numeric(0))
  expect_identical(summary(empty)$estimates$ScoringAlgorithm, character(0))
  expect_match(paste(capture.output(print(summary(empty))), collapse = " "),
               "no valid responses", fixed = TRUE)
})

test_that("public documentation surfaces share the bounded calibration wording", {
  root <- normalizePath(find.package("mfrmr"), winslash = "/")
  readme <- file.path(root, "README.md")
  if (!file.exists(readme)) readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README requires the source tree on this R installation.")
  article_candidates <- c(file.path(root, c(
      "doc/mfrmr-portable-calibration.Rmd",
      "vignettes/mfrmr-portable-calibration.Rmd"
    )), test_path("..", "..", "vignettes", "mfrmr-portable-calibration.Rmd"))
  article <- article_candidates[file.exists(article_candidates)]
  expect_true(file.exists(readme))
  expect_true(file.exists(file.path(root, "NEWS.md")))
  skip_if(!length(article), "Article requires a source tree or an installation with built vignettes.")
  article <- article[1L]

  text <- list(
    README = paste(
      readLines(readme, warn = FALSE), collapse = "\n"
    ),
    NEWS = paste(
      readLines(file.path(root, "NEWS.md"), warn = FALSE), collapse = "\n"
    ),
    article = paste(readLines(article, warn = FALSE), collapse = "\n")
  )

  # The README routes readers to the specialist guide instead of duplicating it.
  expect_match(text$README, "mfrmr-portable-calibration.html", fixed = TRUE)
  core_surfaces <- text["article"]
  for (surface in core_surfaces) {
    expect_match(surface, "RSM", fixed = TRUE)
    expect_match(surface, "PCM", fixed = TRUE)
    expect_match(surface, "MML", fixed = TRUE)
    expect_match(surface, "fixed\\s+standard-normal", perl = TRUE)
    expect_match(surface, "GPCM", fixed = TRUE)
    expect_match(surface, "JML", fixed = TRUE)
  }
  expect_match(text$article, "mfrm_calibration_capabilities()", fixed = TRUE)
  expect_match(
    text$article,
    "score_mfrm_calibration(calibration, new_rows)",
    fixed = TRUE
  )

  public_text <- paste(unlist(text, use.names = FALSE), collapse = "\n")
  expect_false(grepl("mfrmr:::|CORE-[0-9]|G[0-6] exit|PublicAPIAuthorized",
                     public_text, perl = TRUE))
})

test_that("installed public API preserves CSV identities and refusals in a fresh process", {
  package_path <- normalizePath(find.package("mfrmr"), winslash = "/")
  # Check libraries can contain package symlinks into another installation.
  library_packages <- normalizePath(file.path(.libPaths(), "mfrmr"),
                                   winslash = "/", mustWork = FALSE)
  package_is_installed <- package_path %in% library_packages
  skip_if_not(
    package_is_installed,
    "Fresh-process evidence requires the check-installed package."
  )

  fixture <- calibration_public_fixture()
  artifact_path <- tempfile(fileext = ".rds")
  rows_path <- tempfile(fileext = ".csv")
  result_path <- tempfile(fileext = ".rds")
  script_path <- tempfile(fileext = ".R")
  paths <- c(artifact_path, rows_path, result_path, script_path)
  on.exit(unlink(paths), add = TRUE)
  save_mfrm_calibration(fixture$frozen, artifact_path)
  rows <- fixture$rows
  rows$Person <- ifelse(rows$Person == "NEW_LOW", "001", "1")
  literal_na <- rows[rows$Person == "001", , drop = FALSE]
  literal_na$Person <- "NA"
  rows <- rbind(rows, literal_na)
  utils::write.csv(rows, rows_path, row.names = FALSE, na = "")

  script <- c(
    "args <- commandArgs(trailingOnly = TRUE)",
    ".libPaths(c(args[1L], .libPaths()))",
    "library(mfrmr)",
    "calibration <- load_mfrm_calibration(args[2L])",
    'rows <- read.csv(args[3L], colClasses = "character", na.strings = "",',
    '                 check.names = FALSE, fileEncoding = "UTF-8-BOM")',
    'stopifnot(identical(unique(rows$Person), c("001", "1", "NA")))',
    "result <- score_mfrm_calibration(calibration, rows, adaptive_quad_points = c(15L, 31L), cut_scores = c(Low = -0.513456789, High = 0.823456789))",
    "saveRDS(result, args[4L], version = 3)",
    'unknown <- rows; unknown$Rater[1L] <- "unseen-rater"',
    'err <- tryCatch(score_mfrm_calibration(calibration, unknown), mfrm_calibration_error = identity)',
    'stopifnot(inherits(err, "mfrm_calibration_error"), identical(err$code, "SCORING_FACET_LEVEL_UNKNOWN"))',
    'unknown <- rows; unknown$Score[1L] <- "999"',
    'err <- tryCatch(score_mfrm_calibration(calibration, unknown), mfrm_calibration_error = identity)',
    'stopifnot(inherits(err, "mfrm_calibration_error"), identical(err$code, "SCORING_SCORE_UNKNOWN"))'
  )
  expect_false(any(grepl(":::|mfrmr_", script)))
  writeLines(script, script_path, useBytes = TRUE)
  rscript <- file.path(
    R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  output <- suppressWarnings(system2(
    rscript,
    c(
      "--vanilla", shQuote(script_path), shQuote(dirname(package_path)),
      shQuote(artifact_path), shQuote(rows_path), shQuote(result_path)
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L, info = paste(output, collapse = "\n"))
  expect_true(file.exists(result_path))
  result <- readRDS(result_path)
  expect_s3_class(result, "mfrm_calibration_score")
  expect_setequal(result$estimates$Person, c("001", "1", "NA"))
  expect_identical(result$settings$calibration_id, "public-api-rsm")
  expect_identical(result$settings$engine_identity, "artifact_coordinates_v1")
  expect_identical(unique(result$estimates$ScoringAlgorithm),
                   fixture$frozen$scoring_basis$scoring_algorithm)
  expect_identical(unique(result$estimates$IntervalLevel), 0.95)
  expect_equal(result$estimates,
               score_mfrm_calibration(fixture$frozen, rows)$estimates,
               tolerance = 0)
  expect_equal(nrow(result$quadrature_review), nrow(result$estimates) * 2L)
  expect_true(all(result$quadrature_review$Status == "computed"))
  expect_setequal(result$decisions$Person, c("001", "1", "NA"))
  expect_identical(result$settings$cut_scores, c(Low=-0.513456789,High=0.823456789))
  expect_identical(apa_table(result,which="decisions")$table,result$decisions)
  expect_identical(build_summary_table_bundle(summary(result))$tables$decisions,
    result$decisions)
})


test_that("current RSM portable artifacts accept explicit scoring priors without mutation", {
  f <- calibration_public_fixture()
  expect_error(score_mfrm_calibration(f$frozen, f$rows,
    scoring_prior = list(mean = .4, sd = 1.2)), "integration did not pass")
  f$frozen <- freeze_mfrm_calibration(validate_mfrm_calibration(
    extract_mfrm_calibration(f$fit, quadrature_review = f$quadrature_review,
      scoring_quad_points = 121)))
  before <- serialize(f$frozen, NULL)
  prior <- list(mean = .4, sd = 1.2)
  scored <- score_mfrm_calibration(f$frozen, f$rows, scoring_prior = prior)
  native <- predict_mfrm_units(f$fit, f$rows, scoring_prior = prior,
    scoring_quad_points = 121, readiness_policy = "review")
  expect_equal(scored$estimates$Estimate, unname(native$estimates$Estimate[match(scored$estimates$Person, native$estimates$Person)]), tolerance = 1e-5)
  expect_equal(scored$estimates$SD, unname(native$estimates$SD[match(scored$estimates$Person, native$estimates$Person)]), tolerance = 1e-5)
  expect_true(all(scored$estimates$ScoreIntegrationReady))
  expect_identical(serialize(f$frozen, NULL), before)
  expect_true(all(summary(scored)$estimates$PriorSD == prior$sd))
  legacy <- f$frozen
  legacy$scoring_basis$scoring_algorithm <- sub("_v2$", "_v1", legacy$scoring_basis$scoring_algorithm)
  legacy$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(legacy)
  expect_error(score_mfrm_calibration(legacy, f$rows, scoring_prior = prior), "re-extract")
})


test_that("the output guide agrees with conditional portable GPCM capabilities", {
  guide <- mfrmr_output_guide("calibration")
  cuts <- guide$ObjectRole == "conditional ability decisions"
  decisions <- cuts | guide$MainFunction == "classification_table()"
  expect_true(all(guide$GPCMStatus[!decisions] == "supported_with_caveat"))
  expect_true(all(guide$GPCMStatus[decisions] == "unavailable; ordinary RSM/PCM only"))
  expect_identical(guide$GPCMStatus[cuts], "unavailable; ordinary RSM/PCM only")
  expect_true(any(grepl("source checks", guide$UseWhen, fixed = TRUE)))
  expect_true(any(grepl("GPCM", guide$Notes, fixed = TRUE)))
  boundary <- guide$DecisionBoundary[guide$ObjectRole == "portable calibration lifecycle"]
  expect_length(boundary, 1L)
  expect_match(boundary, "GPCM MML", fixed = TRUE)
  expect_match(boundary, "shared-owner GPCM JML", fixed = TRUE)
})

test_that('quadrature review matches retained fits by grid name rather than list position', {
  fixture <- calibration_public_fixture()
  review <- fixture$quadrature_review
  review$fits <- rev(review$fits)
  expect_no_error(extract_mfrm_calibration(fixture$fit,quadrature_review=review))
  duplicate <- review
  names(duplicate$fits) <- rep(names(review$fits)[1],length(review$fits))
  expect_error(extract_mfrm_calibration(fixture$fit,quadrature_review=duplicate),
    'QUADRATURE_REVIEW_INVALID')
  mislabeled <- review
  names(mislabeled$fits) <- rev(names(review$fits))
  expect_error(extract_mfrm_calibration(fixture$fit,quadrature_review=mislabeled),
    'QUADRATURE_REVIEW_INVALID')
})
