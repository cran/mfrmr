ic_test_ready <- function(fit) {
  fit$summary$Converged[1] <- TRUE
  fit$summary$InferenceReady[1] <- TRUE
  fit$summary$ConvergenceCode[1] <- 0L
  fit$summary$ConvergenceStatus[1] <- "converged"
  fit$summary$ConvergenceSeverity[1] <- "pass"
  fit
}

ic_test_refresh_contract <- function(fit) {
  contract <- mfrmr:::build_mfrm_ic_contract(
    loglik = -as.numeric(fit$opt$value[1]),
    npar = length(fit$opt$par),
    prep = fit$prep,
    config = fit$config,
    method = fit$config$method
  )
  for (field in intersect(names(contract), names(fit$summary))) {
    fit$summary[[field]][1] <- contract[[field]]
  }
  fit
}

ic_test_design <- function() {
  data <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = paste0("C", 1:3),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data$Score <- (seq_len(nrow(data)) - 1L) %% 4L
  mfrmr:::prepare_mfrm_data(
    data,
    person_col = "Person",
    facet_cols = c("Rater", "Criterion"),
    score_col = "Score",
    rating_min = 0,
    rating_max = 3
  )
}

ic_test_config <- function(prep,
                           model = "RSM",
                           method = "MML",
                           interaction_specs = list(),
                           anchors = NULL,
                           dummy_facets = character(0),
                           noncenter_facet = "Person") {
  signs <- mfrmr:::build_facet_signs(prep$facet_names)
  mfrmr:::build_estimation_config(
    prep = prep,
    model = model,
    method = method,
    step_facet = if (identical(model, "RSM")) NULL else "Criterion",
    slope_facet = if (identical(model, "GPCM")) "Criterion" else NULL,
    interaction_specs = interaction_specs,
    weight_col = NULL,
    facet_signs = signs$signs,
    positive_facets = signs$positive_facets,
    noncenter_facet = noncenter_facet,
    dummy_facets = dummy_facets,
    anchor_df = anchors,
    group_anchor_df = NULL,
    population = NULL
  )
}

test_that("the common MML information-criterion panel uses Persons exactly", {
  persons <- rep(paste0("P", seq_len(100)), each = 12L)
  prep <- list(data = data.frame(
    Person = persons,
    Weight = rep(1, length(persons)),
    stringsAsFactors = FALSE
  ))
  config <- list(
    weight_col = NULL,
    estimation_control = list(quad_points = 31L)
  )
  contract <- mfrmr:::build_mfrm_ic_contract(
    loglik = -400,
    npar = 12L,
    prep = prep,
    config = config,
    method = "MML"
  )

  expect_identical(contract$ICContractVersion, "mfrmr_ic_person_v2")
  expect_equal(contract$ResponseRows, 1200L)
  expect_equal(contract$WeightedResponseTotal, 1200)
  expect_equal(contract$Persons, 100L)
  expect_equal(contract$ICSampleSize, 100)
  expect_identical(contract$ICSampleSizeBasis, "person_count")
  expect_equal(contract$Deviance, 800)
  expect_equal(contract$AIC, 824)
  expect_equal(contract$BIC, 800 + log(100) * 12, tolerance = 1e-12)
  expect_equal(
    contract$SABIC,
    800 + log((100 + 2) / 24) * 12,
    tolerance = 1e-12
  )
  expect_true(contract$ICEligible)
  expect_true(contract$ICSelectable)
  expect_true(contract$ICIntegrationSelectable)
  expect_identical(contract$ICIntegrationTier, "standard_start")
  expect_identical(contract$ICIntegrationStatus, "selection_start")
  expect_equal(contract$ICQuadraturePoints, 31L)
  expect_true(contract$SABICSelectable)
  expect_identical(contract$AICFormula, "aic_deviance_plus_2k")
  expect_identical(contract$BICFormula, "bic_person_count")
  expect_identical(contract$SABICFormula, "sclove_n_plus_2_over_24")
  expect_identical(contract$IntegrationEvaluationId, "mfrmr_ghq_normal_v1:q=31")

  doubled_rows <- prep
  doubled_rows$data <- rbind(prep$data, prep$data)
  doubled <- mfrmr:::build_mfrm_ic_contract(
    -400, 12L, doubled_rows, config, "MML"
  )
  expect_equal(doubled$ResponseRows, 2400L)
  expect_equal(doubled$Persons, 100L)
  expect_equal(doubled[c("AIC", "BIC", "SABIC")],
               contract[c("AIC", "BIC", "SABIC")])

  missing_layout <- prep
  missing_layout$data <- prep$data[rep(c(TRUE, TRUE, FALSE), 400L), , drop = FALSE]
  missing_contract <- mfrmr:::build_mfrm_ic_contract(
    -400, 12L, missing_layout, config, "MML"
  )
  expect_equal(missing_contract$ResponseRows, 800L)
  expect_equal(missing_contract$Persons, 100L)
  expect_equal(missing_contract[c("AIC", "BIC", "SABIC")],
               contract[c("AIC", "BIC", "SABIC")])
})

test_that("IC eligibility fails closed for JML, non-unit weights, and small-N SABIC", {
  base <- data.frame(
    Person = rep(paste0("P", 1:30), each = 4L),
    Weight = 1,
    stringsAsFactors = FALSE
  )
  config <- list(weight_col = NULL, estimation_control = list(quad_points = 31L))

  explicit <- list(data = base)
  explicit_config <- config
  explicit_config$weight_col <- "W"
  explicit_contract <- mfrmr:::build_mfrm_ic_contract(
    -120, 7L, explicit, explicit_config, "MML"
  )
  expect_identical(explicit_contract$WeightPolicy, "explicit_unit")
  expect_true(explicit_contract$ICEligible)
  expect_true(explicit_contract$ICSelectable)

  constant <- list(data = base)
  constant$data$Weight <- rep(rep(c(1, 2), each = 4L), length.out = nrow(base))
  constant_contract <- mfrmr:::build_mfrm_ic_contract(
    -120, 7L, constant, explicit_config, "MML"
  )
  expect_identical(
    constant_contract$WeightPolicy,
    "nonunit_constant_within_person"
  )
  expect_false(constant_contract$ICEligible)
  expect_true(all(is.na(unlist(
    constant_contract[c("AIC", "BIC", "SABIC")],
    use.names = FALSE
  ))))

  varying <- list(data = base)
  varying$data$Weight <- rep(c(1, 2), length.out = nrow(base))
  varying_contract <- mfrmr:::build_mfrm_ic_contract(
    -120, 7L, varying, explicit_config, "MML"
  )
  expect_identical(varying_contract$WeightPolicy, "nonunit_row_varying")
  expect_false(varying_contract$ICEligible)

  jml_contract <- mfrmr:::build_mfrm_ic_contract(
    -120, 7L, explicit, config, "JML"
  )
  expect_identical(jml_contract$ICStatus, "descriptive_jml")
  expect_false(jml_contract$ICEligible)
  expect_true(all(is.na(unlist(
    jml_contract[c("AIC", "BIC", "SABIC")],
    use.names = FALSE
  ))))
  expect_true(is.finite(jml_contract$LegacyAIC))
  expect_true(is.finite(jml_contract$LegacyBIC))
  jml_text <- paste(
    mfrmr:::mfrm_ic_console_lines(as.data.frame(jml_contract)),
    collapse = " "
  )
  expect_match(jml_text, "not eligible (descriptive_jml)", fixed = TRUE)
  expect_match(jml_text, "Legacy descriptive AIC", fixed = TRUE)
  expect_match(jml_text, "not part of the common MML ranking panel", fixed = TRUE)

  small <- explicit
  small$data <- small$data[small$data$Person %in% paste0("P", 1:22), ]
  small_contract <- mfrmr:::build_mfrm_ic_contract(
    -40, 4L, small, config, "MML"
  )
  expect_true(small_contract$ICEligible)
  expect_true(small_contract$ICSelectable)
  expect_false(small_contract$SABICSelectable)
  expect_identical(small_contract$ICStatus, "ok_sabic_nonpositive")
  expect_equal(small_contract$SABIC, 80, tolerance = 1e-12)
  expect_match(
    paste(mfrmr:::mfrm_ic_console_lines(as.data.frame(small_contract)),
          collapse = " "),
    "automatic SABIC selection is disabled"
  )
})

test_that("IC integration tiers retain raw values but fail closed below q=31", {
  prep <- list(data = data.frame(
    Person = rep(paste0("P", 1:30), each = 4L),
    Weight = 1,
    stringsAsFactors = FALSE
  ))
  expected <- data.frame(
    q = c(7L, 15L, 31L, 61L),
    tier = c(
      "coarse_screening", "intermediate_review", "standard_start",
      "dense_sensitivity"
    ),
    status = c(
      "screening_only", "review_only", "selection_start",
      "dense_sensitivity"
    ),
    selectable = c(FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  contracts <- lapply(expected$q, function(q) {
    mfrmr:::build_mfrm_ic_contract(
      loglik = -120,
      npar = 7L,
      prep = prep,
      config = list(
        weight_col = NULL,
        estimation_control = list(quad_points = q)
      ),
      method = "MML"
    )
  })
  expect_equal(vapply(contracts, `[[`, integer(1), "ICQuadraturePoints"),
               expected$q)
  expect_identical(vapply(contracts, `[[`, character(1), "ICIntegrationTier"),
                   expected$tier)
  expect_identical(vapply(contracts, `[[`, character(1), "ICIntegrationStatus"),
                   expected$status)
  expect_identical(vapply(contracts, `[[`, logical(1), "ICSelectable"),
                   expected$selectable)
  expect_identical(
    vapply(contracts, `[[`, logical(1), "ICIntegrationSelectable"),
    expected$selectable
  )
  expect_true(all(vapply(contracts, `[[`, logical(1), "ICEligible")))
  expect_true(all(vapply(contracts, function(contract) {
    all(is.finite(unlist(contract[c("AIC", "BIC", "SABIC")])))
  }, logical(1))))
  expect_false(contracts[[1]]$SABICSelectable)
  expect_false(contracts[[2]]$SABICSelectable)
  expect_true(contracts[[3]]$SABICSelectable)
  expect_true(contracts[[4]]$SABICSelectable)
  coarse_text <- paste(
    mfrmr:::mfrm_ic_console_lines(as.data.frame(contracts[[1]])),
    collapse = " "
  )
  expect_match(coarse_text, "screening/review only", fixed = TRUE)
  expect_match(coarse_text, "q=7", fixed = TRUE)
  expect_match(coarse_text, "Use q >= 31", fixed = TRUE)
})

test_that("free-parameter fixtures cover constraints and every 0.2.3 family", {
  prep <- ic_test_design()
  interactions <- mfrmr:::build_facet_interaction_specs(
    prep,
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 0,
    interaction_policy = "silent"
  )
  anchored <- data.frame(
    Facet = "Rater", Level = "R1", Anchor = 0,
    stringsAsFactors = FALSE
  )

  configs <- list(
    ICD_001 = ic_test_config(prep),
    ICD_002 = ic_test_config(prep, model = "PCM"),
    ICD_003 = ic_test_config(prep, model = "GPCM"),
    ICD_004 = ic_test_config(prep, interaction_specs = interactions),
    ICD_005 = ic_test_config(prep, anchors = anchored),
    ICD_006 = ic_test_config(prep, dummy_facets = "Rater"),
    ICD_007 = ic_test_config(prep, method = "JML"),
    ICD_009 = ic_test_config(prep, noncenter_facet = "Rater")
  )
  population <- configs$ICD_001$config
  population$population_spec <- list(
    active = TRUE,
    design_matrix = cbind(Intercept = 1, X = c(-1, 0, 1, 2))
  )
  configs$ICD_008 <- list(
    config = population,
    sizes = mfrmr:::build_param_sizes(population)
  )
  configs <- configs[paste0("ICD_00", 1:9)]

  expected <- data.frame(
    theta = c(0, 0, 0, 0, 0, 0, 4, 0, 0),
    Rater = c(1, 1, 1, 1, 0, 0, 1, 1, 2),
    Criterion = rep(2, 9),
    interactions = c(0, 0, 0, 2, 0, 0, 0, 0, 0),
    steps = c(2, 6, 6, 2, 2, 2, 2, 2, 2),
    log_slopes = c(0, 0, 2, 0, 0, 0, 0, 0, 0),
    beta = c(0, 0, 0, 0, 0, 0, 0, 2, 0),
    log_sigma2 = c(0, 0, 0, 0, 0, 0, 0, 1, 0),
    Npar = c(5, 9, 11, 7, 4, 4, 9, 8, 6)
  )

  observed <- do.call(rbind, lapply(configs, function(x) {
    sizes <- x$sizes
    get_size <- function(name) {
      value <- sizes[[name]]
      if (is.null(value)) 0L else as.integer(value)
    }
    values <- vapply(
      c(
        "theta", "Rater", "Criterion", "interactions", "steps",
        "log_slopes", "beta", "log_sigma2"
      ),
      get_size,
      integer(1)
    )
    c(values, Npar = sum(values))
  }))
  storage.mode(observed) <- "double"

  expect_equal(unname(observed), unname(as.matrix(expected)))
  registry_candidates <- c(
    file.path("inst", "validation", "ic-free-dimension-fixtures-0.2.3.csv"),
    testthat::test_path(
      "..", "..", "inst", "validation",
      "ic-free-dimension-fixtures-0.2.3.csv"
    )
  )
  registry_path <- registry_candidates[file.exists(registry_candidates)][1]
  if (length(registry_path) == 1L && !is.na(registry_path)) {
    registry <- utils::read.csv(registry_path, stringsAsFactors = FALSE)
    expect_identical(registry$FixtureId, sub("_", "-", rownames(observed)))
    registry_counts <- as.matrix(registry[, c(
      "ExpectedTheta", "ExpectedRater", "ExpectedCriterion",
      "ExpectedInteractions", "ExpectedSteps", "ExpectedLogSlopes",
      "ExpectedBeta", "ExpectedLogSigma2", "ExpectedNpar"
    )])
    storage.mode(registry_counts) <- "double"
    expect_equal(unname(registry_counts), unname(observed))
  }
  expect_error(
    mfrmr:::build_estimation_summary(
      model = "RSM",
      method = "MML",
      prep = prep,
      config = configs$ICD_001$config,
      sizes = list(theta = 1L),
      opt = list(par = numeric(0), value = 0)
    ),
    "parameter-dimension mismatch"
  )
})

test_that("fitted-object IC comparison audits weights, legacy state, and identities", {
  data <- load_mfrmr_data("example_core")
  keep_persons <- unique(data$Person)[seq_len(24L)]
  data <- data[data$Person %in% keep_persons, , drop = FALSE]
  common <- list(
    data = data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    quad_points = 31,
    maxit = 10
  )

  fit_unweighted <- suppressWarnings(do.call(fit_mfrm, common))
  unit_data <- data
  unit_data$AnalysisWeight <- 1
  fit_unit <- suppressWarnings(do.call(
    fit_mfrm,
    utils::modifyList(common, list(data = unit_data, weight = "AnalysisWeight"))
  ))
  fit_unweighted <- ic_test_ready(fit_unweighted)
  fit_unit <- ic_test_ready(fit_unit)

  unit_comparison <- compare_mfrm(
    Unweighted = fit_unweighted,
    ExplicitUnit = fit_unit
  )
  expect_true(all(unit_comparison$table$ICComparable))
  expect_true(all(unit_comparison$table$ICSelectable))
  expect_identical(
    as.character(unit_comparison$table$WeightPolicy),
    c("unweighted", "explicit_unit")
  )
  expect_equal(
    unit_comparison$table$LogLik[1],
    unit_comparison$table$LogLik[2],
    tolerance = 1e-8
  )
  expect_equal(
    unit_comparison$table$BIC[1],
    unit_comparison$table$BIC[2],
    tolerance = 1e-8
  )

  constant_data <- data
  person_weight <- setNames(rep(c(1, 2), length.out = 24L), keep_persons)
  constant_data$AnalysisWeight <- person_weight[as.character(constant_data$Person)]
  fit_constant <- suppressWarnings(do.call(
    fit_mfrm,
    utils::modifyList(common, list(
      data = constant_data,
      weight = "AnalysisWeight",
      maxit = 6
    ))
  ))
  fit_constant <- ic_test_ready(fit_constant)
  expect_identical(
    as.character(fit_constant$summary$WeightPolicy[1]),
    "nonunit_constant_within_person"
  )
  expect_false(isTRUE(fit_constant$summary$ICEligible[1]))
  expect_true(all(is.na(unlist(
    fit_constant$summary[1, c("AIC", "BIC", "SABIC")],
    use.names = FALSE
  ))))
  expect_warning(
    constant_comparison <- compare_mfrm(A = fit_constant, B = fit_constant),
    "not eligible for the common information-criterion panel"
  )
  expect_false(any(constant_comparison$table$ICComparable))

  varying_data <- data
  varying_data$AnalysisWeight <- ave(
    seq_len(nrow(varying_data)),
    as.character(varying_data$Person),
    FUN = function(i) rep(c(1, 2), length.out = length(i))
  )
  fit_varying <- suppressWarnings(do.call(
    fit_mfrm,
    utils::modifyList(common, list(
      data = varying_data,
      weight = "AnalysisWeight",
      maxit = 6
    ))
  ))
  expect_identical(
    as.character(fit_varying$summary$WeightPolicy[1]),
    "nonunit_row_varying"
  )
  expect_false(isTRUE(fit_varying$summary$ICEligible[1]))

  coarse_fit <- suppressWarnings(do.call(
    fit_mfrm,
    utils::modifyList(common, list(quad_points = 7, maxit = 10))
  ))
  coarse_fit <- ic_test_ready(coarse_fit)
  expect_warning(
    coarse_comparison <- compare_mfrm(A = coarse_fit, B = coarse_fit),
    "screening/review-only"
  )
  expect_true(all(coarse_comparison$table$ICEligible))
  expect_false(any(coarse_comparison$table$ICSelectable))
  expect_false(any(coarse_comparison$table$ICComparable))
  expect_true(all(is.finite(coarse_comparison$table$AIC)))
  expect_true(all(is.na(coarse_comparison$table$Delta_AIC)))
  expect_length(coarse_comparison$preferred, 0L)
  expect_identical(
    unique(as.character(coarse_comparison$table$ICIntegrationTier)),
    "coarse_screening"
  )

  legacy <- fit_unweighted
  legacy$summary$BIC[1] <- legacy$summary$LegacyBIC[1]
  legacy_fields <- setdiff(
    mfrmr:::mfrm_ic_contract_required_fields(),
    c("N", "Persons", "LogLik", "AIC", "BIC")
  )
  legacy$summary[intersect(
    c(legacy_fields, "LegacyAIC", "LegacyBIC", "LegacyICSampleSize",
      "LegacyICSampleSizeBasis"),
    names(legacy$summary)
  )] <- NULL
  expect_warning(
    legacy_comparison <- compare_mfrm(Current = fit_unweighted, Legacy = legacy),
    "legacy, unknown, or incomplete"
  )
  expect_false(any(legacy_comparison$table$ICComparable))
  expect_identical(
    as.character(legacy_comparison$table$ICContractState[2]),
    "legacy_or_unknown"
  )
  expect_identical(
    as.character(legacy_comparison$table$ICStatus[2]),
    "suppressed_legacy_contract"
  )
  expect_true(is.na(legacy_comparison$table$BIC[2]))
  expect_true(is.finite(legacy_comparison$table$LegacyBIC[2]))

  tampered <- fit_unweighted
  tampered$summary$BIC[1] <- tampered$summary$BIC[1] + 1
  expect_warning(
    tampered_comparison <- compare_mfrm(Current = fit_unweighted, Edited = tampered),
    "stored information-criterion value"
  )
  expect_false(any(tampered_comparison$table$ICComparable))
  expect_false(tampered_comparison$table$StoredICConsistent[2])
  expect_match(tampered_comparison$table$StoredICMismatch[2], "BIC")

  different_integration <- fit_unweighted
  different_integration$config$estimation_control$quad_points <- 61L
  different_integration <- ic_test_refresh_contract(different_integration)
  expect_warning(
    integration_comparison <- compare_mfrm(
      Q31 = fit_unweighted,
      Q61 = different_integration
    ),
    "integration-evaluation identity"
  )
  expect_false(any(integration_comparison$table$ICComparable))
  expect_false(integration_comparison$comparison_basis$same_integration_evaluation)

  different_constraint <- fit_unweighted
  different_constraint$config$noncenter_facet <- "Rater"
  expect_warning(
    constraint_comparison <- compare_mfrm(
      PersonCenteredBasis = fit_unweighted,
      RaterCenteredBasis = different_constraint
    ),
    "different centering constraints"
  )
  expect_false(any(constraint_comparison$table$ICComparable))
  expect_false(constraint_comparison$comparison_basis$same_constraint_basis)
})
