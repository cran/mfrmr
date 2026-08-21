shared_recession_problem <- function(fit) {
  config <- fit$config
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  sizes <- mfrmr:::build_param_sizes(config)
  params <- mfrmr:::expand_params(fit$opt$par, sizes, config)
  list(
    prep = fit$prep,
    idx = idx,
    config = config,
    sizes = sizes,
    params = params,
    person_boundary = fit$config$boundary_audit
  )
}

compare_shared_recession_audits <- function(fit) {
  problem <- shared_recession_problem(fit)
  shared <- mfrmr:::mfrmr_jml_recession_shared_geometry(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params
  )
  expect_identical(
    shared$contract_version,
    "mfrmr-jml-recession-shared-geometry-v1"
  )
  expect_identical(shared$state, "ready")
  expect_true(shared$ready)

  structural_legacy <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    shared_geometry = NULL
  )
  structural_shared <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    shared_geometry = shared
  )
  expect_identical(structural_shared, structural_legacy)

  joint_legacy <- mfrmr:::audit_mfrm_jml_joint_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    person_boundary_audit = problem$person_boundary,
    shared_geometry = NULL
  )
  joint_shared <- mfrmr:::audit_mfrm_jml_joint_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    person_boundary_audit = problem$person_boundary,
    shared_geometry = shared
  )
  expect_identical(joint_shared, joint_legacy)
  invisible(shared)
}

make_shared_geometry_data <- function(seed = 254001L) {
  data <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 3,
    score_levels = 5,
    seed = seed,
    model = "PCM",
    step_facet = "Criterion"
  )
  data$Weight <- 1
  data$Weight[seq(1, nrow(data), by = 29)] <- 0
  data$Weight[seq(2, nrow(data), by = 17)] <- 2
  data$Score[seq(7, nrow(data), by = 31)] <- NA
  data
}

test_that("shared geometry preserves anchored weighted PCM audits exactly", {
  skip_if_not_installed("lpSolve")
  data <- make_shared_geometry_data()
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    weight = "Weight",
    model = "PCM",
    method = "JML",
    step_facet = "Criterion",
    anchors = data.frame(
      Facet = "Rater",
      Level = "R1",
      Anchor = 0,
      stringsAsFactors = FALSE
    ),
    group_anchors = data.frame(
      Facet = "Criterion",
      Level = c("C1", "C2"),
      Group = "CG",
      GroupValue = 0,
      stringsAsFactors = FALSE
    ),
    maxit = 40
  ))
  problem <- shared_recession_problem(fit)
  shared <- compare_shared_recession_audits(fit)

  structural_columns <- which(shared$adjacent$map$Block != "Person")
  direct <- mfrmr:::mfrmr_jml_observed_contrast_design(
    shared$adjacent$design[, structural_columns, drop = FALSE],
    problem$idx$score_k,
    shared$adjacent$observation_rows,
    shared$adjacent$transitions
  )
  expect_identical(
    methods::as(direct, "dgCMatrix"),
    methods::as(
      shared$contrast[, structural_columns, drop = FALSE], "dgCMatrix"
    )
  )
})

test_that("shared geometry preserves interaction RSM and bounded GPCM", {
  skip_if_not_installed("lpSolve")
  data <- make_shared_geometry_data()
  retained <- data[data$Weight > 0 & !is.na(data$Score), , drop = FALSE]

  interaction <- suppressWarnings(fit_mfrm(
    retained,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "JML",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 1,
    maxit = 40
  ))
  compare_shared_recession_audits(interaction)

  gpcm_data <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 3,
    score_levels = 5,
    seed = 254002,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion"
  )
  gpcm <- suppressWarnings(fit_mfrm(
    gpcm_data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 40
  ))
  compare_shared_recession_audits(gpcm)
})

test_that("malformed optional shared geometry falls back to legacy work", {
  skip_if_not_installed("lpSolve")
  data <- make_shared_geometry_data()
  fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    weight = "Weight",
    model = "PCM",
    method = "JML",
    step_facet = "Criterion",
    maxit = 20
  ))
  problem <- shared_recession_problem(fit)
  legacy <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params
  )
  malformed <- list(
    contract_version = "mfrmr-jml-recession-shared-geometry-v1",
    state = "ready",
    ready = TRUE,
    adjacent = list(design = Matrix::Matrix(0), map = data.frame()),
    target_system = list(metadata = data.frame(), expansion = Matrix::Matrix(0)),
    contrast = Matrix::Matrix(0)
  )
  fallback <- mfrmr:::audit_mfrm_jml_structural_recession(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params,
    shared_geometry = malformed
  )
  expect_identical(fallback, legacy)

  testthat::local_mocked_bindings(
    mfrmr_jml_recession_shared_geometry = function(...) {
      list(
        contract_version = "mfrmr-jml-recession-shared-geometry-v1",
        state = "not_evaluated_contrast",
        ready = FALSE,
        detail = "injected shared-construction failure"
      )
    },
    .package = "mfrmr"
  )
  fallback_fit <- suppressWarnings(fit_mfrm(
    data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    weight = "Weight",
    model = "PCM",
    method = "JML",
    step_facet = "Criterion",
    maxit = 20
  ))
  expect_identical(fallback_fit$summary, fit$summary)
  expect_identical(fallback_fit$facets, fit$facets)
  expect_identical(fallback_fit$steps, fit$steps)
  expect_identical(fallback_fit$readiness, fit$readiness)
  expect_identical(
    fallback_fit$config$boundary_audit,
    fit$config$boundary_audit
  )
})
