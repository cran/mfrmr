make_contrast_design <- function(n_obs, n_steps, n_parameters,
                                 density, seed) {
  set.seed(seed)
  if (identical(density, 0)) {
    return(Matrix::Matrix(
      0, nrow = n_obs * n_steps, ncol = n_parameters,
      sparse = TRUE
    ))
  }
  Matrix::rsparsematrix(
    nrow = n_obs * n_steps,
    ncol = n_parameters,
    density = density
  )
}

expect_preallocated_contrast_identity <- function(design, score, n_steps) {
  n_obs <- length(score)
  preallocated <- mfrmr:::mfrmr_jml_observed_contrast_design(
    design, score, n_obs, n_steps
  )
  reference <- mfrmr:::mfrmr_jml_observed_contrast_design(
    design, score, n_obs, n_steps,
    implementation = "reference"
  )
  expect_s4_class(preallocated, "dgCMatrix")
  expect_identical(preallocated, reference)
  expect_identical(dim(preallocated), c(n_obs * n_steps, ncol(design)))
}

contrast_problem_from_fit <- function(fit) {
  config <- fit$config
  idx <- mfrmr:::build_indices(
    fit$prep,
    step_facet = config$step_facet,
    slope_facet = config$slope_facet,
    interaction_specs = config$interaction_specs
  )
  sizes <- mfrmr:::build_param_sizes(config)
  params <- mfrmr:::expand_params(fit$opt$par, sizes, config)
  adjacent <- mfrmr:::mfrmr_estimability_adjacent_design(
    prep = fit$prep,
    idx = idx,
    config = config,
    sizes = sizes,
    include_person = TRUE,
    include_population_beta = FALSE
  )
  list(
    prep = fit$prep,
    idx = idx,
    config = config,
    sizes = sizes,
    params = params,
    adjacent = adjacent
  )
}

expect_fit_contrast_identity <- function(fit) {
  problem <- contrast_problem_from_fit(fit)
  preallocated <- mfrmr:::mfrmr_jml_observed_contrast_design(
    problem$adjacent$design,
    problem$idx$score_k,
    problem$adjacent$observation_rows,
    problem$adjacent$transitions
  )
  reference <- mfrmr:::mfrmr_jml_observed_contrast_design(
    problem$adjacent$design,
    problem$idx$score_k,
    problem$adjacent$observation_rows,
    problem$adjacent$transitions,
    implementation = "reference"
  )
  expect_identical(preallocated, reference)

  default_shared <- mfrmr:::mfrmr_jml_recession_shared_geometry(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params
  )
  original <- mfrmr:::mfrmr_jml_observed_contrast_design
  testthat::local_mocked_bindings(
    mfrmr_jml_observed_contrast_design = function(...) {
      original(..., implementation = "reference")
    },
    .package = "mfrmr"
  )
  reference_shared <- mfrmr:::mfrmr_jml_recession_shared_geometry(
    prep = problem$prep,
    idx = problem$idx,
    config = problem$config,
    sizes = problem$sizes,
    params = problem$params
  )
  expect_identical(default_shared, reference_shared)
}

test_that("preallocated contrasts equal the reference across sparse grids", {
  cases <- expand.grid(
    n_obs = c(1L, 2L, 17L, 101L),
    n_steps = c(1L, 2L, 4L, 7L, 10L),
    density = c(0, 0.03, 0.35, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  for (case in seq_len(nrow(cases))) {
    n_obs <- cases$n_obs[case]
    n_steps <- cases$n_steps[case]
    design <- make_contrast_design(
      n_obs = n_obs,
      n_steps = n_steps,
      n_parameters = 13L,
      density = cases$density[case],
      seed = 255000L + case
    )
    patterns <- list(
      rep.int(0L, n_obs),
      rep.int(n_steps, n_obs),
      rep(c(0L, n_steps), length.out = n_obs),
      sample.int(n_steps + 1L, n_obs, replace = TRUE) - 1L
    )
    for (score in patterns) {
      expect_preallocated_contrast_identity(design, score, n_steps)
    }
  }

  empty_columns <- Matrix::Matrix(
    0, nrow = 12L, ncol = 0L, sparse = TRUE
  )
  expect_preallocated_contrast_identity(
    empty_columns, c(0L, 2L, 1L, 2L), n_steps = 3L
  )
})

test_that("preallocated contrasts preserve malformed-input failures", {
  valid <- make_contrast_design(4L, 3L, 5L, 0.2, 255100L)
  invalid <- list(
    list(design = valid, score = c(0L, 1L, NA, 2L), n_obs = 4L,
         n_steps = 3L),
    list(design = valid, score = c(0L, 1L, 2L), n_obs = 4L,
         n_steps = 3L),
    list(design = valid, score = c(0L, 1L, 4L, 2L), n_obs = 4L,
         n_steps = 3L),
    list(design = valid[-1, , drop = FALSE], score = c(0L, 1L, 2L, 3L),
         n_obs = 4L, n_steps = 3L)
  )
  for (case in invalid) {
    for (implementation in c("preallocated", "reference")) {
      expect_error(
        mfrmr:::mfrmr_jml_observed_contrast_design(
          case$design, case$score, case$n_obs, case$n_steps,
          implementation = implementation
        ),
        "malformed adjacent design",
        fixed = TRUE
      )
    }
  }
})

test_that("two-rater sparse anchored PCM preserves full shared geometry", {
  skip_if_not_installed("lpSolve")
  data <- simulate_mfrm_data(
    n_person = 24,
    n_rater = 2,
    n_criterion = 4,
    raters_per_person = 2,
    score_levels = 6,
    seed = 255201L,
    model = "PCM",
    step_facet = "Criterion"
  )
  data$Weight <- 1
  data$Weight[seq(1L, nrow(data), by = 19L)] <- 0
  data$Weight[seq(2L, nrow(data), by = 23L)] <- 2
  data$Score[seq(3L, nrow(data), by = 5L)] <- 0
  data$Score[seq(7L, nrow(data), by = 31L)] <- NA
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
      Facet = "Rater", Level = "R1", Anchor = 0,
      stringsAsFactors = FALSE
    ),
    group_anchors = data.frame(
      Facet = "Criterion", Level = c("C1", "C2"),
      Group = "CriterionGroup", GroupValue = 0,
      stringsAsFactors = FALSE
    ),
    maxit = 30
  ))
  expect_fit_contrast_identity(fit)
})

test_that("interaction RSM and bounded GPCM preserve full shared geometry", {
  skip_if_not_installed("lpSolve")
  rsm_data <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 3,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 5,
    seed = 255202L,
    model = "RSM"
  )
  rsm_data$Score[seq(2L, nrow(rsm_data), by = 4L)] <- 0
  rsm <- suppressWarnings(fit_mfrm(
    rsm_data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "JML",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 1,
    maxit = 30
  ))
  expect_fit_contrast_identity(rsm)

  gpcm_data <- simulate_mfrm_data(
    n_person = 20,
    n_rater = 2,
    n_criterion = 3,
    raters_per_person = 2,
    score_levels = 5,
    seed = 255203L,
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = "Criterion"
  )
  gpcm_data$Score[seq(1L, nrow(gpcm_data), by = 3L)] <- 0
  gpcm <- suppressWarnings(fit_mfrm(
    gpcm_data,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    model = "GPCM",
    method = "JML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    maxit = 30
  ))
  expect_fit_contrast_identity(gpcm)
})
