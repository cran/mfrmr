anchor_population_arguments <- function(model = "RSM", sd = 1.4) {
  ratings <- load_mfrmr_data("example_core")
  list(data = ratings, person = "Person", facets = c("Rater", "Criterion"),
    score = "Score", model = model, method = "MML",
    step_facet = if (model == "PCM") "Criterion" else NULL,
    anchors = data.frame(Facet = "Rater", Level = "R01", Anchor = .3),
    noncenter_facet = "Rater", population_formula = ~ 1,
    person_data = unique(ratings["Person"]), population_sd = sd,
    quad_points = 41L, maxit = 400L, reltol = 1e-11, anchor_policy = "silent")
}

test_that("fixed SD is removed from the free map and reduced gradients are correct", {
  for (model in c("RSM", "PCM")) {
    args <- anchor_population_arguments(model)
    args$jml_penalty <- "none" # The ordinary default preserves fixed-SD MML.
    fit <- suppressWarnings(do.call(fit_mfrm, args))
    cfg <- fit$config; sizes <- build_param_sizes(cfg)
    idx <- build_indices(fit$prep, step_facet = cfg$step_facet)
    q <- gauss_hermite_normal(31L)
    p <- fit$opt$par + seq_along(fit$opt$par) / 500
    params <- expand_params(p, sizes, cfg)
    expect_identical(sizes$log_sigma2, 0L)
    expect_identical(fit$population$fixed_sd, 1.4)
    expect_identical(fit$population$sigma2, 1.4^2)
    expect_identical(params$population$sigma2, 1.4^2)
    expect_equal(collapse_expanded_params(params, cfg), p)
    expect_length(build_initial_param_vector(cfg, sizes), length(p))
    objective <- function(x) mfrm_loglik_mml(x, idx, cfg, sizes, q)
    numerical <- vapply(seq_along(p), function(j) {
      delta <- p * 0; delta[j] <- 1e-5
      (objective(p + delta) - objective(p - delta)) / 2e-5
    }, numeric(1))
    expect_equal(as.numeric(mfrm_grad_mml(p, idx, cfg, sizes, q)), numerical,
      tolerance = 2e-7)
    # A common anchor shift changes the origin, not response probabilities.
    shifted_cfg <- cfg
    spec <- shifted_cfg$facet_specs$Rater
    spec$anchors[is.finite(spec$anchors)] <- spec$anchors[is.finite(spec$anchors)] + .6
    shifted_cfg$facet_specs$Rater <- spec
    shifted <- params
    shifted$facets$Rater <- shifted$facets$Rater + .6
    shifted$population$coefficients <- shifted$population$coefficients + .6
    shifted_p <- collapse_expanded_params(shifted, shifted_cfg)
    expect_equal(mfrm_loglik_mml(shifted_p, idx, shifted_cfg, sizes, q),
      objective(p), tolerance = 1e-10)
    expect_equal(mfrm_grad_mml(shifted_p, idx, shifted_cfg, sizes, q),
      mfrm_grad_mml(p, idx, cfg, sizes, q), tolerance = 1e-10)
    expect_identical(summary(fit)$population_overview$ResidualVarianceStatus, "fixed")
    expect_identical(mfrm_fit_scale_contract(fit)$CoordinateBasis,
      "estimated_population_mean_fixed_sd")
  }
})

test_that("anchor-origin identification requires effective anchors and a full rank design", {
  args <- anchor_population_arguments(sd = 1)
  no_anchor <- args; no_anchor$anchors <- NULL
  expect_error(suppressWarnings(do.call(fit_mfrm, no_anchor)), "retained direct anchor")
  wrong <- args; wrong$anchors$Facet <- "Criterion"
  wrong$anchors$Level <- as.character(args$data$Criterion[1])
  expect_error(suppressWarnings(do.call(fit_mfrm, wrong)), "retained direct anchor")
  unknown <- args; unknown$anchors$Level <- "not-observed"
  expect_error(suppressWarnings(do.call(fit_mfrm, unknown)), "retained direct anchor")
  confounded <- args
  d <- args$data
  keep <- match(d$Rater, unique(d$Rater)) == match(d$Criterion, unique(d$Criterion))
  confounded$data <- d[keep, ]
  expect_error(suppressWarnings(do.call(fit_mfrm, confounded)),
    "[Uu]nidentified|[Ee]stimability|[Aa]lias|[Ii]dentifi|rank")
  full <- args
  full$anchors <- data.frame(Facet = "Rater", Level = unique(d$Rater),
    Anchor = c(.3, -.1, .2, .4))
  fit <- suppressWarnings(do.call(fit_mfrm, full))
  expect_equal(build_param_sizes(fit$config)$Rater, 0L)
  expect_true(is.finite(fit$population$coefficients))
})

test_that("unsupported fixed-SD routes fail rather than silently changing likelihood", {
  args <- anchor_population_arguments()
  for (sd in list(0, -1, NA_real_, Inf, 1i, c(1, 2), 1e300, 1e-300)) {
    invalid <- args; invalid$population_sd <- sd
    expect_error(do.call(fit_mfrm, invalid), "population_sd")
  }
  for (changes in list(list(model = "GPCM"), list(method = "JML"),
      list(method = "CML"), list(mml_engine = "em"),
      list(population_formula = NULL), list(facet_interactions = ~ Rater:Criterion))) {
    expect_error(do.call(fit_mfrm, modifyList(args, changes, keep.null = TRUE)), "population_sd")
  }
  penalized <- args; penalized$jml_penalty <- "jeffreys"
  expect_error(do.call(fit_mfrm, penalized), "Jeffreys")
  jml_origin <- args; jml_origin$centering_scope <- "anchor_origin"
  expect_error(do.call(fit_mfrm, jml_origin), "requires ordinary additive RSM/PCM JML")
  no_sd <- args; no_sd$population_sd <- NULL
  expect_error(do.call(fit_mfrm, no_sd), "identification requires")
  covariates <- args; covariates$person_data$X <- seq_len(nrow(covariates$person_data))
  covariates$population_formula <- ~ X
  expect_error(do.call(fit_mfrm, covariates), "intercept-only")
})

test_that("defaults retain free variance and fixed-SD models have distinct identities", {
  args <- anchor_population_arguments()
  args$noncenter_facet <- "Person"; args$anchors <- NULL
  fixed <- suppressWarnings(do.call(fit_mfrm, args))
  args$population_sd <- NULL
  free <- suppressWarnings(do.call(fit_mfrm, args))
  expect_identical(build_param_sizes(free$config)$log_sigma2, 1L)
  expect_null(free$config$population_spec$fixed_sd)
  expect_null(mfrm_population_design(free)$fixed_sd)
  expect_identical(mfrm_population_design(fixed)$fixed_sd, 1.4)
  expect_false(audit_compare_mfrm_nesting(list(fixed, free), c("fixed", "free"))$eligible)
  changed <- fixed; changed$config$population_spec$fixed_sd <- 1
  expect_false(identical(mfrm_population_design(fixed), mfrm_population_design(changed)))
})
