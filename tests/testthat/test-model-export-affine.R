affine_export_fixture <- function(model = "GPCM", anchored = TRUE, scope = "facet") {
  f <- readRDS(test_path("fixtures", "gmfrm-adaptive-profiles.rds"))$fit
  f$config$n_cat <- 4L
  f$config$rating_max <- f$config$rating_min + 3L
  f$config$score_map <- data.frame(OriginalScore = 0:3, InternalScore = f$config$rating_min + 0:3)
  f$config$model <- model
  f$config$slope_facet <- f$config$facet_names[1L]
  f$config$gpcm_spec <- mfrmr:::build_gpcm_slope_spec(
    f$config$facet_levels[[f$config$slope_facet]], f$config$slope_facet, f$config$step_facet)
  f$config$step_specs <- mfrmr:::mfrmr_default_step_specs(f$config)
  if (anchored) {
    for (facet in f$config$facet_names) {
      levels <- f$config$facet_levels[[facet]]
      groups <- if (length(levels) >= 4L) stats::setNames(rep("g", 2L), levels[2:3]) else NULL
      f$config$facet_specs[[facet]] <- mfrmr:::build_facet_constraint(levels,
        anchors = stats::setNames(.4, levels[1]), groups = groups,
        group_values = if (is.null(groups)) NULL else c(g = -.15),
        centered = TRUE, centering_scope = scope)
    }
    f$config$step_specs <- lapply(f$config$step_specs, function(s)
      mfrmr:::build_step_constraint(s$n_steps, c(`1` = .25), target = .1, scope = s$scope))
  }
  f
}

test_that("affine exports preserve anchored response laws, constraints and prior marginals", {
  for (model in c("RSM", "PCM", "GPCM")) for (scope in c("facet", "remainder")) {
    f <- affine_export_fixture(model, scope = scope)
    f$config$facet_signs[1L] <- 1
    b <- mfrmr:::mfrm_model_export_data(f, 1.3, .7, .4)
    s <- b$stan_data
    sizes <- mfrmr:::build_param_sizes(f$config)
    params <- mfrmr:::expand_params(sin(seq_len(sum(unlist(sizes)))) * .5, sizes, f$config)
    loc <- unlist(params$facets, use.names = FALSE)
    free <- drop(crossprod(s$location_basis, loc - s$location_offset))
    expect_equal(as.vector(s$location_offset) + drop(s$location_basis %*% free), loc, tolerance = 1e-12)
    steps <- if (model == "RSM") matrix(params$steps, 1L) else params$steps_mat
    stacked <- as.vector(t(steps))
    step <- drop(crossprod(s$step_basis, stacked - s$step_offset))
    expect_equal(as.vector(s$step_offset) + drop(s$step_basis %*% step), stacked, tolerance = 1e-12)
    expect_equal(b$levels$PriorSD[!b$levels$Fixed & b$levels$Facet != "Person"],
      rep(1.3, sum(!b$levels$Fixed & b$levels$Facet != "Person")), tolerance = 1e-12)
    expect_equal(b$step_levels$PriorSD[!b$step_levels$Fixed],
      rep(.7, sum(!b$step_levels$Fixed)), tolerance = 1e-12)
    theta <- seq(-2, 2, length.out = s$J)
    slope_free <- if (model == "GPCM") drop(crossprod(s$slope_basis, params$log_slopes)) else numeric()
    a <- exp(drop(s$slope_X %*% slope_free))
    eta <- theta[s$person] + as.vector(s$eta_offset) + drop(s$X %*% free)
    logits <- (outer(eta, 0:(s$C - 1L)) -
      t(apply(cbind(0, steps[s$owner, , drop = FALSE]), 1L, cumsum))) * a
    p <- exp(logits - apply(logits, 1L, max)); p <- p / rowSums(p)
    idx <- mfrmr:::build_indices(f$prep, f$config$step_facet,
      slope_facet = if (model == "GPCM") f$config$slope_facet else NULL,
      gpcm_spec = if (model == "GPCM") f$config$gpcm_spec else NULL)
    eta_native <- mfrmr:::compute_eta(idx, params, f$config, theta_override = theta)
    expected <- mfrmr:::compute_response_probability_bundle(f$config, idx, params, eta_native)$probs
    expect_equal(unname(p), unname(expected), tolerance = 1e-12, ignore_attr = TRUE)
    if (model == "GPCM") {
      expect_equal(sum(params$log_slopes), 0, tolerance = 1e-12)
      expect_equal(unname(sqrt(rowSums(sweep(s$slope_basis, 2L, s$log_slope_sd, `*`)^2))),
        rep(.4, s$A), tolerance = 1e-12)
    }
  }
})

test_that("new export routes remain self-contained and validate before writing", {
  for (anchored in c(FALSE, TRUE)) for (format in c("stan", "brms")) {
    f <- affine_export_fixture(anchored = anchored)
    folder <- tempfile(); on.exit(unlink(folder, recursive = TRUE), add = TRUE)
    out <- export_mfrm_model(f, folder, format, acknowledge_sensitive = TRUE)
    expect_true(all(file.exists(out$files)))
    expect_true(all(c("steps.csv", "slope-levels.csv") %in% names(out$files)))
    b <- readRDS(out$files[["model-data.rds"]])
    expect_identical(b$constraints$facets, f$config$facet_specs)
    if (format == "brms") expect_silent(parse(out$files[["model.R"]]))
  }
  f <- affine_export_fixture()
  f$config$facet_specs[[1]]$anchors[1] <- Inf
  folder <- tempfile()
  expect_error(export_mfrm_model(f, folder), "finite|Inconsistent")
  expect_false(dir.exists(folder))
})

test_that("default GPCM retains its population in Stan and rejects unsupported brms", {
  fit <- readRDS(test_path("fixtures", "anchored-gpcm.rds"))$shared$fit
  before <- unserialize(serialize(fit, NULL))
  folder <- tempfile(); on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  out <- export_mfrm_model(fit, folder, "stan", acknowledge_sensitive = TRUE)
  b <- readRDS(out$files[["model-data.rds"]])
  expect_identical(b$population$mode, "estimated_normal")
  expect_identical(b$population$person_levels, fit$config$person_levels)
  expect_identical(b$constraints$facets, fit$config$facet_specs)
  expect_equal(b$levels$PriorSD[b$levels$Facet == "Person"],
    rep(sqrt(1 + exp(.5)), b$stan_data$J))
  expect_identical(b$schema_version, 2L)
  expect_equal(as.vector(b$stan_data$population_mean_prior), c(0, 1))
  expect_equal(as.vector(b$stan_data$population_log_sd_prior), c(0, .5))
  folder2 <- tempfile()
  expect_error(export_mfrm_model(fit, folder2, "brms"),
    "supports format = 'stan' only.*Keep the original fit")
  expect_false(dir.exists(folder2))
  expect_identical(fit, before)
})

test_that("population export validates priors, design and Person alignment before writing", {
  fit <- readRDS(test_path("fixtures", "anchored-gpcm.rds"))$shared$fit
  for (prior in list(NULL, list(mean = c(0, 1)), list(mean = c(0, -1), log_sd = c(0, .5)),
                     list(mean = c(0, 1), log_sd = c(NA, 1)),
                     list(mean = c(0, 1), other = c(0, 1)))) {
    path <- tempfile()
    expect_error(export_mfrm_model(fit, path, population_prior = prior), "population_prior")
    expect_false(dir.exists(path))
  }
  for (modify in list(
    function(f) { f$config$population_spec$fixed_sd <- 2; f },
    function(f) { f$config$population_spec$design_matrix[, 1] <- 2; f },
    function(f) { f$config$population_spec$person_lookup[1] <- 2L; f })) {
    path <- tempfile()
    expect_error(export_mfrm_model(modify(fit), path), "intercept-only|Person identities")
    expect_false(dir.exists(path))
  }
  # Person-table order is incidental; the explicit lookup preserves identity.
  b <- mfrmr:::mfrm_model_export_data(fit, 1, 1)
  perm <- rev(seq_len(nrow(fit$population$person_table)))
  fit$population$person_table <- fit$population$person_table[perm, , drop = FALSE]
  fit$config$population_spec$design_matrix <- fit$config$population_spec$design_matrix[perm, , drop = FALSE]
  fit$config$population_spec$person_lookup <- match(fit$config$population_spec$person_lookup, perm)
  expect_identical(mfrmr:::mfrm_model_export_data(fit, 1, 1), b)
  custom <- mfrmr:::mfrm_model_export_data(fit, 1, 1,
    population_prior = list(mean = c(.3, .8), log_sd = c(-.2, .4)))
  expect_identical(custom$priors$population, list(mean = c(.3, .8), log_sd = c(-.2, .4)))
  expect_equal(custom$levels$PriorSD[custom$levels$Facet == "Person"],
    rep(sqrt(.8^2 + exp(-.4 + 2 * .4^2)), custom$stan_data$J))
})

test_that("fully fixed blocks and binary steps do not introduce sampled nuisance coordinates", {
  f <- affine_export_fixture("PCM", anchored = FALSE)
  for (facet in f$config$facet_names) {
    levels <- f$config$facet_levels[[facet]]
    f$config$facet_specs[[facet]] <- mfrmr:::build_facet_constraint(levels,
      anchors = stats::setNames(rep(0, length(levels)), levels), centering_scope = "facet")
  }
  f$config$n_cat <- 2L
  f$config$rating_max <- f$config$rating_min + 1L
  f$config$score_map <- f$config$score_map[1:2, ]
  f$prep$data$score_k <- f$prep$data$score_k %% 2L
  f$config$step_specs <- mfrmr:::mfrmr_default_step_specs(f$config)
  b <- mfrmr:::mfrm_model_export_data(f, 1, 1)
  expect_equal(c(b$stan_data$P, b$stan_data$T, b$stan_data$S), c(0, 0, 0))
  expect_true(all(b$levels$Fixed[b$levels$Facet != "Person"]))
  expect_true(all(b$step_levels$Fixed))
})

test_that("owners with different numbers of free steps retain their own constraints", {
  f <- affine_export_fixture("PCM", anchored = FALSE)
  for (i in seq_along(f$config$step_specs)) {
    old <- f$config$step_specs[[i]]
    anchors <- switch((i - 1L) %% 3L + 1L, NULL, c(`1` = .2), c(`1` = .2, `2` = -.1))
    f$config$step_specs[[i]] <- mfrmr:::build_step_constraint(3L, anchors, .3, old$scope)
  }
  b <- mfrmr:::mfrm_model_export_data(f, 1, .6)
  d <- b$stan_data
  z <- cos(seq_len(d$T))
  steps <- as.vector(d$step_offset) + drop(d$step_basis %*% z)
  for (i in seq_along(f$config$step_specs)) {
    spec <- f$config$step_specs[[i]]
    rows <- (i - 1L) * 3L + 1:3
    expect_equal(sum(steps[rows]), .3, tolerance = 1e-12)
    fixed <- !is.na(spec$anchors)
    expect_equal(steps[rows][fixed], unname(spec$anchors[fixed]), tolerance = 1e-12)
    expect_equal(sum(b$step_levels$Fixed[rows]), if (spec$n_params == 0) 3L else sum(fixed))
  }
  expect_equal(d$T, sum(vapply(f$config$step_specs, `[[`, integer(1), "n_params")))
})

test_that("uncentered anchored groups and singleton slopes have the declared priors", {
  spec <- mfrmr:::build_facet_constraint(letters[1:5], c(a = .5),
    c(b = "g", c = "g"), c(g = -.2), centered = FALSE)
  map <- mfrmr:::mfrm_export_affine_map(spec, .8)
  expect_equal(map$offset, c(.5, -.2, -.2, 0, 0))
  expect_equal(sqrt(rowSums(sweep(map$basis, 2L, map$sd, `*`)^2)), c(0, .8, .8, .8, .8))
  native <- mfrmr:::expand_facet_with_constraints(c(.7, -.1, .2), spec)
  free <- drop(crossprod(map$basis, native - map$offset))
  expect_equal(map$offset + drop(map$basis %*% free), unname(native))
  expect_equal(mfrmr:::mfrm_export_basis(1L), matrix(numeric(), 1L, 0L))
})
