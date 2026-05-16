test_that("two-way interaction expansion has zero marginal sums", {
  spec <- list(n_a = 3L, n_b = 4L)
  free <- seq_len((spec$n_a - 1L) * (spec$n_b - 1L)) / 10

  mat <- mfrmr:::expand_two_way_interaction(free, spec)

  expect_equal(rowSums(mat), rep(0, spec$n_a), tolerance = 1e-12)
  expect_equal(colSums(mat), rep(0, spec$n_b), tolerance = 1e-12)

  grad_expanded <- matrix(seq_len(spec$n_a * spec$n_b),
                          nrow = spec$n_a, ncol = spec$n_b)
  projected <- mfrmr:::project_two_way_interaction_gradient(
    grad_expanded,
    spec
  )
  objective <- function(x) {
    sum(grad_expanded * mfrmr:::expand_two_way_interaction(x, spec))
  }
  eps <- 1e-6
  finite_diff <- vapply(seq_along(free), function(i) {
    up <- free
    down <- free
    up[i] <- up[i] + eps
    down[i] <- down[i] - eps
    (objective(up) - objective(down)) / (2 * eps)
  }, numeric(1))

  expect_equal(projected, finite_diff, tolerance = 1e-6)
})

test_that("fit_mfrm estimates explicit non-person facet interactions", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit <- suppressMessages(suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 0,
    maxit = 100
  )))

  tbl <- interaction_effect_table(fit)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 9L)
  expect_equal(unique(tbl$Interaction), "Rater:Criterion")
  expect_true(all(c("Estimate", "N", "WeightedN", "Sparse") %in% names(tbl)))
  expect_false(any(tbl$Sparse))

  mat <- matrix(tbl$Estimate, nrow = 3L, ncol = 3L)
  expect_equal(rowSums(mat), rep(0, 3L), tolerance = 1e-6)
  expect_equal(colSums(mat), rep(0, 3L), tolerance = 1e-6)
  expect_equal(fit$summary$FacetInteractions[1], 1L)
  expect_equal(fit$summary$InteractionParameters[1], 4L)

  sfit <- summary(fit)
  expect_true("interaction_overview" %in% names(sfit))
  expect_equal(sfit$interaction_overview$Cells[1], 9L)
})

test_that("interaction_effect_table is empty for additive fits", {
  fit <- make_toy_fit(method = "JML", maxit = 25, model = "RSM")

  tbl <- interaction_effect_table(fit)

  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 0L)
})

test_that("facet_interactions validates the current modeling boundary", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  base_call <- function(facet_interactions, model = "RSM") {
    suppressMessages(fit_mfrm(
      d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = model,
      step_facet = "Criterion",
      facet_interactions = facet_interactions,
      min_obs_per_interaction = 0,
      maxit = 5
    ))
  }

  expect_error(base_call("Person:Rater"), "non-person")
  expect_error(base_call("Rater:Criterion:Task"), "two-way")
  expect_error(base_call("Rater:Missing"), "not supplied")
  expect_error(base_call("Rater:Rater"), "distinct facets")
  expect_error(base_call("Rater:Criterion", model = "GPCM"), "RSM")
})

test_that("compare_mfrm recognizes additive-to-interaction nesting", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit_add <- suppressMessages(suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 50
  )))
  fit_int <- suppressMessages(suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 0,
    maxit = 50
  )))

  fit_add$summary$Converged[1] <- TRUE
  fit_int$summary$Converged[1] <- TRUE
  comp <- suppressWarnings(compare_mfrm(
    fit_add,
    fit_int,
    labels = c("Additive", "Interaction"),
    nested = TRUE
  ))

  review <- comp$comparison_basis$nesting_review
  expect_true(isTRUE(review$eligible))
  expect_identical(as.character(review$relation), "facet_interaction_extension")
  expect_identical(review$simpler, "Additive")
  expect_identical(review$complex, "Interaction")
  expect_null(comp$lrt)

  fit_rev <- fit_int
  fit_rev$config$facet_interactions <- "Criterion:Rater"
  comp_same <- suppressWarnings(compare_mfrm(
    fit_int,
    fit_rev,
    labels = c("Forward", "Reverse"),
    nested = TRUE
  ))
  expect_false(isTRUE(comp_same$comparison_basis$nesting_review$eligible))
  expect_identical(
    as.character(comp_same$comparison_basis$nesting_review$relation),
    "same_model"
  )
})

test_that("interaction MML falls back from EM to direct optimization", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit <- suppressMessages(suppressWarnings(fit_mfrm(
    d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    model = "RSM",
    facet_interactions = "Rater:Criterion",
    min_obs_per_interaction = 0,
    quad_points = 3,
    maxit = 2,
    mml_engine = "em"
  )))

  expect_equal(fit$summary$MMLEngineRequested[1], "em")
  expect_equal(fit$summary$MMLEngineUsed[1], "direct")
  expect_match(fit$summary$MMLEngineDetail[1], "falling back to direct")
})
