# --------------------------------------------------------------------------
# test-estimation-core.R
# Tests the likelihood and estimation core of the mfrmr package.
# --------------------------------------------------------------------------

# ---- 2.1  RSM category probabilities sum to 1 and are non-negative --------

test_that("category_prob_rsm rows sum to 1 and are non-negative", {
  step_cum <- c(0, -0.5, 0.3, 0.7)
  etas <- c(-3, -1, 0, 1, 3)
  probs <- mfrmr:::category_prob_rsm(etas, step_cum)
  expect_equal(rowSums(probs), rep(1, length(etas)), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- 2.2  PCM category probabilities sum to 1 and are non-negative --------

test_that("category_prob_pcm rows sum to 1 and are non-negative", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  etas <- c(-1, 0, 1, 0.5)
  crit_idx <- c(1L, 2L, 1L, 2L)
  probs <- mfrmr:::category_prob_pcm(etas, step_cum_mat, crit_idx)
  expect_equal(rowSums(probs), rep(1, length(etas)), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- 2.2c  GPCM category probabilities sum to 1 and are non-negative -------

test_that("category_prob_gpcm rows sum to 1 and are non-negative", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  etas <- c(-1, 0, 1, 0.5)
  crit_idx <- c(1L, 2L, 1L, 2L)
  slopes <- c(0.8, 1.35)
  probs <- mfrmr:::category_prob_gpcm(
    eta = etas,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx,
    slopes = slopes
  )
  expect_equal(rowSums(probs), rep(1, length(etas)), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- 2.3a  Monotonicity -- higher eta -> higher expected score (RSM) ------

test_that("higher eta produces monotonically higher expected score under RSM", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  etas <- seq(-3, 3, by = 0.5)
  probs <- mfrmr:::category_prob_rsm(etas, step_cum)
  k_vals <- 0:(ncol(probs) - 1)
  expected <- as.vector(probs %*% k_vals)
  diffs <- diff(expected)
  expect_true(all(diffs >= -1e-10),
              info = "Expected score must be non-decreasing in eta")
})

# ---- 2.3b  Monotonicity -- higher eta -> higher expected score (PCM) ------

test_that("higher eta produces monotonically higher expected score under PCM", {
  step_cum_mat <- matrix(c(0, 0.3, 0.8, 1.5,
                           0, 0.1, 0.6, 1.2), nrow = 2, byrow = TRUE)
  etas <- seq(-3, 3, by = 0.5)
  # All observations belong to criterion 1
  crit_idx <- rep(1L, length(etas))
  probs <- mfrmr:::category_prob_pcm(etas, step_cum_mat, crit_idx)
  k_vals <- 0:(ncol(probs) - 1)
  expected <- as.vector(probs %*% k_vals)
  diffs <- diff(expected)
  expect_true(all(diffs >= -1e-10),
              info = "Expected score must be non-decreasing in eta (PCM)")
})

# ---- 2.3c  Monotonicity -- higher eta -> higher expected score (GPCM) ------

test_that("higher eta produces monotonically higher expected score under GPCM", {
  step_cum_mat <- matrix(c(0, 0.3, 0.8, 1.5,
                           0, 0.1, 0.6, 1.2), nrow = 2, byrow = TRUE)
  etas <- seq(-3, 3, by = 0.5)
  crit_idx <- rep(1L, length(etas))
  probs <- mfrmr:::category_prob_gpcm(
    eta = etas,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx,
    slopes = c(1.4, 0.9)
  )
  k_vals <- 0:(ncol(probs) - 1)
  expected <- as.vector(probs %*% k_vals)
  diffs <- diff(expected)
  expect_true(all(diffs >= -1e-10),
              info = "Expected score must be non-decreasing in eta (GPCM)")
})

# ---- 2.4  Numerical stability for extreme eta values ---------------------

test_that("extreme eta values do not produce NaN or Inf", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  extreme_etas <- c(-50, -10, 0, 10, 50)
  probs <- mfrmr:::category_prob_rsm(extreme_etas, step_cum)
  expect_true(all(is.finite(probs)))
  expect_equal(rowSums(probs), rep(1, length(extreme_etas)), tolerance = 1e-10)
})

# ---- 2.5a  RSM log-likelihood is finite and non-positive -----------------

test_that("RSM log-likelihood is finite and non-positive", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  score_k <- c(0L, 1L, 2L, 3L, 1L)
  eta <- c(-1, 0, 0.5, 1, -0.5)
  ll <- mfrmr:::loglik_rsm(eta, score_k, step_cum)
  expect_true(is.finite(ll))
  expect_lte(ll, 0)
})

# ---- 2.5b  PCM log-likelihood is finite and non-positive -----------------

test_that("PCM log-likelihood is finite and non-positive", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  score_k <- c(0L, 1L, 2L, 3L)
  eta <- c(-1, 0, 0.5, 1)
  crit_idx <- c(1L, 2L, 1L, 2L)
  ll <- mfrmr:::loglik_pcm(eta, score_k, step_cum_mat, crit_idx)
  expect_true(is.finite(ll))
  expect_lte(ll, 0)
})

# ---- 2.5c  GPCM log-likelihood is finite and non-positive ------------------

test_that("GPCM log-likelihood is finite and non-positive", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  score_k <- c(0L, 1L, 2L, 3L)
  eta <- c(-1, 0, 0.5, 1)
  crit_idx <- c(1L, 2L, 1L, 2L)
  ll <- mfrmr:::loglik_gpcm(
    eta = eta,
    score_k = score_k,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx,
    slopes = c(0.75, 1.30)
  )
  expect_true(is.finite(ll))
  expect_lte(ll, 0)
})

# ---- 2.5d  GPCM reduces exactly to PCM when all slopes equal 1 -------------

test_that("GPCM probabilities and log-likelihood reduce exactly to PCM", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  eta <- c(-1.2, -0.4, 0.2, 0.9, 1.4)
  score_k <- c(0L, 1L, 2L, 3L, 1L)
  crit_idx <- c(1L, 2L, 1L, 2L, 1L)
  slopes <- rep(1, nrow(step_cum_mat))

  probs_pcm <- mfrmr:::category_prob_pcm(
    eta = eta,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx
  )
  probs_gpcm <- mfrmr:::category_prob_gpcm(
    eta = eta,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx,
    slopes = slopes
  )
  expect_equal(probs_gpcm, probs_pcm, tolerance = 1e-12)

  ll_pcm <- mfrmr:::loglik_pcm(
    eta = eta,
    score_k = score_k,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx
  )
  ll_gpcm <- mfrmr:::loglik_gpcm(
    eta = eta,
    score_k = score_k,
    step_cum_mat = step_cum_mat,
    criterion_idx = crit_idx,
    slopes = slopes
  )
  expect_equal(ll_gpcm, ll_pcm, tolerance = 1e-12)
})

# ---- 2.6  JML vs MML consistency ----------------------------------------

test_that("JML and MML facet estimates are highly correlated", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit_jml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))
  fit_mml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))

  for (facet in c("Rater", "Task", "Criterion")) {
    est_jml <- fit_jml$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    est_mml <- fit_mml$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    r <- cor(est_jml, est_mml)
    expect_gt(r, 0.9,
              label = paste("JML-MML correlation for", facet))
  }
  expect_true(is.finite(fit_jml$summary$LogLik))
  expect_true(is.finite(fit_mml$summary$LogLik))
})

# ---- 2.7a  Convergence -- optim converges with adequate maxit ------------

test_that("optim converges with adequate maxit", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 100, quad_points = 7
  ))
  expect_true(fit$summary$Converged)
  expect_true(is.finite(fit$summary$LogLik))
  expect_lt(fit$summary$LogLik, 0)
})

# ---- 2.7b  Convergence -- maxit=1 produces warning ----------------------

test_that("maxit=1 produces a convergence warning", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_warning(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", model = "RSM", maxit = 1, quad_points = 7
    ),
    "converge|converg",
    ignore.case = TRUE
  )
})

# ---- 2.7c  Convergence -- more iterations give better LogLik -------------

test_that("more iterations produce equal or better log-likelihood", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_short <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 5, quad_points = 7
  ))
  fit_long <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 100, quad_points = 7
  ))
  # More iterations should give equal or better (less negative) log-likelihood
  expect_gte(fit_long$summary$LogLik, fit_short$summary$LogLik - 0.1)
})

# ---- 2.8  Step parameters are finite ------------------------------------

test_that("RSM step estimates are finite", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true(all(is.finite(fit$steps$Estimate)))
  # After centering, steps should sum to approximately 0
  expect_equal(sum(fit$steps$Estimate), 0, tolerance = 1e-6)
})

# ---- 2.8b  Binary scores work as ordered two-category responses ----------

test_that("binary scores are supported as ordered two-category responses", {
  set.seed(42)
  persons <- paste0("P", 1:60)
  items <- paste0("I", 1:5)
  d <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
  theta <- stats::rnorm(length(persons), 0, 1)
  beta <- seq(-1, 1, length.out = length(items))
  eta <- theta[match(d$Person, persons)] - beta[match(d$Item, items)]
  d$Score <- stats::rbinom(nrow(d), 1, stats::plogis(eta))

  fit_rsm <- suppressWarnings(fit_mfrm(
    d, "Person", "Item", "Score",
    method = "JML", model = "RSM", maxit = 100
  ))
  fit_pcm <- suppressWarnings(fit_mfrm(
    d, "Person", "Item", "Score",
    method = "JML", model = "PCM", step_facet = "Item", maxit = 100
  ))

  expect_true(isTRUE(fit_rsm$summary$Converged[1]))
  expect_true(isTRUE(fit_pcm$summary$Converged[1]))
  expect_equal(fit_rsm$summary$Categories[1], 2)
  expect_equal(fit_pcm$summary$Categories[1], 2)
  expect_equal(nrow(fit_rsm$steps), 1)
  expect_equal(nrow(fit_pcm$steps), length(items))
  expect_true(all(is.finite(fit_pcm$steps$Estimate)))
})

test_that("GPCM requires explicit and aligned step/slope facets", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  expect_error(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      model = "GPCM",
      method = "MML"
    ),
    "requires an explicit `step_facet`",
    fixed = TRUE
  )

  expect_error(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      model = "GPCM",
      method = "MML",
      step_facet = "Criterion",
      slope_facet = "Task"
    ),
    "slope_facet == step_facet",
    fixed = TRUE
  )
})

test_that("GPCM core fits return positive discrimination tables", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_jml <- suppressWarnings(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      model = "GPCM",
      step_facet = "Criterion",
      method = "JML",
      maxit = 30
    )
  )
  fit_mml <- suppressWarnings(
    fit_mfrm(
      d,
      "Person",
      c("Rater", "Task", "Criterion"),
      "Score",
      model = "GPCM",
      step_facet = "Criterion",
      method = "MML",
      quad_points = 5,
      maxit = 20
    )
  )

  for (fit in list(fit_jml, fit_mml)) {
    expect_identical(as.character(fit$summary$Model[1]), "GPCM")
    expect_true(is.data.frame(fit$slopes))
    expect_true(all(c("SlopeFacet", "LogEstimate", "Estimate") %in% names(fit$slopes)))
    expect_true(all(is.finite(fit$slopes$Estimate)))
    expect_true(all(fit$slopes$Estimate > 0))
    expect_equal(exp(mean(log(fit$slopes$Estimate))), 1, tolerance = 1e-6)
    expect_true(is.finite(fit$summary$LogLik[1]))
  }
})

test_that("GPCM config scaffold records slope identification metadata", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  prep <- mfrmr:::prepare_mfrm_data(
    d,
    person_col = "Person",
    facet_cols = c("Rater", "Task", "Criterion"),
    score_col = "Score"
  )
  roles <- mfrmr:::resolve_step_and_slope_facets(
    model = "GPCM",
    step_facet = "Criterion",
    slope_facet = NULL,
    facet_names = prep$facet_names
  )
  signs <- mfrmr:::build_facet_signs(prep$facet_names)
  idx <- mfrmr:::build_indices(
    prep,
    step_facet = roles$step_facet,
    slope_facet = roles$slope_facet
  )
  cfg <- mfrmr:::build_estimation_config(
    prep = prep,
    model = "GPCM",
    method = "MML",
    step_facet = roles$step_facet,
    slope_facet = roles$slope_facet,
    weight_col = NULL,
    facet_signs = signs$signs,
    positive_facets = signs$positive_facets,
    noncenter_facet = "Person",
    dummy_facets = character(0),
    anchor_df = NULL,
    group_anchor_df = NULL,
    population = NULL
  )

  expect_equal(cfg$config$step_facet, "Criterion")
  expect_equal(cfg$config$slope_facet, "Criterion")
  expect_equal(cfg$config$gpcm_spec$identification, "sum_to_zero_log_slopes")
  expect_equal(cfg$config$gpcm_spec$scale_reference, "geometric_mean_one")
  expect_equal(cfg$config$gpcm_spec$reduction_reference, "PCM when all slopes equal 1")
  expect_equal(cfg$sizes$log_slopes, length(prep$levels$Criterion) - 1L)
  expect_identical(idx$step_idx, idx$slope_idx)
})

test_that("GPCM slope scaffold expands positive slopes with geometric mean one", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  prep <- mfrmr:::prepare_mfrm_data(
    d,
    person_col = "Person",
    facet_cols = c("Rater", "Task", "Criterion"),
    score_col = "Score"
  )
  signs <- mfrmr:::build_facet_signs(prep$facet_names)
  cfg <- mfrmr:::build_estimation_config(
    prep = prep,
    model = "GPCM",
    method = "MML",
    step_facet = "Criterion",
    slope_facet = "Criterion",
    weight_col = NULL,
    facet_signs = signs$signs,
    positive_facets = signs$positive_facets,
    noncenter_facet = "Person",
    dummy_facets = character(0),
    anchor_df = NULL,
    group_anchor_df = NULL,
    population = NULL
  )

  start <- mfrmr:::build_initial_param_vector(cfg$config, cfg$sizes)
  params0 <- mfrmr:::expand_params(start, cfg$sizes, cfg$config)
  expect_equal(params0$slopes, rep(1, length(prep$levels$Criterion)), tolerance = 1e-12)
  expect_equal(sum(params0$log_slopes), 0, tolerance = 1e-12)

  slices <- mfrmr:::build_param_slices(cfg$sizes)
  perturbed <- start
  perturbed[slices$log_slopes] <- c(0.25, -0.10)
  params1 <- mfrmr:::expand_params(perturbed, cfg$sizes, cfg$config)

  expect_true(all(params1$slopes > 0))
  expect_equal(sum(params1$log_slopes), 0, tolerance = 1e-12)
  expect_equal(exp(mean(log(params1$slopes))), 1, tolerance = 1e-12)
})

# ---- 2.8c  Latent-regression scaffolding validates person-data contract ---

make_latent_regression_fixture <- function(seed = 2718,
                                          n_person = 90,
                                          n_item = 6) {
  set.seed(seed)
  persons <- paste0("P", sprintf("%03d", seq_len(n_person)))
  items <- paste0("I", seq_len(n_item))
  x <- stats::rnorm(n_person)
  theta <- 0.25 + 0.9 * x + stats::rnorm(n_person, sd = 0.6)
  item_beta <- seq(-1.2, 1.2, length.out = n_item)
  dat <- expand.grid(Person = persons, Item = items, stringsAsFactors = FALSE)
  eta <- theta[match(dat$Person, persons)] - item_beta[match(dat$Item, items)]
  dat$Score <- stats::rbinom(nrow(dat), 1, stats::plogis(eta))
  person_tbl <- data.frame(
    Person = persons,
    X = x,
    stringsAsFactors = FALSE
  )
  list(data = dat, person_data = person_tbl)
}

test_that("latent-regression scaffolding requires person_data and MML", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  persons <- unique(d$Person)
  person_tbl <- data.frame(
    Person = persons,
    Grade = seq_along(persons),
    stringsAsFactors = FALSE
  )

  expect_error(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      population_formula = ~ Grade
    ),
    "person_data"
  )

  expect_error(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML",
      population_formula = ~ Grade,
      person_data = person_tbl
    ),
    "requires `method = 'MML'`"
  )
})

test_that("latent-regression scaffolding catches duplicate or missing persons", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  persons <- unique(d$Person)

  dup_tbl <- data.frame(
    Person = c(persons, persons[1]),
    Grade = seq_len(length(persons) + 1L),
    stringsAsFactors = FALSE
  )
  miss_tbl <- data.frame(
    Person = persons[-1],
    Grade = seq_len(length(persons) - 1L),
    stringsAsFactors = FALSE
  )

  expect_error(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      population_formula = ~ Grade,
      person_data = dup_tbl
    ),
    "Duplicate IDs"
  )

  expect_error(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      population_formula = ~ Grade,
      person_data = miss_tbl
    ),
    "missing .* person"
  )
})

test_that("latent-regression scaffolding expands categorical covariates and fits active MML", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  persons <- unique(d$Person)

  cat_tbl <- data.frame(
    Person = persons,
    Group = rep(c("A", "B"), length.out = length(persons)),
    stringsAsFactors = FALSE
  )
  good_tbl <- data.frame(
    Person = persons,
    Grade = seq_along(persons),
    stringsAsFactors = FALSE
  )

  cat_scaffold <- mfrmr:::prepare_mfrm_population_scaffold(
    data = d,
    person = "Person",
    population_formula = ~ Group,
    person_data = cat_tbl
  )
  expect_true(any(startsWith(colnames(cat_scaffold$design_matrix), "Group")))
  expect_identical(cat_scaffold$xlevels$Group, c("A", "B"))
  expect_true("Group" %in% names(cat_scaffold$contrasts))

  fit <- suppressWarnings(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      population_formula = ~ Grade,
      person_data = good_tbl,
      quad_points = 5,
      maxit = 40
    )
  )
  expect_true(isTRUE(fit$population$active))
  expect_identical(fit$population$posterior_basis, "population_model")
  expect_true(all(is.finite(fit$population$coefficients)))
  expect_true(is.finite(fit$population$sigma2))
  expect_gt(fit$population$sigma2, 0)
})

test_that("latent-regression MML fit supports categorical person covariates", {
  fixture <- make_latent_regression_fixture(n_person = 50, n_item = 4)
  person_data <- transform(
    fixture$person_data,
    Group = ifelse(X >= stats::median(X), "high", "low")
  )[, c("Person", "Group")]

  fit <- suppressWarnings(
    fit_mfrm(
      fixture$data,
      "Person", "Item", "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ Group,
      person_data = person_data,
      quad_points = 5,
      maxit = 60
    )
  )

  expect_true(isTRUE(fit$population$active))
  expect_identical(fit$population$posterior_basis, "population_model")
  expect_true(any(startsWith(fit$population$design_columns, "Group")))
  expect_true("Group" %in% names(fit$population$xlevels))
  expect_equal(sort(fit$population$xlevels$Group), c("high", "low"))
  expect_true("Group" %in% names(fit$population$contrasts))
  expect_true(all(is.finite(fit$population$coefficients)))
  expect_true(is.finite(fit$population$sigma2))
  expect_gt(fit$population$sigma2, 0)

  s <- summary(fit)
  expect_true("population_coding" %in% names(s))
  expect_true(is.data.frame(s$population_coding))
  expect_true(all(c("Variable", "LevelCount", "Levels", "Contrast", "EncodedColumns", "CodingNote") %in% names(s$population_coding)))
  expect_identical(as.character(s$population_coding$Variable[1]), "Group")
  expect_equal(s$population_coding$LevelCount[1], 2L)
  expect_match(s$population_coding$Levels[1], "high|low")
  expect_match(s$population_coding$Contrast[1], "contr.treatment", fixed = TRUE)
  expect_match(s$population_coding$EncodedColumns[1], "Group")
  expect_match(paste(s$next_actions, collapse = " "), "post hoc regression", fixed = TRUE)
  printed <- capture.output(print(s))
  expect_true(any(grepl("Population covariate coding", printed, fixed = TRUE)))
  expect_true(any(grepl("Group", printed, fixed = TRUE)))
})

test_that("population scaffold helper can build a complete-case latent-regression design", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  persons <- unique(d$Person)
  person_tbl <- data.frame(
    Person = persons,
    Grade = seq_along(persons),
    Ses = stats::rnorm(length(persons)),
    stringsAsFactors = FALSE
  )
  person_tbl$Ses[1] <- NA_real_

  scaffold <- mfrmr:::prepare_mfrm_population_scaffold(
    data = d,
    person = "Person",
    population_formula = ~ Grade + Ses,
    person_data = person_tbl,
    population_policy = "omit"
  )

  expect_true(isTRUE(scaffold$active))
  expect_identical(scaffold$posterior_basis, "population_model")
  expect_true(is.matrix(scaffold$design_matrix))
  expect_equal(nrow(scaffold$person_table), length(persons) - 1L)
  expect_equal(nrow(scaffold$person_table_replay), length(persons))
  expect_equal(nrow(scaffold$design_matrix), nrow(scaffold$person_table))
  expect_equal(scaffold$included_persons, persons[-1])
  expect_equal(scaffold$omitted_persons, persons[1])
  expect_equal(scaffold$response_rows_omitted, sum(d$Person == persons[1]))
  expect_equal(scaffold$response_rows_retained, nrow(d) - scaffold$response_rows_omitted)
  expect_true(persons[1] %in% scaffold$person_table_replay$Person)
  expect_true(is.na(scaffold$person_table_replay$Ses[scaffold$person_table_replay$Person == persons[1]][1]))
  expect_identical(scaffold$person_table_replay_scope, "observed_person_subset_pre_omit")
})

test_that("population helpers preserve person alignment and transformed quadrature nodes", {
  gh <- mfrmr:::gauss_hermite_normal(5)
  population <- list(
    active = TRUE,
    posterior_basis = "population_model",
    person_id = "Person",
    person_table = data.frame(
      Person = c("P01", "P02"),
      Grade = c(0, 2),
      stringsAsFactors = FALSE
    ),
    design_matrix = cbind(`(Intercept)` = 1, Grade = c(0, 2)),
    design_columns = c("(Intercept)", "Grade"),
    coefficients = c(0.5, 1.0),
    sigma2 = 4
  )

  spec <- mfrmr:::compact_population_spec(population, person_levels = c("P01", "P02"))
  basis <- mfrmr:::resolve_person_quadrature_basis(
    gh,
    population_spec = spec,
    person_ids = 1:2
  )

  expect_true(isTRUE(spec$active))
  expect_equal(spec$person_lookup, c(1L, 2L))
  expect_equal(spec$design_columns, c("(Intercept)", "Grade"))
  expect_equal(basis$mu, c(0.5, 2.5), tolerance = 1e-10)
  expect_equal(basis$sigma, 2, tolerance = 1e-10)
  expect_equal(basis$nodes[1, ], 0.5 + 2 * gh$nodes, tolerance = 1e-10)
  expect_equal(basis$nodes[2, ], 2.5 + 2 * gh$nodes, tolerance = 1e-10)
  expect_equal(exp(basis$log_weights[1, ]), gh$weights, tolerance = 1e-12)
})

test_that("latent-regression MML fit estimates population coefficients and variance", {
  fixture <- make_latent_regression_fixture()
  fit <- suppressWarnings(
    fit_mfrm(
      fixture$data,
      "Person", "Item", "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ X,
      person_data = fixture$person_data,
      quad_points = 7,
      maxit = 80
    )
  )
  s <- summary(fit)

  expect_true(isTRUE(fit$population$active))
  expect_identical(fit$population$posterior_basis, "population_model")
  expect_true(is.finite(fit$population$sigma2))
  expect_gt(fit$population$sigma2, 0)
  expect_true(all(is.finite(fit$population$coefficients)))
  expect_lt(abs(unname(fit$population$coefficients["(Intercept)"]) - 0.25), 0.15)
  expect_gt(unname(fit$population$coefficients["X"]), 0)
  expect_lt(abs(unname(fit$population$coefficients["X"]) - 0.9), 0.35)
  expect_lt(abs(fit$population$sigma2 - 0.36), 0.20)
  expect_true(isTRUE(s$population_overview$PopulationModel[1]))
  expect_true(is.finite(s$population_overview$ResidualVariance[1]))
  expect_true("population_design" %in% names(s))
  expect_true(is.data.frame(s$population_design))
  expect_equal(nrow(s$population_design), ncol(fit$population$design_matrix))
  expect_true(all(c(
    "Column", "IsIntercept", "PersonRows", "NonMissing",
    "Complete", "Mean", "SD", "Min", "Max", "ZeroVariance"
  ) %in% names(s$population_design)))
  expect_identical(as.character(s$population_design$Column), colnames(fit$population$design_matrix))
  expect_true(isTRUE(s$population_design$IsIntercept[s$population_design$Column == "(Intercept)"]))
  expect_true(isTRUE(s$population_design$ZeroVariance[s$population_design$Column == "(Intercept)"]))
  expect_false(isTRUE(s$population_design$ZeroVariance[s$population_design$Column == "X"]))
  expect_true(all(s$population_design$NonMissing == nrow(fit$population$design_matrix)))
  expect_true(all(s$population_design$Complete))
  expect_true("population_coefficients" %in% names(s))
  expect_true(any(s$population_coefficients$Term == "(Intercept)"))
  expect_true(any(s$population_coefficients$Term == "X"))
  expect_true("caveats" %in% names(s))
  expect_false(any(s$caveats$Area == "population_model"))
})

test_that("optimizer diagnostics distinguish converged, reviewable, and hard warnings", {
  converged <- mfrmr:::build_optimizer_diagnostics(
    opt = list(convergence = 0L, counts = stats::setNames(c(18L, 17L), c("function", "gradient")), message = NULL),
    gradient = c(1e-8, -2e-8),
    reltol = 1e-6,
    maxit = 50L,
    optimizer_method = "BFGS"
  )
  expect_identical(converged$ConvergenceStatus, "converged")
  expect_identical(converged$ConvergenceReason, "tolerance_met")
  expect_identical(converged$ConvergenceSeverity, "pass")
  expect_false(isTRUE(converged$ReviewableWarning))
  expect_equal(converged$FunctionEvaluations, 18L)
  expect_equal(converged$GradientEvaluations, 17L)

  reviewable <- mfrmr:::build_optimizer_diagnostics(
    opt = list(convergence = 1L, counts = stats::setNames(c(50L, 49L), c("function", "gradient")), message = "iteration limit"),
    gradient = c(5e-6, -2e-6),
    reltol = 1e-6,
    maxit = 50L,
    optimizer_method = "BFGS"
  )
  expect_identical(reviewable$ConvergenceStatus, "reviewable_warning")
  expect_identical(reviewable$ConvergenceReason, "iteration_limit_small_gradient")
  expect_identical(reviewable$ConvergenceSeverity, "review")
  expect_true(isTRUE(reviewable$ReviewableWarning))

  hard_fail <- mfrmr:::build_optimizer_diagnostics(
    opt = list(convergence = 1L, counts = stats::setNames(c(50L, 49L), c("function", "gradient")), message = "iteration limit"),
    gradient = c(5e-2, -2e-2),
    reltol = 1e-6,
    maxit = 50L,
    optimizer_method = "BFGS"
  )
  expect_identical(hard_fail$ConvergenceStatus, "iteration_limit")
  expect_identical(hard_fail$ConvergenceReason, "iteration_limit_large_gradient")
  expect_identical(hard_fail$ConvergenceSeverity, "fail")
  expect_false(isTRUE(hard_fail$ReviewableWarning))

  em_converged <- mfrmr:::build_optimizer_diagnostics(
    opt = list(convergence = 0L, counts = stats::setNames(c(120L, 95L), c("function", "gradient")), message = "EM converged"),
    gradient = c(1e-4, -8e-5),
    reltol = 1e-6,
    maxit = 50L,
    optimizer_method = "EM",
    convergence_basis = "relative_loglik"
  )
  expect_identical(em_converged$ConvergenceBasis, "relative_loglik")
  expect_identical(em_converged$ConvergenceStatus, "converged")
  expect_identical(em_converged$ConvergenceReason, "relative_loglik_tolerance_met")
  expect_false(isTRUE(em_converged$ReviewableWarning))
})

test_that("MML engine planner falls back only for unsupported combinations", {
  supported <- mfrmr:::resolve_mml_engine_plan(
    method = "MML",
    model = "RSM",
    requested = "hybrid",
    population_active = FALSE
  )
  expect_identical(supported$Used, "hybrid")
  expect_false(isTRUE(supported$Fallback))

  latent_regression <- mfrmr:::resolve_mml_engine_plan(
    method = "MML",
    model = "RSM",
    requested = "em",
    population_active = TRUE
  )
  expect_identical(latent_regression$Used, "direct")
  expect_true(isTRUE(latent_regression$Fallback))

  gpcm <- mfrmr:::resolve_mml_engine_plan(
    method = "MML",
    model = "GPCM",
    requested = "em",
    population_active = FALSE
  )
  expect_identical(gpcm$Used, "direct")
  expect_true(isTRUE(gpcm$Fallback))
})

test_that("EM and hybrid MML engines are wired for RSM/PCM", {
  toy <- mfrmr:::sample_mfrm_data(seed = 444)
  toy <- toy[toy$Person %in% unique(toy$Person)[1:8] &
               toy$Task %in% unique(toy$Task)[1:2] &
               toy$Criterion %in% unique(toy$Criterion)[1:2], , drop = FALSE]

  fit_em <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      model = "RSM",
      quad_points = 5,
      maxit = 8,
      mml_engine = "em"
    )
  )
  expect_identical(fit_em$summary$MMLEngineRequested[1], "em")
  expect_identical(fit_em$summary$MMLEngineUsed[1], "em")
  expect_identical(fit_em$summary$OptimizerMethod[1], "EM")
  expect_identical(fit_em$summary$ConvergenceBasis[1], "relative_loglik")
  expect_true(is.finite(fit_em$summary$EMIterations[1]))
  expect_gte(fit_em$summary$EMIterations[1], 1)

  fit_hybrid <- suppressWarnings(
    fit_mfrm(
      toy,
      "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "MML",
      model = "PCM",
      step_facet = "Criterion",
      quad_points = 5,
      maxit = 8,
      mml_engine = "hybrid"
    )
  )
  expect_identical(fit_hybrid$summary$MMLEngineRequested[1], "hybrid")
  expect_identical(fit_hybrid$summary$MMLEngineUsed[1], "hybrid")
  expect_identical(fit_hybrid$summary$OptimizerMethod[1], "BFGS")
  expect_true(is.finite(fit_hybrid$summary$EMIterations[1]))
  expect_gte(fit_hybrid$summary$EMIterations[1], 1)
})

test_that("EM and hybrid MML engines stay numerically close to direct MML", {
  toy <- mfrmr:::sample_mfrm_data(seed = 777)
  toy <- toy[toy$Person %in% unique(toy$Person)[1:10] &
               toy$Task %in% unique(toy$Task)[1:2] &
               toy$Criterion %in% unique(toy$Criterion)[1:2], , drop = FALSE]

  compare_engines <- function(model, step_facet = NULL) {
    fits <- lapply(c("direct", "em", "hybrid"), function(engine) {
      suppressWarnings(
        fit_mfrm(
          toy,
          "Person", c("Rater", "Task", "Criterion"), "Score",
          method = "MML",
          model = model,
          step_facet = step_facet,
          quad_points = 5,
          maxit = 15,
          mml_engine = engine
        )
      )
    })
    names(fits) <- c("direct", "em", "hybrid")

    direct_loglik <- fits$direct$summary$LogLik[1]
    other_loglik <- vapply(fits[c("em", "hybrid")], function(f) f$summary$LogLik[1], numeric(1))
    expect_true(all(abs(other_loglik - direct_loglik) < 0.01))

    direct_tbl <- fits$direct$facets$others[, c("Facet", "Level", "Estimate")]
    for (engine in c("em", "hybrid")) {
      engine_tbl <- fits[[engine]]$facets$others[, c("Facet", "Level", "Estimate")]
      merged <- merge(
        direct_tbl,
        engine_tbl,
        by = c("Facet", "Level"),
        suffixes = c("_direct", paste0("_", engine)),
        all = FALSE
      )
      diff <- abs(merged$Estimate_direct - merged[[paste0("Estimate_", engine)]])
      expect_true(all(diff < 0.01))
    }
  }

  compare_engines("RSM")
  compare_engines("PCM", step_facet = "Criterion")
})

test_that("latent-regression omit policy excludes incomplete persons from estimation", {
  fixture <- make_latent_regression_fixture()
  fixture$person_data$X[1] <- NA_real_
  omitted_id <- fixture$person_data$Person[1]

  fit <- suppressWarnings(
    fit_mfrm(
      fixture$data,
      "Person", "Item", "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ X,
      person_data = fixture$person_data,
      population_policy = "omit",
      quad_points = 7,
      maxit = 80
    )
  )
  complete_person_data <- fixture$person_data[!is.na(fixture$person_data$X), , drop = FALSE]
  complete_data <- fixture$data[fixture$data$Person %in% complete_person_data$Person, , drop = FALSE]
  complete_fit <- suppressWarnings(
    fit_mfrm(
      complete_data,
      "Person", "Item", "Score",
      method = "MML",
      model = "RSM",
      population_formula = ~ X,
      person_data = complete_person_data,
      quad_points = 7,
      maxit = 80
    )
  )

  expect_true(isTRUE(fit$population$active))
  expect_equal(fit$population$omitted_persons, omitted_id)
  expect_equal(fit$population$response_rows_omitted, length(unique(fixture$data$Item)))
  expect_equal(nrow(fit$facets$person), nrow(fixture$person_data) - 1L)
  expect_false(omitted_id %in% fit$facets$person$Person)
  expect_setequal(fit$facets$person$Person, complete_fit$facets$person$Person)
  expect_equal(fit$population$coefficients, complete_fit$population$coefficients, tolerance = 1e-8)
  expect_equal(fit$population$sigma2, complete_fit$population$sigma2, tolerance = 1e-8)
  expect_true(is.data.frame(fit$population$person_table_replay))
  expect_equal(nrow(fit$population$person_table_replay), nrow(fixture$person_data))
  expect_true(omitted_id %in% fit$population$person_table_replay$Person)
  expect_true(is.na(fit$population$person_table_replay$X[fit$population$person_table_replay$Person == omitted_id][1]))
  expect_identical(fit$population$person_table_replay_scope, "observed_person_subset_pre_omit")

  s <- summary(fit)
  expect_true(any(s$caveats$Condition == "population_complete_case_omission"))
  expect_true(any(grepl("Latent-regression fit omitted", s$key_warnings, fixed = TRUE)))
})

test_that("population caveat collector covers all structured warning conditions", {
  population_overview <- data.frame(
    PopulationModel = TRUE,
    ResidualVariance = 0,
    OmittedPersons = 2L,
    OmittedRows = 8L,
    stringsAsFactors = FALSE
  )
  population_design <- data.frame(
    Column = c("(Intercept)", "X_zero", "X_incomplete"),
    IsIntercept = c(TRUE, FALSE, FALSE),
    ZeroVariance = c(TRUE, TRUE, FALSE),
    Complete = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  caveats <- mfrmr:::collect_mfrm_population_caveats(
    population_overview = population_overview,
    population_design = population_design,
    population_coefficients = data.frame()
  )

  expect_setequal(
    caveats$Condition,
    c(
      "population_complete_case_omission",
      "population_design_zero_variance",
      "population_design_incomplete",
      "population_residual_variance_unstable",
      "population_coefficients_missing"
    )
  )
  expect_true(all(caveats$Area == "population_model"))
})

# ---- 2.9a  EAP person estimates are all finite (MML) --------------------

test_that("all EAP person estimates are finite", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true(all(is.finite(fit$facets$person$Estimate)))
})

# ---- 2.9b  EAP SDs are all positive (MML) -------------------------------

test_that("all EAP person SDs are positive", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true("SD" %in% names(fit$facets$person))
  expect_true(all(fit$facets$person$SD > 0))
})

# ---- 2.10a  ZSTD: mnsq = 1.0 -> ZSTD near 0 ----------------------------

test_that("ZSTD is near zero when mnsq equals 1", {
  z <- mfrmr:::zstd_from_mnsq(1.0, df = 30)
  expect_equal(z, 0, tolerance = 0.3)
})

# ---- 2.10b  ZSTD: mnsq > 1 -> ZSTD > 0 ---------------------------------

test_that("ZSTD is positive when mnsq exceeds 1", {
  z <- mfrmr:::zstd_from_mnsq(1.5, df = 30)
  expect_gt(z, 0)
})

# ---- 2.10c  ZSTD: mnsq < 1 -> ZSTD < 0 ---------------------------------

test_that("ZSTD is negative when mnsq is below 1", {
  z <- mfrmr:::zstd_from_mnsq(0.7, df = 30)
  expect_lt(z, 0)
})

# ---- 2.11a  Gauss-Hermite: weights sum to 1 -----------------------------

test_that("Gauss-Hermite weights sum to 1", {
  for (n in c(5, 11, 21)) {
    gh <- mfrmr:::gauss_hermite_normal(n)
    expect_equal(sum(gh$weights), 1, tolerance = 1e-10,
                 label = paste("GH weights sum for n =", n))
  }
})

# ---- 2.11b  Gauss-Hermite: nodes are symmetric around 0 -----------------

test_that("Gauss-Hermite nodes are symmetric around 0", {
  for (n in c(5, 11, 21)) {
    gh <- mfrmr:::gauss_hermite_normal(n)
    nodes_sorted <- sort(gh$nodes)
    nodes_neg_sorted <- sort(-gh$nodes)
    expect_equal(nodes_sorted, nodes_neg_sorted, tolerance = 1e-10,
                 label = paste("GH node symmetry for n =", n))
  }
})

# ---- 2.11c  Gauss-Hermite: E[X] approx 0, E[X^2] approx 1 for N(0,1) ---

test_that("Gauss-Hermite recovers N(0,1) moments", {
  gh <- mfrmr:::gauss_hermite_normal(21)
  ex <- sum(gh$weights * gh$nodes)
  ex2 <- sum(gh$weights * gh$nodes^2)
  expect_equal(ex, 0, tolerance = 0.01,
               label = "E[X] should be ~0 under standard normal")
  expect_equal(ex2, 1, tolerance = 0.01,
               label = "E[X^2] should be ~1 under standard normal")
})

# ---- 2.12  LogLik non-deterioration across optimization ------------------

test_that("log-likelihood does not deteriorate with more iterations (MML)", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_short <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 5, quad_points = 7
  ))
  fit_long <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 100, quad_points = 7
  ))
  expect_gte(fit_long$summary$LogLik, fit_short$summary$LogLik - 0.1)
})

# ---- Additional: logsumexp stability ------------------------------------

test_that("logsumexp handles large values without overflow", {
  x <- c(1000, 1001, 999)
  result <- mfrmr:::logsumexp(x)
  expect_true(is.finite(result))
  # Should be close to log(exp(1000) + exp(1001) + exp(999)) ~= 1001.41
  expect_equal(result, log(exp(0) + exp(1) + exp(-1)) + 1000, tolerance = 1e-10)
})

test_that("logsumexp handles very negative values", {
  x <- c(-1000, -999, -1001)
  result <- mfrmr:::logsumexp(x)
  expect_true(is.finite(result))
})

# ---- Additional: RSM probability boundary behavior -----------------------

test_that("RSM at extreme positive eta concentrates on highest category", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  probs <- mfrmr:::category_prob_rsm(50, step_cum)
  # At very high eta, highest category should dominate

  expect_gt(probs[1, ncol(probs)], 0.99)
})

test_that("RSM at extreme negative eta concentrates on lowest category", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  probs <- mfrmr:::category_prob_rsm(-50, step_cum)
  # At very low eta, lowest category should dominate
  expect_gt(probs[1, 1], 0.99)
})
