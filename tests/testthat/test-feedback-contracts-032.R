test_that("GPCM local bias uses fitted slopes, weights and independent cell likelihoods", {
  d <- load_mfrmr_data("example_core")
  d$Weight <- rep(c(.5, 1, 2), length.out = nrow(d))
  for (formula in c("full_predictor", "ability_only")) {
    f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "GPCM",
      step_facet = "Criterion", slope_facet = c("Rater", "Criterion"), weight = "Weight",
      spec = list(gpcm_mml_identification = "fixed_standard_normal", slope_type = formula))
    b <- estimate_bias(f, interaction_facets = c("Rater", "Criterion"))
    st <- .mfrm_fit_state(f)
    native <- mfrm_cpp_fast_bias_pairwise(st$enc, st$theta, st$theta_exclude,
      st$beta_list, st$step_cum, 0L, 1L, TRUE, st$alpha)$cells[[1L]]
    expect_true(all(is.na(b$table$BiasSE)))
    expect_true(all(grepl("nonunit_observation_weights", b$table$InferenceReason)))
    dd <- f$encoding$retained_data
    expected <- predict(f)
    for (i in c(1L, 5L, 10L)) {
      row <- b$table[i, ]
      idx <- which(dd$Rater == row$LevelA & dd$Criterion == row$LevelB)
      w <- dd$Weight[idx]
      expect_equal(row$ExpScore, sum(w * expected[idx]), tolerance = 1e-10)
      theta <- f$person$Estimate[match(dd$Person[idx], f$person$Person)]
      beta <- sum(f$facets$Estimate[f$facets$Level %in% c(row$LevelA, row$LevelB)])
      alpha <- prod(f$slopes$Estimate[f$slopes$Level %in% c(row$LevelA, row$LevelB)])
      cum <- f$step_cumulative[match(row$LevelB, f$encoding$step_owners), ]
      k <- 0:(f$encoding$n_cat - 1L)
      y <- dd$Score[idx] - f$rating_min
      nll <- function(delta) {
        eta <- if (formula == "full_predictor") alpha * (theta - beta - delta) else alpha * theta - beta - delta
        lp <- outer(eta, k) - rep(if (formula == "full_predictor") alpha * cum else cum, each = length(idx))
        lp <- lp - apply(lp, 1L, max)
        -sum(w * (lp[cbind(seq_along(idx), y + 1L)] - log(rowSums(exp(lp)))))
      }
      optimum <- optimize(nll, c(-8, 8), tol = 1e-10)$minimum
      h <- 1e-3
      curvature <- (nll(optimum + h) - 2 * nll(optimum) + nll(optimum - h)) / h^2
      expect_equal(row$BiasSize, optimum, tolerance = 2e-6)
      expect_equal(native$bias_se[i], 1 / sqrt(curvature), tolerance = 2e-6)
    }
    pw <- b$pairwise
    expect_equal(pw$ContrastSE, sqrt(pw$SE1^2 + pw$SE2^2), tolerance = 1e-12)
  }
})

test_that("joint bias changes only the declared interaction and uses checked LRTs", {
  d <- load_mfrmr_data("example_core")
  d$Weight <- rep(c(1, 2), length.out = nrow(d))
  args <- list(data = d, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
    method = "MML", population_formula = ~ 1, weight = "Weight", control = list(quad_points = 41L))
  f <- do.call(fit_mfrm, args)
  expect_warning(b <- estimate_bias(f, interaction_facets = c("Rater", "Criterion"), method = "joint"), "eligible")
  ref <- do.call(fit_mfrm, c(args, list(facet_interactions = list(c("Rater", "Criterion")))))
  expect_equal(b$joint_fit$loglik, ref$loglik, tolerance = 1e-8)
  expect_null(b$lrt)
  expect_warning(comparison <- anova(f, ref), "eligible")
  expect_null(comparison$lrt)
  expect_true(all(is.na(b$table$pValue)))
  expect_identical(b$joint_fit$control, f$control)
  for (key in c("weight", "pop_X_mat", "pop_mode", "fixed_sigma", "custom_step_B", "step_cat_min", "step_cat_max"))
    expect_identical(b$joint_fit$encoding[[key]], f$encoding[[key]])
  expect_error(estimate_bias(b$joint_fit, interaction_facets = c("Criterion", "Rater"), method = "joint"), "already")
  f$converged <- FALSE
  expect_warning(bad <- estimate_bias(f, method = "joint"), "eligible")
  expect_null(bad$lrt)
  expect_true(all(is.na(bad$table$pValue)))
  expect_true(all(as.character(plot(bad)$data$Status) == "Unavailable"))
})

test_that("DFF retains weights and population covariates in both likelihood fits", {
  d <- load_mfrmr_data("example_core")
  d$Weight <- rep(c(1, 2), length.out = nrow(d))
  pd <- unique(d[c("Person", "Group")])
  pd$X <- seq_len(nrow(pd)) / nrow(pd)
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "MML",
    weight = "Weight", population_formula = ~ X, person_data = pd, control = list(quad_points = 41L))
  z <- analyze_dff(f, facet = "Rater", group = "Group")
  base <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "MML",
    weight = "Weight", population_formula = ~ X + Group, person_data = pd, control = list(quad_points = 41L))
  ds <- d
  focal <- z$focal; ref <- z$reference
  ds$Rater <- as.character(ds$Rater)
  ds$Rater[d$Rater == "R01"] <- paste0("R01_", d$Group[d$Rater == "R01"])
  split <- fit_mfrm(ds, "Person", c("Rater", "Criterion"), "Score", method = "MML",
    weight = "Weight", population_formula = ~ X + Group, person_data = pd, control = list(quad_points = 41L))
  expect_equal(z$loglik_base, base$loglik, tolerance = 1e-8)
  levels <- split$facets$Level
  focal <- split$facets$Estimate[levels == paste0("R01_", z$focal)]
  reference <- split$facets$Estimate[levels == paste0("R01_", z$reference)]
  expect_equal(z$table$Contrast[z$table$Level == "R01"], focal - reference, tolerance = 1e-7)
  expect_true(all(is.na(z$table$Statistic)))
  expect_true(all(grepl("nonunit_observation_weights", z$table$InferenceReason)))
  shuffled <- d[rev(seq_len(nrow(d))), ]
  expect_equal(analyze_dff(f, facet = "Rater", group = "Group", data = shuffled)$table, z$table, tolerance = 1e-10)
  shuffled$Score[1] <- (shuffled$Score[1] %% 4) + 1
  expect_error(analyze_dff(f, facet = "Rater", group = "Group", data = shuffled), "match")
  g <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "GPCM", step_facet = "Criterion")
  expect_error(analyze_dff(g, facet = "Rater", group = "Group"), "cannot be dropped")
})

test_that("DFF checks both refits and preserves unavailable screens", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score")
  original <- mfrmr:::mfrm_cpp_fast_batch_warm_list
  local_mocked_bindings(mfrm_cpp_fast_batch_warm_list = function(...) {
    out <- original(...)
    out[[1L]]$converged <- FALSE
    out
  }, .package = "mfrmr")
  z <- analyze_dff(f, facet = "Rater", group = "Group")
  expect_true(is.na(z$table$pValue[1L]))
  expect_true(is.na(z$table$Flagged[1L]))
  expect_match(z$table$InferenceReason[1L], "optimizer_not_converged")
  expect_true(all(is.finite(z$table$pValue[-1L])))
  expect_equal(summary(z)$overview$ElementsUnavailable, 1)
  expect_true("Unavailable" %in% as.character(plot(z)$data$Status))
})

test_that("joint refits retain anchored and constrained encodings", {
  d <- load_mfrmr_data("example_core")
  ga <- data.frame(Facet = "Rater", Level = "R01", Group = "fixed", GroupValue = 0)
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", group_anchors = ga)
  b <- estimate_bias(f, interaction_facets = c("Rater", "Criterion"), method = "joint")
  ref <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", group_anchors = ga,
    facet_interactions = list(c("Rater", "Criterion")))
  expect_equal(b$joint_fit$loglik, ref$loglik, tolerance = 1e-8)
  expect_equal(b$joint_fit$encoding$anchor_fixed[1:2], f$encoding$anchor_fixed)
  expect_equal(b$joint_fit$encoding$anchor_val[1:2], f$encoding$anchor_val)
  expect_error(analyze_dff(f, facet = "Rater", group = "Group"), "cannot be dropped")
})

test_that("Raju-Welch DFF uses the fixed GPCM calibration with interactions", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "GPCM",
    step_facet = "Criterion", facet_interactions = list(c("Rater", "Criterion")))
  local_mocked_bindings(mfrm_cpp_fast_fit_spec = function(...) stop("unexpected refit"), .package = "mfrmr")
  z <- analyze_dff(f, facet = "Rater", group = "Group", method = "raju_welch")
  expect_equal(nrow(z$table), length(unique(d$Rater)))
  expect_true(all(is.finite(z$table$Contrast)))
  expect_true(all(z$table$Level %in% d$Rater))
})

test_that("structural boundaries survive diagnostics, DFF and saved-object displays", {
  d <- expand.grid(Person = paste0("P", 1:20), Rater = paste0("R", 1:3), stringsAsFactors = FALSE)
  d$Score <- ifelse(d$Rater == "R1", 3, rep(0:2, length.out = nrow(d)))
  d$Group <- ifelse(as.integer(sub("P", "", d$Person)) %% 2L, "A", "B")
  f <- suppressWarnings(fit_mfrm(d, "Person", "Rater", "Score"))
  z <- suppressWarnings(analyze_dff(f, facet = "Rater", group = "Group"))
  expect_true(all(is.na(z$table$pValue)))
  expect_true(all(is.na(z$table$Flagged)))
  expect_true(is.na(summary(z)$overview$Flagged))
  expect_equal(summary(z)$overview$ElementsTested, 0)
  expect_equal(summary(z)$overview$ElementsUnavailable, 3)
  dg <- diagnose_mfrm(f)
  expect_false(summary(dg)$readiness$ready)
  expect_true(all(is.na(summary(dg)$reliability$pValue)))
  expect_match(paste(capture.output(print(summary(dg))), collapse = " "), "structural_likelihood_boundary")
  expect_match(paste(attr(plot(dg), "mfrm_notes"), collapse = " "), "Intervals unavailable")
  old <- dg; old$readiness <- NULL
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(list(diagnostics = old, dff = z), path)
  saved <- readRDS(path)
  expect_false(summary(saved$diagnostics)$readiness$ready)
  expect_true(all(as.character(plot(saved$dff)$data$Status) == "Unavailable"))
})

test_that("G and Phi distinguish zero total variance, valid zero and negative components", {
  d <- load_mfrmr_data("example_core")
  d$Score <- d$Other <- 2
  g <- mfrm_g_study(d, scores = c("Score", "Other"), person = "Person", facets = c("Rater", "Criterion"), method = "urgenova")
  expect_true(all(is.na(g$g_study_summary[c("G", "Phi")])))
  expect_true(all(is.na(g$univariate[c("G", "Phi")])))
  expect_true(all(g$d_study$GStatus == "Zero total variance"))
  expect_true(all(is.na(g$d_study$MeetsG80)))
  expect_equal(mfrm_d_study(g, g$d_study[g$facets]), mfrm_d_study(g))
  # Controlled components isolate algebraic states without relying on rounded estimates.
  g$variance_components$Estimate[] <- 0
  g$variance_components$Estimate[g$variance_components$Component == "Person:Rater:Criterion,e"] <- 1
  ds <- mfrm_d_study(g)
  expect_true(all(ds$G == 0 & ds$Phi == 0))
  expect_true(all(ds$GStatus == "Available"))
  g$variance_components$Estimate[g$variance_components$Component == "Person"] <- 1
  g$variance_components$Estimate[g$variance_components$Component == "Person:Rater:Criterion,e"] <- -1
  ds <- mfrm_d_study(g)
  expect_true(all(is.na(ds$G) & is.na(ds$RelSEM)))
  expect_true(all(ds$GStatus == "Negative projected error variance"))
  old <- mfrm_g_study(d, scores = "Score", person = "Person", facets = c("Rater", "Criterion"), method = "urgenova")
  old$g_study_summary$G <- old$g_study_summary$Phi <- 0
  old$d_study$G <- old$d_study$Phi <- 0
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(old, path); old <- readRDS(path)
  expect_true(is.na(summary(old)$g_study_summary$G))
  expect_true(all(is.na(plot(old)$data$Coef)))
  expect_match(paste(attr(plot(old), "mfrm_notes"), collapse = " "), "Zero total variance")
  expect_true(all(is.na(confint(old, component = "coefficients"))))
})
