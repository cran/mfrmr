test_that("slope intervals use log-scale covariance and remain positive", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "GPCM",
                step_facet = "Criterion", slope_facet = "Criterion",
                spec = list(gpcm_mml_identification = "fixed_standard_normal"))
  idx <- grep("^LogSlope:", names(f$xi))
  for (kind in c("model", "sandwich")) {
    V <- vcov(f, se_type = kind)[idx, idx, drop = FALSE]
    # Last log-slope is the negative sum of the free coordinates.
    B <- rbind(diag(length(idx)), -1)
    se <- sqrt(diag(B %*% V %*% t(B)))
    expected <- exp(cbind(log(f$slopes$Estimate) - qnorm(.975) * se,
                          log(f$slopes$Estimate) + qnorm(.975) * se))
    ci <- confint(f, component = "slopes", se_type = kind)
    expect_equal(unname(ci), unname(expected), tolerance = 1e-10)
    expect_true(all(ci > 0))
    expect_equal(unname(confint(f, component = "slopes", se_type = kind, parm = 2)),
                 unname(expected[2, , drop = FALSE]))
  }
  p <- plot(f, type = "slopes")
  expect_equal(unname(as.matrix(p$data[c("Low", "High")])),
               unname(confint(f, component = "slopes")))
  f$converged <- FALSE
  for (type in c("slopes", "persons", "raters")) {
    p <- plot(f, type = type)
    expect_true(all(is.na(p$data$Low) & is.na(p$data$High)))
    expect_match(paste(attr(p, "mfrm_notes"), collapse = " "), "Intervals unavailable.*optimizer_not_converged")
    expect_gt(nrow(p$data), 0L)
  }
})

test_that("G-Study covariance agrees with Gaussian quadratic forms on unbalanced and replicated data", {
  set.seed(319)
  d <- expand.grid(Person = paste0("P", 1:8), Rater = paste0("R", 1:4))
  d$Score <- rep(rnorm(8, sd = 2), 4) + rep(rnorm(4), each = 8) + rnorm(32, sd = .5)
  unbalanced <- d[-c(1, 2, 3, 9, 10, 17, 25, 26, 27, 28), ]
  replicated <- rbind(unbalanced, transform(unbalanced[1:6, ], Score = Score + .2))
  for (x in list(d, unbalanced, replicated, replicated[sample(nrow(replicated)), ])) {
    g <- mfrm_g_study(x, scores = "Score", person = "Person", facets = "Rater", method = "urgenova")
    C <- lapply(c("Person", "Rater"), function(nm) outer(x[[nm]], x[[nm]], "==") * 1)
    C[[3L]] <- C[[1L]] * C[[2L]]
    if (anyDuplicated(x[c("Person", "Rater")])) C[[4L]] <- diag(nrow(x))
    A <- lapply(C, function(z) z / rowSums(z) - 1 / nrow(x))
    K <- sapply(C, function(z) vapply(A, function(a) sum(a * z), numeric(1)))
    inverse <- solve(K)
    raw <- drop(inverse %*% vapply(A, function(a) drop(crossprod(x$Score, a %*% x$Score)), numeric(1)))
    V <- Reduce(`+`, Map(function(z, v) z * v, C, pmax(0, raw)))
    Q <- outer(seq_along(A), seq_along(A), Vectorize(function(i, j) 2 * sum((A[[i]] %*% V) * t(A[[j]] %*% V))))
    expected <- inverse %*% Q %*% t(inverse)
    expect_equal(g$variance_components$RawEstimate, raw, tolerance = 1e-10)
    expect_equal(unname(vcov(g)), unname(expected), tolerance = 1e-10)
    expect_equal(g$variance_components$SE, sqrt(diag(expected)), tolerance = 1e-10)
  }
})

test_that("G-Study rejects unidentified designs, missing-error input and invalid D-Study counts", {
  d <- data.frame(Person = 1:8, Rater = 1:8, Score = 1:8)
  expect_error(mfrm_g_study(d, "Score", "Person", facets = "Rater", method = "urgenova"), "not identifiable")
  d <- expand.grid(Person = 1:8, Rater = 1:4)
  d$Score <- sin(seq_len(nrow(d))) + d$Person
  missing <- d; missing$Score[1] <- NA_real_
  expect_error(mfrm_g_study(missing, "Score", "Person", facets = "Rater", method = "urgenova", missing = "error"), "missing|Missing")
  g <- mfrm_g_study(d, "Score", "Person", facets = "Rater", method = "urgenova")
  for (n in c(0, -1, 1.5, Inf, NA_real_)) {
    expect_error(mfrm_d_study(g, data.frame(Rater = n)), "positive integer")
    expect_error(mfrm_g_study(d, "Score", "Person", facets = "Rater", method = "urgenova", d_study_grid = data.frame(Rater = n)), "positive integer")
  }
})

test_that("missing codes and documented category support reach the unified estimator", {
  d <- load_mfrmr_data("example_core")
  x <- d; x$Score[1] <- -1
  a <- fit_mfrm(x, "Person", c("Rater", "Criterion"), "Score", missing_codes = -1, rating_min = 1, rating_max = 4)
  b <- fit_mfrm(d[-1, ], "Person", c("Rater", "Criterion"), "Score")
  expect_equal(a$n_obs, nrow(d) - 1L)
  expect_equal(a$xi, b$xi)
  x <- d; x$Rater[1] <- NA
  expect_equal(fit_mfrm(x, "Person", c("Rater", "Criterion"), "Score")$n_obs, nrow(d) - 1L)
  support <- data.frame(Level = sort(unique(d$Criterion)), Min = 1, Max = 4)
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion", category_support = support)
  expect_equal(f$encoding$step_cat_max, rep(3L, 4))
})

test_that("nested comparisons require matching responses, weights and usable likelihood fits", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score")
  p <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion")
  expect_identical(anova(f, p)$comparison_basis$lrt_status, "computed")
  set.seed(42)
  x <- d; x$Score <- sample(x$Score)
  other <- fit_mfrm(x, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion")
  stopped <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion", maxit = 1)
  weighted <- p; weighted$encoding$weight[1] <- 2
  for (bad in list(other, stopped, weighted)) {
    expect_warning(z <- anova(f, bad), "identical observations")
    expect_identical(z$comparison_basis$lrt_status, "ineligible")
    expect_null(z$lrt)
  }
  worse <- p; worse$loglik <- f$loglik - 1
  expect_warning(z <- anova(f, worse), "lower maximized likelihood")
  expect_null(z$lrt)
})

test_that("structural likelihood boundaries propagate into classification and saved fits", {
  d <- expand.grid(Person = paste0("P", 1:20), Rater = paste0("R", 1:3))
  d$Score <- ifelse(d$Rater == "R1", 3, rep(0:2, length.out = nrow(d)))
  expect_warning(f <- fit_mfrm(d, "Person", "Rater", "Score", rating_min = 0, rating_max = 3), "No finite structural MML")
  expect_true(f$boundary_audit$detected)
  expect_false(f$summary$InferenceReady)
  expect_error(classification_table(f), "review|ready|readiness")
  expect_error(confint(f), "not ready")
  expect_false(classification_table(f, readiness_policy = "review")$source$ready)
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(f, path)
  expect_identical(readRDS(path)$boundary_audit, f$boundary_audit)
  old <- f; old$boundary_audit <- old$readiness <- NULL
  expect_false(summary(old)$readiness$ready)
  expect_true("structural_likelihood_boundary" %in% summary(old)$readiness$reason_codes)
})

test_that("EAP reliability and WLE diagnostics use their stated score definitions", {
  d <- load_mfrmr_data("example_operational")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score")
  r <- diagnose_mfrm(f)$reliability
  v <- var(f$person$EAP); s <- mean(f$person$PosteriorSD^2)
  expect_equal(r$ModelReliability[r$Facet == "Person"], v / (v + s))
  expect_equal(f$population$eap_reliability, v / (v + s))
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "PCM", step_facet = "Criterion", person_score = "WLE")
  pr <- predict(f, type = "probs")
  expected <- drop(pr %*% 1:4)
  variance <- drop(pr %*% (1:4)^2) - expected^2
  residual <- d$Score - expected
  for (p in seq_len(f$n_person)) {
    ii <- f$encoding$person0 + 1L == p
    expect_equal(f$person$ExpAvg[p], mean(expected[ii]), tolerance = 1e-10)
    expect_equal(f$person$InfitMnSq[p], sum(residual[ii]^2) / sum(variance[ii]), tolerance = 1e-10)
  }
  # On this centered balanced PCM, neutral facet locations give the fair score.
  neutral <- f; neutral$facets$Estimate[] <- 0
  fair <- predict(neutral)
  expected_fair <- vapply(seq_len(f$n_person), function(p) mean(fair[f$encoding$person0 + 1L == p]), numeric(1))
  expect_equal(f$person$FairM, expected_fair, tolerance = 1e-10)
  expect_equal(f$person$FairZ, expected_fair, tolerance = 1e-10)
  for (nm in unique(f$facets$Facet)) for (lev in unique(f$facets$Level[f$facets$Facet == nm])) {
    ii <- d[[nm]] == lev
    row <- f$facets[f$facets$Facet == nm & f$facets$Level == lev, ]
    expect_equal(row$InfitMnSq, sum(residual[ii]^2) / sum(variance[ii]), tolerance = 1e-10)
  }
})

test_that("prediction delta SE matches numerical perturbations for slopes and interactions", {
  d <- load_mfrmr_data("example_core")
  cases <- list(
    list(model = "PCM", step_facet = "Criterion"),
    list(model = "GPCM", step_facet = "Criterion", slope_facet = c("Criterion", "Rater"), spec = list(slope_reference = "Criterion")),
    list(model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion", spec = list(slope_type = "ability_only")),
    list(model = "RSM", facet_interactions = list(c("Rater", "Criterion")))
  )
  for (args in cases) {
    f <- do.call(fit_mfrm, c(list(data = d, person = "Person", facets = c("Rater", "Criterion"), score = "Score"), args))
    nd <- expand.grid(theta = c(-1, 0, 1), Rater = c("R01", "R04"), Criterion = c("Content", "Language"))
    # Perturb the physical facet, step and log-slope parameters along each free coordinate.
    perturb <- function(j, h) {
      z <- f; e <- f$encoding; l <- f$xi_layout
      for (k in seq_along(e$facets)) {
        pos <- match(j, l$facet_start[k] + seq_len(l$facet_len[k]))
        if (!is.na(pos)) {
          delta <- mfrmr:::.mfrm_facet_basis(e, k)[, pos] * h
          if (k <= length(e$base_facets)) {
            ii <- z$facets$Facet == e$facets[k]
            z$facets$Estimate[ii] <- z$facets$Estimate[ii] + delta
          } else z$interactions$Estimate[z$interactions$Interaction == e$facets[k]] <- z$interactions$Estimate[z$interactions$Interaction == e$facets[k]] + delta
        }
      }
      pos <- match(j, l$step_start + seq_len(l$step_len))
      if (!is.na(pos)) {
        B <- mfrmr:::.mfrm_step_basis(e, l)
        K <- e$n_cat - 1L
        for (g in seq_len(e$n_step_groups)) z$step_cumulative[g, ] <- z$step_cumulative[g, ] + c(0, cumsum(B[(g - 1L) * K + seq_len(K), pos] * h))
      }
      for (s in seq_len(e$n_slope_facets)) {
        pos <- match(j, l$slope_offset[s] + seq_len(e$slope_q[s]))
        if (!is.na(pos)) {
          ii <- z$slopes$Facet == e$slope_facets[s]
          z$slopes$Estimate[ii] <- z$slopes$Estimate[ii] * exp(e$slope_B[[s]][, pos] * h)
        }
      }
      predict(z, nd)
    }
    G <- vapply(seq_along(f$xi), function(j) (perturb(j, 1e-5) - perturb(j, -1e-5)) / 2e-5, numeric(nrow(nd)))
    expected <- sqrt(rowSums((G %*% f$xi_covariance_model) * G))
    expect_true(all(is.finite(expected)))
    expect_equal(predict(f, nd, se.fit = TRUE)$se.fit, expected, tolerance = 1e-8)
  }
})

test_that("prediction SE includes estimated steps when every facet is anchored", {
  d <- load_mfrmr_data("example_core")
  anchors <- do.call(rbind, lapply(c("Rater", "Criterion"), function(nm) data.frame(Facet = nm, Level = unique(d[[nm]]), Anchor = 0)))
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", anchors = anchors)
  nd <- data.frame(theta = c(-1, 0, 1), Rater = d$Rater[1], Criterion = d$Criterion[1])
  grad <- vapply(seq_len(nrow(f$steps)), function(j) {
    plus <- minus <- f
    cols <- (j + 1L):ncol(f$step_cumulative)
    plus$step_cumulative[1, cols] <- plus$step_cumulative[1, cols] + 1e-5
    minus$step_cumulative[1, cols] <- minus$step_cumulative[1, cols] - 1e-5
    (predict(plus, nd) - predict(minus, nd)) / 2e-5
  }, numeric(nrow(nd)))
  V <- vcov(f, component = "steps", se_type = "model")
  expected <- sqrt(rowSums((grad %*% V) * grad))
  expect_equal(predict(f, nd, se.fit = TRUE)$se.fit, expected, tolerance = 1e-8)
  expect_true(all(expected > 0))
})

test_that("common reference plans affect scores and retain normalized weights", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score")
  panel <- data.frame(Rater = c("R01", "R04"), Criterion = c("Content", "Language"), Weight = c(1, 3))
  a <- fair_average_table(f, reference_plan = panel)
  b <- fair_average_table(f, reference_plan = panel[1, ])
  expect_equal(a$reference_plan$NormalizedWeight, c(.25, .75))
  expected <- vapply(f$person$Estimate, function(theta) sum(predict(f, transform(panel, theta = theta)) * c(.25, .75)), numeric(1))
  expect_equal(a$reference_scores$Estimate, expected)
  expect_gt(max(abs(a$reference_scores$Estimate - b$reference_scores$Estimate)), .01)
  expect_true(all(a$reference_scores$Units == "original_score"))
  expect_error(fair_average_table(f, reference_plan = transform(panel, Rater = "unknown")), "Unknown")
  expect_error(fair_average_table(f, reference_plan = transform(panel, Weight = 0)), "positive")
})

test_that("agreement rejects repeated ratings instead of taking the last row", {
  d <- data.frame(Person = rep(1:4, each = 2), Rater = rep(c("A", "B"), 4), Score = rep(1:4, each = 2))
  x <- rbind(d, transform(d[1, ], Score = 4))
  expect_error(interrater_agreement(x, metric = "exact"), "Repeated ratings")
  expect_error(interrater_agreement(x[nrow(x):1, ], metric = "exact"), "Repeated ratings")
  x <- rbind(transform(d, Occasion = 1), transform(d, Occasion = 2))
  expect_no_error(interrater_agreement(x, unit = c("Person", "Occasion"), metric = "exact"))
})


test_that("G-Study weights never override the requested estimator", {
  d <- load_mfrmr_data("example_core")
  for (method in c("anova", "minque0")) {
    expect_error(mfrm_g_study(d, "Score", facets = c("Rater", "Criterion"),
      method = method, weights = 1), "requires method = 'urgenova'", fixed = TRUE)
  }
  expect_error(mfrm_g_study(d, "Score", facets = c("Rater", "Criterion"), weights = 1),
    "requires method = 'urgenova'", fixed = TRUE)
  a <- mfrm_g_study(d, "Score", facets = c("Rater", "Criterion"), method = "urgenova", weights = 1)
  b <- mfrm_g_study(d, "Score", facets = c("Rater", "Criterion"), method = "urgenova")
  expect_equal(a$variance_components, b$variance_components)
})

test_that("small Satterthwaite df is visible in intervals, saved tidy tables and reports", {
  d <- load_mfrmr_data("example_core")
  wf <- step_mfrm_gtheory(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score"))
  g <- wf$gtheory
  ci <- confint(g, component = "variance_components")
  v <- g$variance_components
  unstable <- which(is.finite(v$df) & v$df < 1)
  expect_length(unstable, 1L)
  expect_identical(v$Component[unstable], "Person:Criterion")
  expect_true(all(is.na(ci[unstable, ])))
  expect_identical(unname(attr(ci, "status")[unstable]), "unstable_satterthwaite_df")
  stable <- which(is.finite(v$df) & v$df >= 1)
  expected <- v$df[stable] * v$Estimate[stable] / qchisq(.025, v$df[stable])
  expect_equal(as.numeric(ci[stable, 2L]), expected, ignore_attr = TRUE)
  expect_identical(unname(attr(confint(g, parm = "Person:Criterion"), "status")),
    "unstable_satterthwaite_df")
  expect_error(confint(g, level = 1), "strictly between")
  expect_error(confint(g, method = "mls"), "must be empty")
  saved <- unserialize(serialize(wf, NULL))
  tab <- tidy(saved, "g_study")
  expect_true(is.na(tab$conf.high[unstable]))
  expect_identical(tab$interval_status[unstable], "unstable_satterthwaite_df")
  expect_identical(tidy(g), tidy(unserialize(serialize(g, NULL))))
  for (format in c("markdown", "html")) {
    report <- render_mfrm_report(saved, format = format)
    expect_match(report, "unstable_satterthwaite_df", fixed = TRUE)
  }
})

test_that("shared tidy generics dispatch workflows and workflow fits have an explicit contract", {
  d <- load_mfrmr_data("example_core")
  wf <- step_mfrm_measure(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score"))
  for (fun in c("tidy", "glance", "augment"))
    expect_identical(getExportedValue("mfrmr", fun), getExportedValue("generics", fun))
  expect_equal(generics::tidy(wf, "facets"), tidy(wf, "facets"))
  expect_equal(generics::glance(wf), glance(wf))
  expect_equal(generics::augment(wf), augment(wf))
  expect_equal(generics::tidy(wf$fit), tidy(wf$fit))
  expect_s3_class(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score", fit = wf$fit), "tidymfrm")
  canonical <- make_toy_fit()
  canonical$schema_version <- NULL
  expect_error(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score", fit = canonical),
    "Refit older saved objects")
})

test_that("calibrations and their diagnostics share person-fit tables, summaries, flags and plots", {
  d <- load_mfrmr_data("example_core")
  results <- lapply(c(TRUE, FALSE), function(use_diagnostics) {
    f <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "JML"))
    z <- compute_person_fit_indices(if (use_diagnostics) diagnose_mfrm(f) else f)
    expect_s3_class(z, "data.frame")
    expect_identical(z$table, as.data.frame(z))
    expect_identical(summary(z)$overview$Persons, nrow(z))
    expect_identical(z$Flagged, z$ReportFlag)
    expect_equal(z$ReportFlag, abs(z$ReportValue) > qnorm(.975))
    p <- plot(f, type = "person_fit")
    expect_equal(p$data$ReportValue, z$ReportValue)
    expect_true(all(z$lz_star_status == "computed_jml_conditional_calibration"))
    list(fit = f, result = z)
  })
  expect_identical(names(results[[1]]$result), names(results[[2]]$result))
  f <- results[[1]]$fit
  z <- results[[1]]$result
  pr <- predict(f, type = "probs")
  # Independent category-by-category moments of log probability and ability score.
  ii <- which(f$encoding$person0 == 0L)
  moments <- t(vapply(ii, function(i) {
    p <- unname(pr[i, ]); k <- seq_along(p) - 1L
    lp <- log(p); m <- sum(p * lp)
    u <- f$alpha_obs[i] * (k - sum(p * k))
    c(E = m, V = sum(p * (lp - m)^2), I = sum(p * u^2), C = sum(p * (lp - m) * u),
      U = u[f$encoding$score0[i] + 1L], L = lp[f$encoding$score0[i] + 1L])
  }, numeric(6)))
  m <- colSums(moments)
  expected <- (m['L'] - m['E'] - m['C'] / m['I'] * m['U']) / sqrt(m['V'] - m['C']^2 / m['I'])
  expect_equal(z$lz_star[match(f$encoding$person_levels[1], z$Person)], unname(expected))
  obs <- .mfrm_person_fit_moments(f)
  obs$PrObserved[ii[1]] <- NA_real_
  incomplete <- compute_person_fit_indices(list(obs = obs))
  row <- incomplete[incomplete$Person == obs$Person[ii[1]], ]
  expect_equal(row$UnavailableN, 1L)
  expect_equal(row$N, length(ii) - 1L)
  expect_true(is.na(row$lz))
  dg <- diagnose_mfrm(f)
  mismatch <- f; mismatch$person$Estimate[1] <- mismatch$person$Estimate[1] + .1
  expect_error(compute_person_fit_indices(dg, mismatch), "source calibration")
  for (scoring in c("EAP", "WLE", "JML_xtreme", "JML_jeffreys")) {
    changed <- f; changed$score_method <- scoring
    cz <- compute_person_fit_indices(changed)
    expect_true(all(is.na(cz$lz_star)))
    expect_true(all(grepl("^not_applicable", cz$lz_star_status)))
  }
  changed <- f; changed$converged <- FALSE
  expect_true(all(compute_person_fit_indices(changed)$lz_star_status == "fit_not_converged"))
})

test_that("GMFRM alias preserves model encoding and probability score derivatives", {
  d <- load_mfrmr_data("example_core")
  for (form in c("full_predictor", "ability_only")) {
    args <- list(data = d, person = "Person", facets = c("Rater", "Criterion"), score = "Score",
      step_facet = "Criterion", slope_facet = c("Rater", "Criterion"),
      spec = list(slope_type = form, gpcm_mml_identification = "fixed_standard_normal"))
    f <- do.call(fit_mfrm, c(args, list(model = "GMFRM")))
    g <- do.call(fit_mfrm, c(args, list(model = "GPCM")))
    expect_identical(f$model, "GPCM")
    expect_equal(f$xi, g$xi)
    obs <- .mfrm_person_fit_moments(f)
    plus <- minus <- f
    h <- 1e-5
    plus$person$Estimate <- plus$person$Estimate + h
    minus$person$Estimate <- minus$person$Estimate - h
    gradient <- (log(predict(plus, type = "probs")) - log(predict(minus, type = "probs"))) / (2*h)
    pr <- predict(f, type = "probs")
    expect_equal(obs$ScoreInformation, rowSums(pr * gradient^2), tolerance = 1e-7)
    expect_equal(obs$ObservedScoreDerivative,
      gradient[cbind(seq_len(nrow(d)), f$encoding$score0 + 1L)], tolerance = 1e-7)
  }
})

test_that("explicit G-Study weights preserve the composite scale and names", {
  d <- expand.grid(Person = 1:8, Rater = 1:4)
  d$A <- d$Person + sin(seq_len(nrow(d)))
  d$B <- .7 * d$Person + cos(seq_len(nrow(d))) + .2*d$Rater
  g <- function(weights) mfrm_g_study(d, c("A", "B"), "Person", facets = "Rater",
                                     method = "urgenova", weights = weights)
  average <- g(NULL); total <- g(c(1,1))
  direct <- transform(d, Score = A + B)
  sumfit <- mfrm_g_study(direct, "Score", "Person", facets = "Rater", method = "urgenova")
  expect_equal(total$variance_components$Estimate, sumfit$variance_components$Estimate)
  expect_equal(total$variance_components$Estimate, 4 * average$variance_components$Estimate)
  expect_equal(total$variance_components$SE, 4 * average$variance_components$SE)
  expect_equal(total$g_study_summary[c("RelSEM", "AbsSEM")], 2 * average$g_study_summary[c("RelSEM", "AbsSEM")])
  expect_equal(total$g_study_summary[c("G", "Phi")], average$g_study_summary[c("G", "Phi")])
  expect_equal(g(c(B = .2, A = .8))$g_study_summary, g(c(.8,.2))$g_study_summary)
  expect_error(g(c(X = .2, A = .8)), "match.*exactly")
  expect_error(g(c(".5", ".5")), "numeric")
})

test_that("readiness separates population-only links, effective anchors and observation weights", {
  d <- load_mfrmr_data("example_core")
  panel <- d[(as.integer(factor(d$Person)) <= 24) == (as.integer(factor(d$Rater)) <= 2), ]
  f <- suppressWarnings(fit_mfrm(panel, "Person", c("Rater", "Criterion"), "Score"))
  expect_false(f$readiness$ready)
  expect_true("population_assumption_linked" %in% f$readiness$reason_codes)
  expect_error(confint(f), "population_assumption_linked")
  # Direct design-rank reference with explicit free Person columns.
  enc <- f$encoding
  B <- do.call(cbind, lapply(seq_along(enc$facets), function(k)
    enc$facet_signs[k] * .mfrm_facet_basis(enc,k)[enc$facet_mat0[,k]+1L,,drop=FALSE]))
  P0 <- model.matrix(~ factor(enc$person0) - 1)
  P <- P0 %*% rbind(diag(ncol(P0) - 1L), -1)
  expect_equal(f$conditional_location_rank$rank, qr(cbind(P,B))$rank - ncol(P))
  anc <- data.frame(Facet = "Rater", Level = sort(unique(panel$Rater)), Anchor = 0)
  anchored <- suppressWarnings(fit_mfrm(panel, "Person", c("Rater", "Criterion"), "Score", anchors = anc))
  expect_true(anchored$readiness$ready)
  d$Weight <- rep(c(.5, 1, 2), length.out = nrow(d))
  weighted <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", weight = "Weight")
  expect_true(all(is.finite(weighted$xi)))
  expect_false(weighted$readiness$ready)
  expect_true("nonunit_observation_weights_inference_unvalidated" %in% weighted$readiness$reason_codes)
  old <- weighted; old$conditional_location_rank <- NULL
  expect_identical(.mfrm_readiness(old), weighted$readiness)
})

test_that("report table text round-trips delimiters without introducing markup", {
  values <- c("A|B", "|ZSTD| >= 2", "a\\|b", "a\\\\|b", "<script>x</script>", "**plain**", "a\nb", "")
  d <- data.frame(values, check.names = FALSE)
  names(d) <- "Criterion|Name"
  md <- .tidymfrm_fmt_df_md(d)
  html <- .tidymfrm_md_to_html(md, "Title <x>")
  expect_match(md, "A\\|B", fixed = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr("<td>", html, fixed = TRUE))), length(values))
  expect_match(html, "<th>Criterion|Name</th>", fixed = TRUE)
  expect_match(html, "<td>a\\|b</td>", fixed = TRUE)
  expect_match(html, "<td>a\\\\|b</td>", fixed = TRUE)
  expect_match(html, "&lt;script&gt;x&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "<td>**plain**</td>", fixed = TRUE)
  expect_match(html, "<td>a b</td>", fixed = TRUE)
  expect_false(grepl("_No rows", .tidymfrm_md_to_html(.tidymfrm_fmt_df_md(d[FALSE,,drop=FALSE]), "Empty"), fixed = TRUE))
})

test_that("rater-card central tendency is within criterion and halo has no universal cutoff", {
  wf <- make_core_workflow()
  rater <- as.character(wf$diagnostics$raters$Level[1L])
  card <- .mfrm_rater_card(wf, rater, "Criterion", 5L)
  table <- card$tables[["Central tendency screen"]]
  obs <- wf$diagnostics$observations
  obs <- obs[as.character(obs$Rater) == rater, ]
  expect_equal(nrow(table), length(unique(obs$Criterion)))
  for (i in seq_len(nrow(table))) {
    z <- obs[obs$Criterion == table$Criterion[i], ]
    expect_equal(table$ResidualExpectedCorrelation[i], cor(z$Residual, z$Expected))
  }
  changed <- wf
  delta <- 10 * as.integer(factor(changed$diagnostics$observations$Criterion))
  changed$diagnostics$observations$Expected <- changed$diagnostics$observations$Expected + delta
  changed$diagnostics$observations$Residual <- changed$diagnostics$observations$Residual - delta
  other <- .mfrm_rater_card(changed, rater, "Criterion", 5L)
  expect_equal(other$tables[["Central tendency screen"]], table)
  expect_match(paste(card$notes, collapse = " "), "Neither a zero cutoff nor -1/(K-1) is calibrated", fixed = TRUE)
})
