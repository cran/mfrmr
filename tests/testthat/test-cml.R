cml_example <- function() {
  patterns <- as.matrix(expand.grid(rep(list(0:2), 4L)))
  patterns <- patterns[rep(seq_len(nrow(patterns)), 1L + 2L * (patterns[, 1] == 2) +
    (patterns[, 2] == 0)), , drop = FALSE]
  data.frame(Person = rep(sprintf("P%03d", seq_len(nrow(patterns))), each = 4L),
    Rater = rep(c("A", "A", "B", "B"), nrow(patterns)),
    Item = rep(c("X", "Y", "X", "Y"), nrow(patterns)), Score = as.vector(t(patterns)))
}

cml_problem_for <- function(data, model = "RSM", owner = NULL, positive = NULL, interactions = NULL) {
  prep <- prepare_mfrm_data(data, "Person", c("Rater", "Item"), "Score",
    rating_min = 0, rating_max = 2, keep_original = TRUE)
  specs <- build_facet_interaction_specs(prep, interactions,
    min_obs_per_interaction = 1, interaction_policy = "silent")
  built <- build_estimation_config(prep, model, "CML", owner, NULL,
    interaction_specs = specs, weight_col = NULL,
    facet_signs = build_facet_signs(prep$facet_names, positive)$signs,
    positive_facets = positive, noncenter_facet = "Person", dummy_facets = NULL,
    anchor_df = NULL, group_anchor_df = NULL)
  idx <- build_indices(prep, owner, interaction_specs = specs)
  list(problem = mfrmr_cml_problem(prep, idx, built$config, built$sizes),
    prep = prep, idx = idx, config = built$config, sizes = built$sizes)
}

cml_enumerate <- function(features, par) {
  patterns <- as.matrix(expand.grid(lapply(features, function(x) seq_len(nrow(x)))))
  stats <- t(apply(patterns, 1L, function(row) {
    Reduce(`+`, lapply(seq_along(row), function(j) features[[j]][row[j], ]))
  }))
  total <- rowSums(patterns - 1L)
  lapply(0:max(total), function(s) {
    T <- stats[total == s, , drop = FALSE]
    eta <- drop(T %*% par); w <- exp(eta - max(eta)); w <- w / sum(w)
    mu <- drop(crossprod(w, T))
    list(logz = max(eta) + log(sum(exp(eta - max(eta)))), means = mu,
      covariance = crossprod(sweep(T, 2L, mu) * sqrt(w)))
  })
}

test_that("CML retains a stationary solution when conditional information is ill-conditioned", {
  # Three binary ratings at total one give a three-cell multinomial likelihood.
  # Aggregate integer multiplicities exercise conditioning without allocating N rows.
  N <- 1e9
  features <- list(rbind(c(0,0), c(0,0)), rbind(c(0,0), c(1,1)), rbind(c(0,0), c(1,-1)))
  problem <- list(parameters=2L, observed=c(N/2, N/2-2),
    groups=list(list(features=features, counts=c(0, N, 0, 0))))
  logits <- log(c(1-2/N, 2/N))
  exact <- c(sum(logits)/2, (logits[1]-logits[2])/2)
  expect_warning(fit <- mfrmr_cml_optimize(problem, exact,
    list(optimizer="BFGS", maxit=300L, reltol=1e-10)), "point estimates retained.*ill-conditioned")
  expect_equal(fit$opt$par, exact, tolerance=1e-7)
  expect_equal(fit$opt$convergence, 0L)
  expect_lt(fit$checks$gradient_max, 1e-6)
  expect_identical(fit$checks$information_status, "unavailable")
  expect_true(all(is.na(fit$covariance)))
  expect_true(all(is.finite(fit$information)))
  expect_lt(min(fit$checks$information_eigenvalues), 1e-8*max(fit$checks$information_eigenvalues))
})

test_that("information-only failures preserve CML points, reasons and public consumers", {
  skip_if_not_installed("lpSolve")
  # Equal frequencies of every pattern make zero the exact conditional root.
  patterns <- as.matrix(expand.grid(rep(list(0:2), 4L)))
  d <- data.frame(Person=rep(seq_len(nrow(patterns)), each=4L),
    Rater=rep(c("A","A","B","B"), nrow(patterns)),
    Item=rep(c("X","Y","X","Y"), nrow(patterns)), Score=as.vector(t(patterns)))
  evaluate <- mfrmr_cml_evaluate
  mode <- "evaluation"; information_calls <- 0L
  local_mocked_bindings(mfrmr_cml_evaluate=function(problem, par, information=FALSE) {
    if (information) information_calls <<- information_calls+1L
    # Preserve information needed for Newton refinement; inject only after
    # the actual score has passed the unchanged stationarity criterion.
    fail <- information && information_calls > 1L &&
      max(abs(evaluate(problem, par)$gradient)) <= 1e-6
    if (fail && mode == "evaluation")
      stop("Controlled conditional-information evaluation failure.")
    out <- evaluate(problem, par, information)
    if (fail && mode != "evaluation")
      out$information[,] <- if (mode == "singular") 0 else NA_real_
    out
  })
  for (mode in c("evaluation", "singular", "nonfinite")) {
    information_calls <- 0L
    expect_warning(fit <- fit_mfrm(d,"Person",c("Rater","Item"),"Score",
      method="CML", category_policy="preserve", facet_interactions=list(c("Rater","Item")),
      reltol=1e-12), "point estimates retained; uncertainty unavailable")
    expect_equal(unname(coef(fit)), rep(0, length(coef(fit))), tolerance=1e-7)
    expect_true(is.finite(as.numeric(logLik(fit))))
    expect_identical(fit$checks$information_status, "unavailable")
    expect_true(nzchar(fit$checks$information_reason))
    s <- summary(fit)
    for (tab in s[c("locations","steps","interactions")]) {
      expect_true(all(is.finite(tab$Estimate)))
      expect_true(all(is.na(tab$SE) & is.na(tab$Lower) & is.na(tab$Upper)))
      expect_false(any(tab$CIEligible))
      expect_true(all(grepl(fit$checks$information_reason, tab$InferenceReview, fixed=TRUE)))
    }
    expect_output(print(fit), "information_unavailable")
    ci <- mfrm_facet_intervals(fit,"Rater")
    expect_true(all(is.finite(ci$table$Estimate)))
    expect_true(all(ci$table$Status == "information_unavailable"))
    expect_identical(plot_data(plot(ci,draw=FALSE),"table"), ci$table)
    expect_true(all(is.na(apa_table(ci)$table$SE)))
    V <- vcov(fit); limits <- confint(fit)
    expect_true(all(is.na(V)) && all(is.na(limits)))
    expect_identical(attr(V,"unavailable_reason"), fit$checks$information_reason)
    expect_true(all(grepl(fit$checks$information_reason,
      attr(limits,"diagnostics")$InferenceReview, fixed=TRUE)))
    fixed <- mfrm_facet_intervals(fit,"Rater",
      matrix(c(1,1),1,dimnames=list("sum",c("A","B"))))$table
    expect_identical(fixed$Status, "fixed")
    expect_equal(fixed$SE, 0)
    expect_true(is.na(fixed$Lower))
    expect_identical(summary(unserialize(serialize(fit,NULL))), s)
  }
})

test_that("CML ESF, score and information agree with enumeration and finite differences", {
  d <- cml_example()
  for (model in c("RSM", "PCM")) for (positive in list(NULL, "Rater")) {
    b <- cml_problem_for(d, model, if (model == "PCM") "Item", positive)
    prob <- b$problem; par <- seq(-.4, .3, length.out = prob$parameters)
    x <- prob$groups[[1L]]$features
    actual <- mfrmr_cml_esf(x, par, TRUE); exact <- cml_enumerate(x, par)
    expect_equal(actual$logz, vapply(exact, `[[`, numeric(1), "logz"), tolerance = 1e-12)
    expect_equal(actual$means, do.call(rbind, lapply(exact, `[[`, "means")), tolerance = 1e-12)
    for (s in seq_along(exact)) expect_equal(actual$covariance[, , s], exact[[s]]$covariance, tolerance = 1e-12)
    at <- mfrmr_cml_evaluate(prob, par, TRUE)
    step <- diag(prob$parameters) * 1e-5
    gd <- vapply(seq_along(par), function(j)
      (mfrmr_cml_evaluate(prob, par + step[, j])$value -
         mfrmr_cml_evaluate(prob, par - step[, j])$value) / 2e-5, numeric(1))
    hd <- vapply(seq_along(par), function(j)
      (mfrmr_cml_evaluate(prob, par + step[, j])$gradient -
         mfrmr_cml_evaluate(prob, par - step[, j])$gradient) / 2e-5, numeric(length(par)))
    expect_equal(at$gradient, gd, tolerance = 1e-7)
    expect_equal(at$information, hd, tolerance = 1e-7)
    # Independent native category logits, rather than a second use of the map.
    expanded <- expand_params(par, b$sizes, b$config)
    for (j in seq_along(x)) {
      # Cell order is the native Rater:Item integer order, sorted by roster.
      at_row <- which(b$idx$facets$Rater == c(1, 1, 2, 2)[j] &
        b$idx$facets$Item == c(1, 2, 1, 2)[j])[1L]
      eta <- sum(vapply(b$config$facet_names, function(f)
        b$config$facet_signs[[f]] * expanded$facets[[f]][b$idx$facets[[f]][at_row]], numeric(1)))
      steps <- if (model == "RSM") expanded$steps else expanded$steps_mat[b$idx$step_idx[at_row], ]
      expect_equal(drop(x[[j]] %*% par), c(0, cumsum(eta - steps)), tolerance = 1e-12)
    }
  }
  # Log recursion stays finite where exponentiating category weights would overflow.
  stable <- mfrmr_cml_esf(list(matrix(c(0, 1, 2), ncol = 1), matrix(c(0, -1, 2), ncol = 1)), 1000, TRUE)
  expect_true(all(is.finite(unlist(stable))))
})

test_that("one-coordinate CML matches the paired binary analytic likelihood", {
  skip_if_not_installed("lpSolve")
  paired <- function(a, b) {
    patterns <- rbind(c(1, 0), c(0, 1), c(0, 0), c(1, 1))
    x <- patterns[rep(seq_len(4L), c(a, b, 2L, 3L)), , drop = FALSE]
    data.frame(Person = rep(seq_len(nrow(x)), each = 2L),
      Rater = rep(c("A", "B"), nrow(x)), Score = as.vector(t(x)))
  }
  fit <- function(d, model = "RSM") fit_mfrm(d, "Person", "Rater", "Score",
    method = "CML", model = model, step_facet = if (model == "PCM") "Rater",
    category_policy = "preserve")
  # For severity (d, -d), P(10 | total = 1) = plogis(-2*d).
  # Two discordant responses still have a finite MLE; information, rather
  # than total Person count (which includes extremes), determines the SE.
  for (model in c("RSM", "PCM")) for (counts in list(c(3L, 7L), c(1L, 1L), c(1L, 999L))) {
    a <- counts[1L]; b <- counts[2L]; n <- a + b
    f <- fit(paired(a, b), model)
    estimate <- log(b / a) / 2
    information <- 4 * a * b / n
    expect_length(f$opt$par, 1L)
    expect_equal(unname(f$opt$par), estimate, tolerance = 1e-6)
    expect_equal(as.numeric(f$information), information, tolerance = 1e-6)
    expect_equal(as.numeric(vcov(f)), 1 / information, tolerance = 1e-6)
    expect_equal(as.numeric(logLik(f)), a * log(a / n) + b * log(b / n), tolerance = 1e-10)
    expect_equal(sum(f$person_totals$Extreme), 5L)
    expect_identical(f$checks$boundary$state, "no_recession_detected")
    ci <- mfrm_facet_intervals(f, "Rater",
      matrix(c(1, -1), 1L, dimnames = list("A-B", c("A", "B"))))$table
    expect_equal(ci$Estimate, 2 * estimate, tolerance = 1e-6)
    expect_equal(ci$SE, sqrt(1 / a + 1 / b), tolerance = 1e-6)
    expect_equal(ci$Lower, 2 * estimate - qnorm(.975) * sqrt(1 / a + 1 / b), tolerance = 1e-6)
  }
  expect_error(fit(paired(0L, 10L)), "recession_detected")
  expect_error(fit(paired(10L, 0L)), "recession_detected")
  expect_error(fit(paired(0L, 0L)), "conditioned CML model is unidentified")
})

test_that("CML counts independent repeated-cell patterns with their multiplicity", {
  skip_if_not_installed("lpSolve")
  # Two distinct A events and one B event per Person, all at total one.
  # P(A total = 1 | total = 1) = 2*exp(-d)/(2*exp(-d) + exp(d));
  # each of the two A patterns has half that probability.
  a <- 6L; b <- 7L
  patterns <- diag(3L)[rep(seq_len(3L), c(a / 2L, a / 2L, b)), , drop = FALSE]
  d <- data.frame(Person = rep(seq_len(nrow(patterns)), each = 3L),
    Rater = rep(c("A", "A", "B"), nrow(patterns)), Score = as.vector(t(patterns)))
  fit <- function(data) fit_mfrm(data, "Person", "Rater", "Score", method = "CML",
    repeated_ratings = "independent", category_policy = "preserve")
  f <- fit(d)
  expect_equal(unname(f$opt$par), log(2 * b / a) / 2, tolerance = 1e-6)
  expect_equal(as.numeric(f$information), 4 * a * b / (a + b), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f)),
    a * log(a / (2 * (a + b))) + b * log(b / (a + b)), tolerance = 1e-10)
  # Duplicating independent Persons multiplies information, not support;
  # it differs from adding observations to each Person's conditioning roster.
  additional <- d; additional$Person <- additional$Person + nrow(patterns)
  replicated <- fit(rbind(d, additional))
  expect_equal(replicated$opt$par, f$opt$par, tolerance = 1e-6)
  expect_equal(as.numeric(replicated$information), 2 * as.numeric(f$information), tolerance = 1e-6)
  expect_equal(as.numeric(vcov(replicated)), as.numeric(vcov(f)) / 2, tolerance = 1e-6)
  expect_equal(as.numeric(logLik(replicated)), 2 * as.numeric(logLik(f)), tolerance = 1e-10)
})

test_that("native CML fits expose conditional covariance and reusable facet intervals", {
  skip_if_not_installed("lpSolve")
  d <- cml_example()
  for (model in c("RSM", "PCM")) {
    fit <- fit_mfrm(d, "Person", c("Rater", "Item"), "Score", method = "CML",
      model = model, step_facet = if (model == "PCM") "Item", category_policy = "preserve")
    expect_s3_class(fit, "mfrm_cml_fit")
    expect_false(inherits(fit, "mfrm_fit"))
    expect_equal(fit$checks$boundary$state, "no_recession_detected")
    expect_lt(fit$checks$gradient_max, 1e-4)
    expect_equal(unname(fit$information %*% vcov(fit)), diag(length(fit$opt$par)), tolerance = 1e-9)
    expect_equal(as.numeric(logLik(fit)), -fit$opt$value)
    s <- summary(fit)
    expect_true(all(is.finite(s$locations$SE)))
    expect_true(all(is.finite(s$steps$SE)))
    C <- matrix(c(1, -1), 1, dimnames = list("A_minus_B", c("A", "B")))
    ci <- mfrm_facet_intervals(fit, "Rater", C)
    expect_error(mfrm_results(fit), "does not support CML fits.*summary\\(fit\\)")
    expect_error(mfrm_results(fit, intervals = ci), "plot\\(ci\\).*saveRDS")
    expect_equal(ci$table$Estimate, unname(2 * fit$facets$others$Estimate[1L]))
    expect_equal(ci$covariance, ci$target_jacobian %*% vcov(fit) %*% t(ci$target_jacobian))
    expect_s3_class(plot(ci, draw = FALSE), "mfrm_plot_data")
    notes <- plot(ci, draw = FALSE)$data$subtitle
    expect_match(notes, "Conditional information")
    expect_false(grepl("Weak information", notes))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      graphic <- ggplot2::ggplot_build(as_ggplot(plot(ci, draw = FALSE)))
      expect_true(all(!is.na(graphic$data[[1L]]$colour)))
      expect_true(all(!is.na(graphic$data[[1L]]$linetype)))
    }
    path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
    saveRDS(list(fit = fit, ci = ci), path)
    saved <- readRDS(path)
    expect_identical(summary(saved$fit), s)
    expect_identical(mfrm_facet_intervals(saved$fit, "Rater", C)$table, ci$table)
    expect_error(mfrm_facet_intervals(fit, "Rater", C, scale = "standardized"), "does not estimate a population SD")
    expect_error(mfrm_facet_intervals(fit, "Rater", method = "sandwich"), "require method")
    expect_error(mfrm_facet_intervals(fit, "Rater", level = 1), "0 < level")
    for (level in list(matrix(.95), structure(.95, class = "confidence")))
      expect_error(mfrm_facet_intervals(fit, "Rater", level = level), "0 < level")
    expect_identical(mfrm_facet_intervals(fit, "Rater", C, level = c(confidence = .95)), ci)
    high_level <- 1 - .Machine$double.eps / 2
    high <- mfrm_facet_intervals(fit, "Rater", C, level = high_level)$table
    critical <- (high$Upper - high$Estimate) / high$SE
    expect_equal(pnorm(critical, lower.tail = FALSE, log.p = TRUE), log1p(-high_level) - log(2))
    expect_true(is.finite(high$Lower) && is.finite(high$Upper))
    scaled <- mfrm_facet_intervals(fit, "Rater", C * 1e-100)$table
    expect_equal(scaled$SE / 1e-100, ci$table$SE, tolerance = 1e-12)
    expect_equal(scaled$Lower / 1e-100, ci$table$Lower, tolerance = 1e-12)
    expect_error(mfrm_facet_intervals(fit, "Rater", C * 1e-200), "positive representable variance")
    fixed <- mfrm_facet_intervals(fit, "Rater", C * 0 + 1)$table
    expect_equal(fixed$Status, "fixed")
    expect_equal(fixed$SE, 0)
    expect_true(is.na(fixed$Lower))
    expect_error(mfrm_facet_intervals(fit, "Rater", C * 1e308), "must be finite")
  }
})

test_that("CML refuses unidentified, boundary, unsupported and unconverged models", {
  skip_if_not_installed("lpSolve")
  d <- cml_example()
  fit <- function(data = d, ...) fit_mfrm(data, "Person", c("Rater", "Item"), "Score",
    method = "CML", category_policy = "preserve", ...)
  expect_error(fit(model = "GPCM"), "RSM/PCM")
  expect_error(fit(population_formula = ~1), "does not accept.*population_formula")
  expect_error(fit(anchors = NULL), "does not accept.*anchors")
  expect_error(fit(fit_scope = "comparisons"), "Supply the comparison_review")
  expect_error(fit(model = "PCM"), "requires .step_facet.")
  expect_error(fit(maxit = 1), "stationary")
  weighted <- d; weighted$W <- 2
  expect_error(fit(weighted, weight = "W"), "unit weights")
  repeated <- rbind(d, d)
  expect_error(suppressWarnings(fit(repeated)), "resolution of repeated rating cells")
  independent <- fit(repeated, repeated_ratings = "independent")
  expect_s3_class(independent, "mfrm_cml_fit")
  expect_equal(nrow(independent$prep$data), nrow(repeated))
  nested <- d; nested$Rater <- ifelse(as.integer(factor(nested$Person)) %% 2 == 0, "A", "B")
  expect_error(suppressWarnings(fit(nested, repeated_ratings = "independent")), "conditioned CML model is unidentified")
  separated <- d; separated$Score[separated$Rater == "A"] <- 2
  expect_error(fit(separated), "recession_detected")
  extreme <- d; extreme$Score <- 2 * (as.integer(factor(d$Person)) %% 2)
  expect_error(fit(extreme, rating_min = 0, rating_max = 2), "conditioned CML model is unidentified")
})

test_that("CML preserves missing rosters, score origins, interactions and extreme accounting", {
  skip_if_not_installed("lpSolve")
  d <- cml_example()
  fit <- function(data, ...) fit_mfrm(data, "Person", c("Rater", "Item"), "Score",
    method = "CML", category_policy = "preserve", ...)
  missing <- d
  missing$Score[seq(4L, nrow(d), by = 28L)] <- NA
  a <- fit(missing, model = "PCM", step_facet = "Rater")
  b <- fit(missing[nrow(missing):1L, ], model = "PCM", step_facet = "Rater")
  expect_equal(a$opt$par, b$opt$par, tolerance = 1e-8)
  shifted <- missing; shifted$Score <- shifted$Score + 1
  c <- fit(shifted, model = "PCM", step_facet = "Rater", rating_min = 1, rating_max = 3)
  expect_equal(a$opt$par, c$opt$par, tolerance = 1e-8)
  expect_equal(a$covariance, c$covariance, tolerance = 1e-8)
  extreme <- with(a$person_totals, Person[Extreme])
  removed <- fit(missing[!missing$Person %in% extreme, ], model = "PCM", step_facet = "Rater")
  expect_equal(a$opt$par, removed$opt$par, tolerance = 1e-8)
  expect_equal(a$covariance, removed$covariance, tolerance = 1e-8)
  interaction <- fit(d, facet_interactions = list(c("Rater", "Item")))
  expect_true(all(is.finite(summary(interaction)$interactions$SE)))
  expect_equal(sum(interaction$interactions$effects$Estimate), 0, tolerance = 1e-12)
  expect_equal(sort(a$person_totals$Person), sort(unique(d$Person)))
  if (requireNamespace("eRm", quietly = TRUE)) {
    full <- cml_problem_for(d, "PCM", "Rater")
    W <- do.call(rbind, lapply(full$problem$groups[[1L]]$features, function(x) x[-1L, , drop = FALSE]))
    reference <- eRm::LPCM(matrix(missing$Score, ncol = 4L, byrow = TRUE), W = W, sum0 = FALSE)
    expect_lt(abs(-a$opt$value - reference$loglik), 1e-7)
    expect_lt(max(abs(a$opt$par - reference$etapar)), 1e-5)
    expect_lt(max(abs(a$covariance - solve(reference$hessian))), 1e-6)
  }
})

test_that("conditional DAG recession checks agree with enumerated support", {
  skip_if_not_installed("lpSolve")
  d <- cml_example()
  separated <- d; separated$Score[separated$Rater == "A"] <- 2
  for (data in list(d, separated)) {
    prob <- cml_problem_for(data)$problem
    zero <- mfrmr_cml_evaluate(prob, numeric(prob$parameters), TRUE)
    actual <- mfrmr_cml_boundary(prob, zero$gradient)
    features <- prob$groups[[1L]]$features
    patterns <- as.matrix(expand.grid(rep(list(0:2), 4L)))
    stats <- t(apply(patterns, 1, function(y)
      Reduce(`+`, lapply(seq_along(y), function(j) features[[j]][y[j] + 1L, ]))))
    obs <- matrix(data$Score, ncol = 4L, byrow = TRUE)
    constraints <- do.call(rbind, lapply(seq_len(nrow(obs)), function(i) {
      z <- Reduce(`+`, lapply(seq_len(4L), function(j) features[[j]][obs[i, j] + 1L, ]))
      sweep(stats[rowSums(patterns) == sum(obs[i, ]), , drop = FALSE], 2, z, function(a, b) b - a)
    }))
    reference <- lpSolve::lp("max", c(-zero$gradient, zero$gradient),
      rbind(cbind(constraints, -constraints), 1),
      c(rep(">=", nrow(constraints)), "<="), c(rep(0, nrow(constraints)), 1))
    expect_equal(reference$status, 0)
    expect_equal(actual$objective, reference$objval, tolerance = 1e-8)
    expect_equal(actual$state, if (reference$objval > 1e-7) "recession_detected" else "no_recession_detected")
  }
})

test_that("matched eRm LPCM independently reproduces multifacet CML fits and information", {
  skip_if_not_installed("lpSolve"); skip_if_not_installed("eRm")
  d <- cml_example()
  for (model in c("RSM", "PCM")) {
    fit <- fit_mfrm(d, "Person", c("Rater", "Item"), "Score", method = "CML",
      model = model, step_facet = if (model == "PCM") "Item", category_policy = "preserve")
    b <- cml_problem_for(d, model, if (model == "PCM") "Item")
    W <- do.call(rbind, lapply(b$problem$groups[[1L]]$features, function(x) x[-1L, , drop = FALSE]))
    X <- matrix(d$Score, ncol = 4L, byrow = TRUE)
    reference <- eRm::LPCM(X, W = W, sum0 = FALSE)
    expect_equal(unname(reference$W), W)
    expect_equal(-fit$opt$value, reference$loglik, tolerance = 1e-8)
    expect_lt(max(abs(fit$opt$par - reference$etapar)), 1e-5)
    # eRm 1.0-10 uses nlm's finite-difference Hessian; the exact information
    # above is separately checked against enumeration and score derivatives.
    expect_lt(max(abs(vcov(fit) - solve(reference$hessian))), 1e-6)
  }
})
