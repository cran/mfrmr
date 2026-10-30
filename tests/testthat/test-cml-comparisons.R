cml_comparison_example <- function() {
  do.call(rbind, lapply(c("A", "B"), function(g) {
    raters <- if (g == "A") c("R1", "R2", "R3") else c("R1", "R2")
    patterns <- as.matrix(expand.grid(rep(list(0:2), length(raters))))
    copies <- 1L + if (g == "A") 2L * (patterns[, 1L] == 2) else (patterns[, 2L] == 2)
    patterns <- patterns[rep(seq_len(nrow(patterns)), copies), , drop = FALSE]
    data.frame(Person = rep(paste0(g, seq_len(nrow(patterns))), each = length(raters)),
      Rater = rep(raters, nrow(patterns)), Group = g, Score = as.vector(t(patterns)))
  }))
}

cml_comparison_args <- local({
  cache <- list()
  function(model = "RSM", owner = NULL, positive = NULL) {
  key <- paste(model, paste(owner, collapse = "|"), paste(positive, collapse = "|"), sep = ":")
  if (!is.null(cache[[key]])) return(cache[[key]])
  d <- cml_comparison_example()
  W <- rbind(panel = c(1, -.5, -.5), pair = c(1, -1, 0), twice = c(2, -2, 0))
  colnames(W) <- c("R1", "R2", "R3")
  review <- describe_mfrm_data(d, person = "Person", facets = c("Rater", "Group"),
    score = "Score", rating_min = 0, rating_max = 2, category_policy = "preserve",
    comparisons = list(facet = "Rater", group = "Group", groups = c("A", "B"),
      contrasts = W, model = model, method = "CML", step_facet = owner,
      positive_facets = positive))
  cache[[key]] <<- c(list(data = d, fit_scope = "comparisons", comparison_review = review,
    interaction_policy = "silent"), review$comparison_review$fit_arguments)
  cache[[key]]
  }
})

# Explicitly enumerate support at the observed totals; no ESF or DP recursion.
cml_comparison_patterns <- function(problem) lapply(problem$groups, function(g) {
  patterns <- as.matrix(expand.grid(lapply(g$features, function(x) 0:(nrow(x) - 1L))))
  statistics <- t(apply(patterns, 1L, function(y)
    Reduce(`+`, lapply(seq_along(y), function(j) g$features[[j]][y[j] + 1L, ]))))
  lapply(which(g$counts > 0), function(s)
    list(statistics = statistics[rowSums(patterns) == s - 1L, , drop = FALSE], count = g$counts[s]))
})

test_that("CML comparison points survive an unavailable joint covariance", {
  skip_if_not_installed("lpSolve")
  a <- cml_comparison_args()
  reference <- do.call(fit_mfrm, a)
  optimize <- mfrmr_cml_optimize
  local_mocked_bindings(mfrmr_cml_optimize=function(...) {
    out <- optimize(...)
    out$covariance[,] <- NA_real_
    out$checks$information_status <- "unavailable"
    out$checks$information_reason <- "Controlled comparison-information failure."
    out
  })
  fit <- do.call(fit_mfrm, a)
  expect_equal(coef(fit), coef(reference))
  expect_identical(fit$comparisons$EstimateStatus, c("unavailable", "conditional_point", "conditional_point"))
  expect_false(any(fit$comparisons$FormalInferenceEligible))
  expect_true(all(is.na(fit$comparisons$SE)))
  expect_true(all(is.na(fit$covariance)))
  expect_true(all(fit$comparisons$NumericalState == "review"))
  expect_true(all(grepl("Controlled comparison-information failure", fit$comparisons$EstimateReason[-1])))
  ci <- mfrm_facet_intervals(fit, level=.9)
  expect_equal(ci$table$Estimate, reference$comparisons$Estimate)
  expect_false(any(ci$table$CIEligible))
  expect_true(all(is.na(ci$table$SE) & is.na(ci$table$Lower) & is.na(ci$table$Upper)))
  expect_identical(plot_data(plot(ci,draw=FALSE),"table"), ci$table)
  expect_identical(apa_table(ci)$table$CIEligible, ci$table$CIEligible)
  expect_true(all(is.na(vcov(fit))))
  limits <- confint(fit)
  expect_true(all(is.na(limits)))
  expect_true(all(grepl("Controlled comparison-information failure",
    attr(limits,"diagnostics")$InferenceReview[-1])))
  expect_true(grepl("unidentified", attr(limits,"diagnostics")$InferenceReview[1]))
  expect_identical(mfrm_facet_intervals(unserialize(serialize(fit,NULL))), mfrm_facet_intervals(fit))
})

test_that("conditional support DP agrees with enumerated patterns, including restricted totals", {
  set.seed(618)
  for (totals in list(0L, 1L, 2L, 3L, 6L, c(1L, 4L))) {
    features <- lapply(1:3, function(i) rbind(0, matrix(sample(-2:2, 8, TRUE), 2, 4)))
    counts <- numeric(7); counts[totals + 1L] <- 1
    p <- list(parameters = 4L, groups = list(list(features = features, counts = counts)))
    blocks <- unlist(cml_comparison_patterns(p), recursive = FALSE)
    exact <- do.call(rbind, lapply(blocks, function(b) sweep(b$statistics, 2, b$statistics[1, ])))
    actual <- mfrmr_cml_support_design(p)
    expect_equal(qr(actual)$rank, qr(exact)$rank)
    expect_equal(qr(rbind(actual, exact))$rank, qr(exact)$rank)
  }
  expect_null(mfrmr_cml_support_design(list(parameters = 1500L)))
  # Only categories 0/1 can occur at total 1; category 2 statistics add no rank.
  features <- rep(list(rbind(c(0, 0), c(1, 0), c(2, 1))), 3)
  D <- mfrmr_cml_support_design(list(parameters = 2L,
    groups = list(list(features = features, counts = c(0, 3, rep(0, 5))))))
  expect_equal(D, matrix(0, 1, 2))
})

test_that("partial CML contrasts and joint covariance agree with an independent QR and enumeration fit", {
  skip_if_not_installed("lpSolve")
  for (setting in list(list("RSM", NULL), list("PCM", "Rater"), list("PCM", "Group"))) {
    directions <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) list(NULL, "Rater") else
      if (identical(setting[[2]], "Rater")) list("Rater") else list(NULL)
    for (positive in directions) {
      a <- cml_comparison_args(setting[[1]], setting[[2]], positive)
      f <- do.call(fit_mfrm, a); tab <- interaction_effect_table(f)
      expect_identical(tab$ConditionalState, c("target_unidentified", rep("target_identified_only", 2)))
      expect_identical(tab$FormalInferenceEligible, c(FALSE, TRUE, TRUE))
      expect_true(all(is.na(unlist(tab[1, c("Estimate", "SE", "Lower", "Upper")]))))
      expect_true(all(is.na(f$covariance[1, ])) && all(is.na(f$covariance[, 1])))
      expect_equal(tab$Estimate[3], 2 * tab$Estimate[2], tolerance = 1e-12)
      expect_equal(f$covariance[2:3, 2:3], outer(c(1, 2), c(1, 2)) * tab$SE[2]^2,
        ignore_attr = TRUE, tolerance = 1e-12)
      expect_identical(summary(f), tab)
      expect_identical(interaction_effect_table(unserialize(serialize(f, NULL))), tab)
      expect_output(print(f), "pointwise 95% Wald")
      cfg <- f$config; sizes <- build_param_sizes(cfg)
      idx <- build_indices(f$prep, cfg$step_facet, interaction_specs = cfg$interaction_specs)
      p <- mfrmr_cml_problem(f$prep, idx, cfg, sizes)
      blocks <- unlist(cml_comparison_patterns(p), recursive = FALSE)
      D <- do.call(rbind, lapply(blocks, function(b) sweep(b$statistics, 2, b$statistics[1, ])))
      qrD <- qr(D, tol = 1e-10)
      kept <- qrD$pivot[seq_len(qrD$rank)]
      ev <- function(z) {
        value <- -sum(p$observed[kept] * z); gradient <- -p$observed[kept]
        H <- matrix(0, length(z), length(z))
        for (b in blocks) {
          X <- b$statistics[, kept, drop = FALSE]
          eta <- drop(X %*% z); w <- exp(eta - max(eta)); w <- w / sum(w)
          mu <- drop(crossprod(w, X))
          value <- value + b$count * (max(eta) + log(sum(exp(eta - max(eta)))))
          gradient <- gradient + b$count * mu
          H <- H + b$count * crossprod(sweep(X, 2, mu) * sqrt(w))
        }
        list(value = value, gradient = gradient, information = H)
      }
      opt <- stats::optim(numeric(length(kept)), function(z) ev(z)$value,
        function(z) ev(z)$gradient, method = "BFGS", control = list(reltol = 1e-12))
      z <- opt$par - solve(ev(opt$par)$information, ev(opt$par)$gradient)
      C <- f$comparison_review$target_map[2:3, kept, drop = FALSE]
      expect_equal(tab$Estimate[2:3], drop(C %*% z), tolerance = 1e-7)
      expect_equal(f$covariance[2:3, 2:3], C %*% solve(ev(z)$information, t(C)),
        ignore_attr = TRUE, tolerance = 1e-7)
      expect_equal(f$optimization$opt$value, ev(z)$value, tolerance = 1e-10)
      # A conditional null direction changes the unavailable target but no
      # conditional probability, even though unconditional cell logits change.
      null <- svd(D, nv = ncol(D))$v[, qrD$rank + 1L]
      expect_gt(abs(drop(f$comparison_review$target_map[1, ] %*% null)), .1)
      for (b in blocks) expect_lt(diff(range(drop(b$statistics %*% null))), 1e-10)
    }
  }
})

test_that("CML review uses actual totals and unavailable fits retain their reasons", {
  skip_if_not_installed("lpSolve")
  a <- cml_comparison_args()
  a$data$Score <- ifelse(a$data$Group == "A", 0L, 2L)
  a$comparison_review$comparison_review$source$rows <- mfrmr_comparison_rows(transform(a$data, Weight = 1), a$facets)
  f <- do.call(fit_mfrm, a)
  expect_true(all(f$comparisons$ConditionalState == "target_unidentified"))
  expect_true(all(is.na(f$comparisons$Estimate)))
  expect_null(f$optimization)
  a <- cml_comparison_args(); a$maxit <- 1L
  f <- do.call(fit_mfrm, a)
  expect_true(all(is.na(f$comparisons$Estimate)))
  expect_true(all(f$comparisons$NumericalState == "review"))
  expect_false(any(f$comparisons$FormalInferenceEligible))
  a <- cml_comparison_args()
  a$data$Score[a$data$Rater == "R1"] <- 2L
  a$comparison_review$comparison_review$source$rows <- mfrmr_comparison_rows(transform(a$data, Weight = 1), a$facets)
  f <- do.call(fit_mfrm, a)
  expect_true(all(f$comparisons$BoundaryState == "recession_detected"))
  expect_true(all(is.na(f$comparisons$SE)))
  a <- cml_comparison_args(); a$population_formula <- ~1
  expect_error(do.call(fit_mfrm, a), "does not accept")
  a <- cml_comparison_args(); a$comparison_review$comparison_review$specification$population_formula <- ~1
  expect_error(do.call(fit_mfrm, a), "no population model")
  a <- cml_comparison_args(); a$data$Score[1] <- 1L
  expect_error(do.call(fit_mfrm, a), "retained responses")
})

test_that("retained incomplete school examples and reversed origins give reproducible CML inference", {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  skip_if_not_installed("lpSolve")
  cases <- readRDS(test_path("fixtures", "partial-comparisons.rds"))
  for (name in c("RSM-JML", "PCM-JML")) {
    x <- cases[[name]]; r <- x$review; r$comparison_review$specification$method <- "CML"
    a <- x$review$comparison_review$fit_arguments
    a$method <- "CML"; a$population_formula <- NULL
    a <- c(list(data = x$data, fit_scope = "comparisons", comparison_review = r,
      interaction_policy = "silent"), a)
    f <- do.call(fit_mfrm, a)
    expect_identical(f$comparisons$FormalInferenceEligible, c(FALSE, TRUE))
    expect_gt(f$comparisons$SE[2], 0)
    a$data <- a$data[rev(seq_len(nrow(a$data))), ]
    a$comparison_review$comparison_review$specification$groups <- c("B", "A")
    a$comparison_review$comparison_review$target_map[,] <- 0
    reverse <- do.call(fit_mfrm, a)
    expect_equal(reverse$comparisons$Estimate, -f$comparisons$Estimate, tolerance = 1e-7)
    expect_equal(reverse$comparisons$Lower, -f$comparisons$Upper, tolerance = 1e-7)
    expect_equal(reverse$covariance, f$covariance, tolerance = 1e-7)
  }
})


test_that("comparison identification is invariant to tiny and large contrast units", {
  a <- cml_comparison_args()
  prep <- prepare_mfrm_data(a$data, a$person, a$facets, a$score,
    a$rating_min, a$rating_max, keep_original = TRUE)
  spec <- a$comparison_review$comparison_review$specification
  base <- mfrmr_build_comparison_review(prep, spec)
  for (scale in c(1e-200, 1e200)) {
    other <- spec; other$contrasts <- spec$contrasts * scale
    review <- mfrmr_build_comparison_review(prep, other)
    expect_identical(review$table$ConditionalState, base$table$ConditionalState)
    expect_equal(review$table$NullFraction, base$table$NullFraction, tolerance = 1e-12)
  }
  other <- spec; other$contrasts[,] <- 1e308
  expect_error(mfrmr_build_comparison_review(prep, other), "summing to zero")
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    other <- spec; other$contrasts <- spec$contrasts * 1e-200
    a$comparison_review <- mfrmr_build_comparison_review(prep, other)
    expect_error(do.call(fit_mfrm, a), "positive representable variance")
  }
})
