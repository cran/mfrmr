small_sample_ratings <- function() data.frame(
  Person = rep(paste0("P", 1:4), each = 2), Rater = rep(c("A", "B"), 4),
  Score = c(1, 0, 1, 0, 1, 0, 0, 1))

# Independent continuous N(0,1) integral for the one free severity b, with
# severities (b,-b). Differentiate probabilities inside the integral.
small_sample_normal_reference <- function(b) {
  patterns <- rbind(c(1, 0), c(0, 1)); frequency <- c(3, 1)
  rows <- t(apply(patterns, 1, function(y) {
    integrand <- function(theta, order) {
      a <- plogis(theta - b); other <- plogis(theta + b)
      probability <- (if (y[1] == 1) a else 1-a) * (if (y[2] == 1) other else 1-other)
      score <- -y[1] + y[2] + a - other
      multiplier <- switch(order, 1, score, score^2 - a*(1-a) - other*(1-other))
      probability * multiplier * dnorm(theta)
    }
    value <- vapply(1:3, function(k) integrate(integrand, -Inf, Inf,
      order=k, rel.tol=1e-11, abs.tol=1e-12)$value, numeric(1))
    c(NLL=-log(value[1]), Gradient=-value[2]/value[1],
      Information=(value[2]/value[1])^2-value[3]/value[1])
  }))
  drop(crossprod(frequency, rows))
}

test_that("eight ratings retain identifiable MML and CML points and uncertainty", {
  skip_if_not_installed("lpSolve")
  d <- small_sample_ratings()
  root <- uniroot(function(b) small_sample_normal_reference(b)["Gradient"],
    c(-5, 5), tol=1e-10)$root
  reference <- small_sample_normal_reference(root)
  left <- small_sample_normal_reference(root-1e-4)
  right <- small_sample_normal_reference(root+1e-4)
  expect_equal(unname((right["NLL"]-left["NLL"])/2e-4),
    unname(reference["Gradient"]), tolerance=1e-7)
  expect_equal(unname((right["Gradient"]-left["Gradient"])/2e-4),
    unname(reference["Information"]), tolerance=1e-7)
  for (model in c("RSM", "PCM")) {
    args <- list(data=d, person="Person", facets="Rater", score="Score",
      model=model, rating_min=0, rating_max=1, category_policy="preserve",
      maxit=300, reltol=1e-10)
    if (model == "PCM") args$step_facet <- "Rater"
    fit <- suppressWarnings(do.call(fit_mfrm, c(args, list(method="MML", quad_points=61))))
    expect_equal(nobs(fit, "ratings"), 8)
    expect_equal(unname(coef(fit)), root, tolerance=1e-7)
    expect_equal(fit$opt$value, unname(reference["NLL"]), tolerance=1e-10)
    C <- matrix(c(1,-1), 1, dimnames=list("A-B", c("A","B")))
    expect_warning(ci <- mfrm_facet_intervals(fit, "Rater", C), "occur only once")
    expect_true(all(ci$checks$Passed))
    expect_equal(ci$table$SE, 2/sqrt(unname(reference["Information"])), tolerance=1e-6)
    expect_true(ci$table$CIEligible)
    expect_false(ci$table$CoverageQualified)
    expect_false(fit$readiness$fit$InferenceReady)
    expect_match(ci$table$InferenceCaution, "small samples")
    expect_warning(V <- vcov(fit), "occur only once")
    expect_equal(as.numeric(V), 1/unname(reference["Information"]), tolerance=1e-6)
    expect_warning(bounds <- confint(fit), "occur only once")
    expect_equal(unname(bounds[1,]), root+c(-1,1)*qnorm(.975)*sqrt(as.numeric(V)),
      tolerance=1e-7)
    expect_identical(ci$fit, fit)
    result <- mfrm_results(fit, intervals=list(raters=ci), compute="never")
    expect_identical(result$tables$facet_raters_intervals$InferenceCaution, ci$table$InferenceCaution)
    expect_identical(mfrm_report(result)$tables$facet_raters_intervals$SE, ci$table$SE)
    expect_identical(unserialize(serialize(ci, NULL)), ci)
    cml <- suppressWarnings(do.call(fit_mfrm, c(args, list(method="CML"))))
    expect_equal(unname(coef(cml)), -log(3)/2, tolerance=1e-7)
    expect_equal(as.numeric(vcov(cml)), 1/3, tolerance=1e-7)
    expect_equal(mfrm_facet_intervals(cml, "Rater", C)$table$SE, sqrt(4/3), tolerance=1e-7)
  }
})

test_that("sparse counts cannot waive identification, source or information failures", {
  d <- small_sample_ratings()
  f <- suppressWarnings(fit_mfrm(d,"Person","Rater","Score", quad_points=61,
    rating_min=0,rating_max=1,category_policy="preserve",reltol=1e-10))
  for (kind in c("identification", "optimizer", "parameters", "empty_category")) {
    bad <- f
    if (kind == "identification") bad$readiness$fit$EstimabilityState <- "structurally_unidentified"
    if (kind == "optimizer") bad$opt$optimizer_diagnostics$ConvergenceSeverity <- "fail"
    if (kind == "parameters") bad$opt$par <- bad$opt$par+.1
    if (kind == "empty_category") bad$prep$data$Score[bad$prep$data$Rater=="A"] <- 0L
    ci <- suppressWarnings(mfrm_facet_intervals(bad,"Rater"))
    expect_true(all(ci$table$Status=="numerical_review_failed"), info=kind)
    expect_true(all(is.na(ci$table$SE) & is.na(ci$table$Lower)), info=kind)
    expect_false(any(ci$table$CIEligible), info=kind)
    expect_true(all(nzchar(ci$table$InferenceReview)), info=kind)
  }
  local_mocked_bindings(compute_mml_parameter_covariance=function(...)
    list(status="unavailable", detail="Independent information failure"), .package="mfrmr")
  ci <- mfrm_facet_intervals(f,"Rater")
  expect_equal(ci$table$Estimate, unname(f$facets$others$Estimate))
  expect_true(all(is.na(ci$table$SE)))
  expect_match(ci$table$InferenceReview[1], "Independent information failure")
  expect_error(vcov(f), "Joint covariance is unavailable")
})
