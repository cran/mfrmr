predictive_comparison_fixture <- local({
  saved <- NULL
  function() {
    if (!is.null(saved)) return(saved)
    d <- expand.grid(Person = paste0("P", 1:12), Rater = c("A", "B"), Criterion = c("X", "Y"))
    d$Score <- c(0,1,1,0,1,0,1,1,0,0,1,0, 1,0,0,1,0,1,1,0,1,0,0,1,
      1,0,1,0,1,0,1,0,1,1,0,1, 0,1,0,1,0,0,1,1,0,1,0,0)
    o <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", quad_points = 61)
    saved <<- list(ordinary = o, data = d)
    saved
  }
})

predictive_zero_testlet <- function(o, data) {
  a <- mfrm_ordinary_response_input(o)
  input <- mfrm_testlet_data(data, "Person", "Score", "Rater", c("Rater", "Criterion"), a$score_levels, "omit")
  beta <- as.vector(qr.solve(input$X, a$offset))
  par <- c(beta, a$steps - a$mean, 0, a$sd^2)
  fixed <- o$facets$others[c("Facet", "Level", "Estimate")]
  fixed$Parameter <- "Fixed facet"
  structure(list(input = input, parameters = par,
    calibration = list(beta = beta, steps = a$steps - a$mean, variance = 0,
      person_variance = a$sd^2, person_sd = a$sd),
    calibration_table = fixed, checks = data.frame(NumericalReady = TRUE, InformationPositive = TRUE),
    loglik = o$summary$LogLik,
    settings = list(quad_points = 61L, fixed_person_sd = if (!isTRUE(o$population$active)) 1,
      method = "Fixed reference calibration")), class = "mfrm_testlet")
}

test_that("ordinary response moments agree with independent continuous integration", {
  f <- predictive_comparison_fixture(); o <- f$ordinary
  for (estimated in c(FALSE, TRUE)) {
    if (estimated) {
      o$population <- list(active = TRUE, design_columns = "(Intercept)", coefficients = .6, sigma2 = .81)
      o$config$population_spec$design_matrix <- matrix(1, 12, 1)
    }
    a <- mfrm_ordinary_response_input(o)
    out <- mfrm_response_diagnostics(o, rows = c(1, 13), group_by = "Rater")
    ids <- which(a$person == a$person[1])
    # Deliberately no package response kernel or quadrature in this oracle.
    integral <- function(replicate = FALSE, within = FALSE) integrate(function(theta) vapply(theta, function(t) {
      p <- plogis(t - a$offset[ids] - a$steps)
      likelihood <- prod(ifelse(a$y[ids] == 1, p, 1 - p))
      likelihood * dnorm(t, a$mean, a$sd) * if (within) p[1] * (1 - p[1]) else if (replicate) p[1] else 1
    }, numeric(1)), -Inf, Inf, rel.tol = 1e-10)$value
    mass <- integral(); mu <- integral(TRUE) / mass
    expect_equal(out$probabilities[1, 2], unname(mu), tolerance = 1e-8)
    expect_equal(out$rows$PredictiveVariance[1], mu * (1 - mu), tolerance = 1e-8)
    expect_gt(out$rows$PredictiveVariance[1], integral(within = TRUE) / mass)
    all <- mfrm_response_diagnostics(o, group_by = "Rater")
    expect_equal(out$probabilities, all$probabilities[c(1,13), ])
    e <- predictive_zero_testlet(o, f$data)
    other <- mfrm_response_diagnostics(e, group_by = "Rater")
    expect_equal(other$probabilities, all$probabilities, tolerance = 1e-8)
  }
})

test_that("original categories, omitted positions and duplicated events survive matching", {
  f <- predictive_comparison_fixture(); d <- f$data
  d <- rbind(d, d[1, ], d[1, ]); d$Score[c(1, 49)] <- NA
  o <- f$ordinary; o$prep <- prepare_mfrm_data(d, "Person", c("Rater", "Criterion"), "Score")
  e <- predictive_zero_testlet(o, d)
  a <- mfrm_response_diagnostics(o, group_by = "Rater")
  b <- mfrm_response_diagnostics(e, group_by = "Rater")
  expect_equal(a$rows$InputRow, seq_len(50))
  expect_equal(which(a$rows$Status == "missing_score"), c(1,49))
  x <- compare_mfrm(o, e, response_diagnostics = list(a, b))
  expect_equal(nrow(x$responses$rows), 50)
  expect_equal(sum(x$responses$rows$Status == "missing_score"), 2)
  expect_equal(x$responses$measures$Missing, c(2,0))
  expect_true(all(is.na(x$responses$probabilities$Difference[x$responses$probabilities$Status == "missing_score"])))
  # Reversing source order also moves missing/duplicate input-row identities.
  rev <- o; rev$prep <- prepare_mfrm_data(d[50:1, ], "Person", c("Rater", "Criterion"), "Score")
  aa <- mfrm_response_diagnostics(rev, group_by = "Rater")
  y <- compare_mfrm(rev, e, response_diagnostics = list(aa, b))
  expect_equal(y$responses$rows$ExpectedScoreReference, x$responses$rows$ExpectedScoreReference)
  expect_equal(y$responses$probabilities$Difference, x$responses$probabilities$Difference)
  expect_error(compare_mfrm(rev,e,response_diagnostics=list(a,b)), "exact source roster")
  shifted <- o; shifted$prep$score_map$OriginalScore <- shifted$prep$score_map$OriginalScore + 5
  moved <- mfrm_response_diagnostics(shifted, group_by = "Rater")
  expect_equal(moved$rows$Score, a$rows$Score + 5)
  expect_equal(moved$rows$ExpectedScore, a$rows$ExpectedScore + 5)
  expect_equal(moved$rows$PredictiveVariance, a$rows$PredictiveVariance)
  expect_equal(moved$measures$Infit, a$measures$Infit)
  old <- o; old$prep$omitted_input_rows <- NULL
  expect_error(mfrm_response_diagnostics(old), "older fit")
})

test_that("predictive comparisons check selected events, target and source calibration", {
  f <- predictive_comparison_fixture(); o <- f$ordinary; e <- predictive_zero_testlet(o, f$data)
  a <- mfrm_response_diagnostics(o, rows = 1:12, group_by = "Rater")
  b <- mfrm_response_diagnostics(e, rows = 1:12, group_by = "Rater")
  x <- compare_mfrm(o,e,response_diagnostics=list(a,b))
  reverse <- compare_mfrm(e,o,response_diagnostics=list(b,a))
  expect_equal(reverse$responses$rows$ExpectedScoreDifference, -x$responses$rows$ExpectedScoreDifference)
  expect_equal(reverse$responses$measures$InfitDifference, -x$responses$measures$InfitDifference)
  expect_equal(rowsum(x$responses$probabilities$Difference, x$responses$probabilities$MatchedRow), matrix(0,12,1), ignore_attr=TRUE, tolerance=1e-12)
  expect_error(compare_mfrm(o,e,response_diagnostics=list(a,mfrm_response_diagnostics(e,rows=2:13,group_by="Rater"))), "Selected rating events")
  bad <- b; bad$settings$target <- "Another target"
  expect_error(compare_mfrm(o,e,response_diagnostics=list(a,bad)), "same posterior predictive target")
  bad <- b; bad$settings$group_by <- "Person"
  expect_error(compare_mfrm(o,e,response_diagnostics=list(a,bad)), "group_by")
  bad <- o; bad$steps$Estimate <- bad$steps$Estimate + .1
  expect_error(compare_mfrm(bad,e,response_diagnostics=list(a,b)), "matching calibration")
  expect_error(compare_mfrm(o,o,response_diagnostics=list(a,a)), "one ordinary RSM MML")
  # A failed integral cannot be silently removed from a row or group difference.
  bad <- b; bad$rows$Status[1] <- "unavailable"; bad$rows$Reason[1] <- "Unresolved integral"
  bad$probabilities[1,] <- NA; bad$measures$Status <- "unavailable"; bad$measures$Available <- 11L
  failed <- compare_mfrm(o,e,response_diagnostics=list(a,bad))
  expect_true(is.na(failed$responses$rows$ExpectedScoreDifference[1]))
  expect_true(all(is.na(failed$responses$measures$InfitDifference)))
  expect_equal(failed$responses$measures$ComparisonAvailable, 11)
})

test_that("predictive plots and reports reuse saved quantities without cutoffs", {
  f <- predictive_comparison_fixture(); o <- f$ordinary; e <- predictive_zero_testlet(o, f$data)
  a <- mfrm_response_diagnostics(o, group_by="Rater"); b <- mfrm_response_diagnostics(e, group_by="Rater")
  x <- compare_mfrm(o,e,response_diagnostics=list(a,b))
  local_mocked_bindings(mfrm_response_diagnostics=function(...) stop("Unexpected integration"),
    mfrm_ordinary_response_probabilities=function(...) stop("Unexpected integration"),
    mfrm_testlet_response_probabilities=function(...) stop("Unexpected integration"), .package="mfrmr")
  for (metric in c("expected_score", "variance", "infit", "outfit", "probability")) {
    p <- plot(x, metric=metric, style="difference", draw=FALSE, show_labels=FALSE, show_title=FALSE, show_notes=FALSE)
    expect_equal(p$data$table$Y, p$data$table$Difference)
    expect_true(all(is.finite(p$data$table$Y)))
    expect_s3_class(as_ggplot(p), "ggplot")
    expect_match(p$data$alt_text, "no model ranking or fit cutoffs")
  }
  expect_error(plot(x,metric="variance",category=0,draw=FALSE), "only to probability")
  expect_error(plot(x,metric="probability",category=9,draw=FALSE), "numeric score")
  p <- plot(x,metric="probability",category=1,draw=FALSE)
  expect_equal(unique(p$data$table$Score), 1)
  expect_equal(p$data$xlim, c(0,1))
  expect_false(p$data$display$show_labels)
  r <- mfrm_results(o,response_diagnostics=a,compute="never")
  expect_identical(r$response_diagnostics, a)
  expect_identical(mfrm_report(r)$tables$response_measures, a$measures)
  expect_identical(r$tables$response_overview$Status, "available")
  expect_match(mfrm_report(r)$markdown, r$tables$response_overview$Detail, fixed=TRUE)
  old <- r; old$tables$response_overview <- NULL
  expect_identical(mfrm_report(old)$tables$response_overview, r$tables$response_overview)
  expect_identical(summary(old)$status, summary(r)$status)
  expect_s3_class(plot(r,type="response_diagnostics",draw=FALSE), "mfrm_plot_data")
  er <- mfrm_results(e,response_diagnostics=b,comparison=x,compute="never")
  expect_null(e$calibration_table[["SE"]])
  ci <- confint(e)
  expect_true(all(is.na(ci)))
  expect_true(all(attr(ci, "diagnostics")$IntervalStatus == "invalid_estimate_or_se"))
  expect_false(any(attr(ci, "diagnostics")$CIEligible))
  expect_identical(er$tables$comparison_response_probabilities,x$responses$probabilities)
  expect_identical(mfrm_report(er)$tables$comparison_response_rows,x$responses$rows)
  expect_s3_class(plot(er,type="response_comparison",draw=FALSE), "mfrm_plot_data")
  expect_error(mfrm_results(o,diagnostics=a,response_diagnostics=a), "only once")
  exported <- export_mfrm_results(r,output_dir=tempfile(),include="replay",acknowledge_sensitive=TRUE)
  expect_true(any(exported$written_files$Format == "rds"))
  script <- exported$written_files$Path[exported$written_files$Component == "replay_code"]
  expect_true(any(grepl("readRDS",readLines(script))))
})

test_that("ordinary predictive diagnostics refuse unsupported model and population contracts", {
  o <- predictive_comparison_fixture()$ordinary
  for (model in c("PCM", "GPCM")) {
    bad <- o; bad$config$model <- model
    expect_error(mfrm_response_diagnostics(bad), "RSM MML")
  }
  bad <- o; bad$config$method <- "JML"
  expect_error(mfrm_response_diagnostics(bad), "RSM MML")
  bad <- o; bad$prep$data$Weight[1] <- 2
  expect_error(mfrm_response_diagnostics(bad), "unit weights")
  bad <- o; bad$config$facet_signs[1] <- 1
  expect_error(mfrm_response_diagnostics(bad), "severity facets")
  bad <- o; bad$population <- list(active=TRUE,design_columns=c("(Intercept)","Group"))
  expect_error(mfrm_response_diagnostics(bad), "intercept-only")
  bad <- o; bad$readiness$fit$FitReadiness <- "blocked"
  expect_error(mfrm_response_diagnostics(bad), "blocked source checks")
})
