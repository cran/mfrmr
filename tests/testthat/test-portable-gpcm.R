# Fixed synthetic coordinates; no training responses or optimizer runs needed.
portable_gpcm_fixture <- function(owner) {
  readRDS(test_path("fixtures", paste0("calibration-gpcm-", owner, ".rds")))
}

portable_gpcm_oracle <- function(x, rows, interval, prior = list(mean = x$scoring_basis$prior_mean, sd = x$scoring_basis$prior_sd)) {
  coords <- x$parameters$coordinates
  mu <- prior$mean; sd <- prior$sd
  # Fixed row coordinates need resolving only once, independently of theta.
  responses <- lapply(seq_len(nrow(rows)), function(i) {
    location <- sum(vapply(x$model$facet_names, function(f) {
      x$model$facet_signs[[f]] * coords$Value[coords$ParameterClass == 'facet' &
        coords$OwnerFacet == f & coords$Level == rows[[f]][i]]
    }, numeric(1)))
    a <- coords$Value[coords$ParameterClass == 'slope' &
      coords$Level == rows[[x$model$slope_owner]][i]]
    st <- coords[coords$ParameterClass == 'owned_step' &
      coords$Level == rows[[x$model$step_owner]][i], ]
    list(location = location, slope = a, steps = st$Value[order(as.integer(st$Step))],
      score = match(rows$Score[i], x$response$score_map$OriginalScore))
  })
  likelihood <- function(theta) vapply(theta, function(t) {
    logp <- vapply(responses, function(r) {
      logits <- c(0, cumsum(r$slope * (t + r$location - r$steps)))
      logits[r$score] - max(logits) - log(sum(exp(logits - max(logits))))
    }, numeric(1))
    exp(sum(logp)) * dnorm(t, mu, sd)
  }, numeric(1))
  int <- function(fun, lo = -Inf, hi = Inf) integrate(fun, lo, hi,
    rel.tol = 1e-10, subdivisions = 1000L)$value
  mass <- int(likelihood)
  mean <- int(function(t) t * likelihood(t)) / mass
  variance <- int(function(t) (t - mean)^2 * likelihood(t)) / mass
  c(Estimate = mean, SD = sqrt(variance),
    IntervalMass = int(likelihood, interval[1], interval[2]) / mass)
}

test_that("qualified GPCM artifacts retain owners, evidence and native scores on replay", {
  for (owner in c("shared", "separate")) {
    f <- portable_gpcm_fixture(owner); x <- f$calibration
    expect_identical(x$header$schema_version, 2L)
    expect_identical(x$model$step_owner == x$model$slope_owner, owner == "shared")
    expect_identical(x$eligibility$source_readiness_status, "conditional")
    expect_false(x$eligibility$source_scoring_evidence$inference_ready)
    expect_identical(x$eligibility$source_scoring_evidence$source_audit_states[["boundary"]], "not_evaluated")
    expect_length(mfrmr:::mfrmr_calibration_find_prohibited(x), 0L)
    expect_equal(nrow(review_mfrm_calibration(x)), 0)
    scores <- score_mfrm_calibration(x, f$rows)
    cols <- c("Estimate", "SD", "Lower", "Upper")
    expect_equal(unname(as.matrix(scores$estimates[cols])),
      unname(as.matrix(f$native[cols])), tolerance = 1e-8)
    expect_true(all(scores$estimates$ScoreIntegrationReady))
    expect_true(any(scores$estimates$Disposition == "scored"))
    expect_false(any(grepl("GPCM_SCORING_REQUIRES_REVIEW", scores$estimates$ReasonCodes)))
    for (i in seq_len(nrow(scores$estimates))) {
      e <- scores$estimates[i, ]
      ref <- portable_gpcm_oracle(x, f$rows[f$rows$Person == e$Person, ], c(e$Lower, e$Upper))
      expect_equal(e$Estimate, unname(ref["Estimate"]), tolerance = 1e-6)
      expect_equal(e$SD, unname(ref["SD"]), tolerance = 1e-6)
      expect_equal(unname(ref["IntervalMass"]), .95, tolerance = 1e-7)
    }
    path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
    save_mfrm_calibration(x, path)
    expect_identical(load_mfrm_calibration(path), x)
    expect_identical(score_mfrm_calibration(load_mfrm_calibration(path), f$rows), scores)
    expect_s3_class(summary(scores), "summary.mfrm_calibration_score")
    expect_match(paste(capture.output(print(x)), collapse = " "), "Scoring prior: normal")
    for (type in c("interval", "precision", "edge_mass")) {
      figure <- plot(scores, type = type, draw = FALSE)
      expect_s3_class(figure, "mfrm_plot_data")
      expect_s3_class(as_ggplot(figure), "ggplot")
    }
  }
})

test_that("portable scoring priors preserve frozen identity and output assumptions", {
  f <- portable_gpcm_fixture("separate"); x <- f$calibration
  before <- serialize(x, NULL)
  prior <- list(mean = x$scoring_basis$prior_mean + .5, sd = 1.5 * x$scoring_basis$prior_sd)
  scores <- score_mfrm_calibration(x, f$rows, scoring_prior = prior)
  for (i in seq_len(nrow(scores$estimates))) {
    e <- scores$estimates[i, ]
    ref <- portable_gpcm_oracle(x, f$rows[f$rows$Person == e$Person, ], c(e$Lower, e$Upper), prior)
    expect_equal(e$Estimate, unname(ref["Estimate"]), tolerance = 1e-6)
    expect_equal(e$SD, unname(ref["SD"]), tolerance = 1e-6)
    expect_equal(unname(ref["IntervalMass"]), .95, tolerance = 1e-7)
  }
  expect_identical(serialize(x, NULL), before)
  expect_identical(scores$settings$semantic_components, x$integrity$semantic_components)
  expect_true(all(scores$estimates$PriorSource == "user_supplied"))
  summarized <- summary(scores, digits = 0)
  csv <- tempfile(fileext = ".csv"); on.exit(unlink(csv), add = TRUE)
  utils::write.csv(summarized$estimates, csv, row.names = FALSE)
  exported <- utils::read.csv(csv)
  fields <- c("PriorMean", "PriorSD", "RetainedPriorMean", "RetainedPriorSD")
  expect_identical(exported$Person, scores$estimates$Person)
  expect_equal(as.list(exported[fields]), as.list(scores$estimates[fields]), tolerance = 1e-14)
  expect_true(all(summarized$estimates$PriorSD == prior$sd))
  expect_true(all(summarized$estimates$RetainedPriorSD == x$scoring_basis$prior_sd))
  expect_identical(summarized$settings$source_scoring_evidence, x$eligibility$source_scoring_evidence)
  figure <- plot(scores, draw = FALSE)
  expect_identical(figure$data$settings$prior_identity, scores$settings$prior_identity)
  expect_identical(figure$data$settings$source_scoring_evidence, scores$settings$source_scoring_evidence)
  expect_true(all(figure$data$data$PriorSD == prior$sd))
  expect_silent(capture.output(print(summarized)))
  text <- paste(capture.output(print(summarized)), collapse = " ")
  expect_match(text, "user-supplied normal")
  bad <- scores; bad$settings$retained_prior_identity$mean <- 99
  expect_error(summary(bad), "prior or numerical-check records")
  expect_error(plot(bad, draw = FALSE), "prior or numerical-check records")
  bad <- scores; bad$settings$source_scoring_evidence <- NULL
  expect_error(summary(bad), "prior or numerical-check records")
  bad <- summarized; bad$settings$score_integration_review <- NULL
  expect_error(print(bad), "prior or numerical-check records")
  empty <- f$rows; empty$Score <- NA_real_
  omitted <- score_mfrm_calibration(x, empty, missing_response = "omit")
  expect_equal(nrow(omitted$estimates), 0)
  expect_true(all(omitted$person_dispositions$Disposition == "not_scored"))
  expect_s3_class(summary(omitted), "summary.mfrm_calibration_score")
})

test_that("invalid evidence, parameters and coarse scoring do not become portable qualifications", {
  f <- portable_gpcm_fixture("separate"); x <- f$calibration
  mutations <- list(
    function(z) {z$model$slope_owner <- "Unknown"; z},
    function(z) {z$model$slope_action <- "ability_only"; z},
    function(z) {z$parameters$coordinates$Value[z$parameters$coordinates$ParameterClass == "slope"][1] <- 0; z},
    function(z) {z$scoring_basis$prior_sd <- 1e300; z},
    function(z) {z$eligibility$source_readiness_status <- "eligible"; z},
    function(z) {z$eligibility$source_scoring_evidence$local_calibration_review$eligible <- FALSE; z},
    function(z) {z$eligibility$source_scoring_evidence$local_calibration_review <- 1; z},
    function(z) {z$eligibility$source_scoring_evidence$local_calibration_review$integration$NLLChange <- .01; z},
    function(z) {z$eligibility$source_scoring_evidence$source_audit_states["boundary"] <- "has_exclusions"; z})
  for (mutate in mutations) {
    bad <- mutate(x)
    bad$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(bad)
    expect_gt(nrow(review_mfrm_calibration(bad)), 0)
    expect_error(score_mfrm_calibration(bad, f$rows), class = "mfrm_calibration_error")
  }
  changed <- x; changed$scoring_basis$prior_mean <- changed$scoring_basis$prior_mean + 1
  expect_true("IDENTITY_COMPONENT_MISMATCH" %in% review_mfrm_calibration(changed)$Code)
  coarse <- x
  quad <- mfrmr:::gauss_hermite_normal(2L)
  coarse$scoring_basis$quadrature_order <- 2L
  coarse$scoring_basis$nodes <- sort(quad$nodes)
  coarse$scoring_basis$weights <- quad$weights[order(quad$nodes)]
  coarse$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(coarse)
  expect_equal(nrow(review_mfrm_calibration(coarse)), 0)
  expect_error(score_mfrm_calibration(coarse, f$rows), "integration did not pass")
  expect_error(score_mfrm_calibration(x, f$rows, scoring_prior = list(mean = 0, sd = -1)), "positive finite")
  weighted <- f$rows; weighted$Weight <- 2
  expect_error(score_mfrm_calibration(x, weighted, weight = "Weight"), "unit observation weights")
})

test_that("public extraction qualifies native GPCM without a blanket review waiver", {
  f <- readRDS(test_path("fixtures", "mfrm-conditional-scoring-gpcm.rds"))
  # Version 2 compares stored values without ALTREP materialization state.
  before <- serialize(f$fit, NULL, version = 2)
  draft <- extract_mfrm_calibration(f$fit)
  expect_identical(draft$eligibility$source_readiness_status, "conditional")
  expect_identical(serialize(f$fit, NULL, version = 2), before)
  expect_equal(nrow(review_mfrm_calibration(draft)), 0)
  failed <- f$fit; failed$readiness$fit$NumericalState <- "failed"
  expect_error(extract_mfrm_calibration(failed), "passing conditional calibration checks")
  failed <- f$fit; failed$population$sigma2 <- 0
  expect_error(extract_mfrm_calibration(failed), "passing conditional calibration checks")
})

test_that("GPCM scores replay in a new R process with only artifact and scoring rows", {
  f <- portable_gpcm_fixture("separate")
  dir <- tempfile(); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  calibration <- file.path(dir, "calibration.rds"); rows <- file.path(dir, "rows.rds")
  output <- file.path(dir, "scores.rds"); worker <- file.path(dir, "worker.R")
  save_mfrm_calibration(f$calibration, calibration); saveRDS(f$rows, rows)
  root <- normalizePath(find.package("mfrmr"))
  writeLines(c("args <- commandArgs(TRUE)", "if (file.exists(file.path(args[1], 'R', 'api-calibration.R'))) pkgload::load_all(args[1], quiet=TRUE, compile=FALSE) else library(mfrmr, lib.loc=args[5])",
    "x <- load_mfrm_calibration(args[2])", "s <- score_mfrm_calibration(x, readRDS(args[3]))",
    "saveRDS(list(scores=s, package_path=normalizePath(find.package('mfrmr'))), args[4])"), worker)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(worker), shQuote(root), shQuote(calibration), shQuote(rows), shQuote(output), shQuote(dirname(find.package("mfrmr")))),
    stdout = TRUE, stderr = TRUE)
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  expect_true(file.exists(output))
  if (file.exists(output)) {
    replay <- readRDS(output)
    expect_identical(replay$package_path, root)
    expect_identical(replay$scores, score_mfrm_calibration(f$calibration, f$rows))
  }
})

test_that("adaptive score summaries distinguish reported-score checks and preserve unscored people", {
  f <- portable_gpcm_fixture("separate")
  rows <- f$rows
  missing <- rows[1, , drop = FALSE]; missing$Person <- "NA"; missing$Score <- NA
  rows <- rbind(rows, missing)
  scores <- score_mfrm_calibration(f$calibration, rows, missing_response = "omit")
  s <- summary(scores)
  text <- paste(capture.output(print(s)), collapse = " ")
  expect_match(text, "Reported-score integration check")
  expect_match(text, "4 of 4 scored Persons passed")
  expect_match(text, "fixed grid in this table did not produce")
  expect_true(all(s$settings$score_integration_review$Passed))
  # The unrelated fixed grid can disagree while the actually reported adaptive
  # scores agree with their references. Do not present these as one decision.
  expect_gt(max(s$quadrature_overview$MaxAbsSDChange), 1e-5)
  expect_identical(s$review$Disposition[s$review$Person == "NA"], "not_scored")
  figure <- plot(scores, draw = FALSE)
  expect_true("NA" %in% figure$data$unplotted_dispositions$Person)
  expect_match(figure$data$subtitle, "1 not scored")
  expect_identical(figure$data$reference_lines$label, "Scale origin")
  empty <- score_mfrm_calibration(f$calibration, missing, missing_response = "omit")
  expect_match(paste(capture.output(print(summary(empty))), collapse = " "), "No scored Persons")
})

portable_gpcm_jml_fixture <- local({
  cache <- list()
  function(owner = "Criterion") {
    if (is.null(cache[[owner]])) {
      skip_if_not_installed("lpSolve")
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:14], ]
      fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
        method = "JML", model = "GPCM", step_facet = owner, slope_facet = owner,
        maxit = 600, reltol = 1e-11))
      rows <- expand.grid(Person = c("001", "1", "NA", "sparse", "missing"),
        Rater = unique(d$Rater), Criterion = unique(d$Criterion), stringsAsFactors = FALSE)
      rows$Score <- 3
      rows$Score[rows$Person == "001"] <- min(d$Score)
      rows$Score[rows$Person == "1"] <- max(d$Score)
      rows$Score[rows$Person == "missing"] <- NA_real_
      sparse <- which(rows$Person == "sparse"); rows <- rows[-sparse[-1], ]
      cache[[owner]] <<- list(fit = fit, rows = rows)
    }
    cache[[owner]]
  }
})

test_that("GPCM JML portable EAP retains incomplete global audits and reference priors", {
  for (owner in c("Criterion", "Rater")) {
    f <- portable_gpcm_jml_fixture(owner)
    before <- serialize(f$fit, NULL)
    a <- freeze_mfrm_calibration(validate_mfrm_calibration(
      extract_mfrm_calibration(f$fit, scoring_quad_points = 141)))
    expect_identical(a$header$schema_version, 4L)
    expect_identical(a$model$slope_owner, owner)
    expect_identical(a$scoring_basis$type, "post_hoc_standard_normal")
    evidence <- a$eligibility$source_scoring_evidence
    expect_false(evidence$inference_ready)
    expect_identical(unname(evidence$source_audit_states[1:2]), rep("not_evaluated", 2))
    expect_identical(evidence$local_calibration_review$basis, "local_jml_and_curvature_v1")
    expect_false(any(c("integration", "covariance") %in% names(evidence$local_calibration_review)))
    expect_error(confint(f$fit, parm = "slopes"), "MML")
    for (prior in list(NULL, list(mean = .25, sd = 1.2))) {
      scored <- score_mfrm_calibration(a, f$rows, missing_response = "omit", scoring_prior = prior)
      native <- predict_mfrm_units(f$fit, f$rows[!is.na(f$rows$Score), ],
        scoring_quad_points = 141, scoring_prior = prior)
      fields <- c("Estimate", "SD", "Lower", "Upper")
      pos <- match(scored$estimates$Person, native$estimates$Person)
      expect_equal(unname(as.matrix(scored$estimates[fields])),
        unname(as.matrix(native$estimates[pos, fields])), tolerance = 1e-10)
      expect_identical(native$settings$source_scoring_status, "conditional")
      expect_false(native$settings$source_inference_ready)
      expect_true(all(scored$estimates$ScoreIntegrationReady))
      expect_s3_class(summary(native), "summary.mfrm_unit_prediction")
      broken <- native; broken$settings$local_calibration_review$checks$known_boundary_absent <- FALSE
      expect_error(summary(broken), "inconsistent numerical")
      expect_s3_class(summary(scored), "summary.mfrm_calibration_score")
      expect_identical(scored$settings$source_scoring_evidence, evidence)
      expect_true(any(grepl("remain incomplete", capture.output(print(scored)), fixed = TRUE)))
      expect_identical(scored$person_dispositions$Disposition[scored$person_dispositions$Person == "missing"], "not_scored")
      expect_true(all(scored$person_dispositions$Disposition[scored$person_dispositions$Person %in% c("001", "1", "sparse")] == "scored_review"))
      target <- scored$estimates[scored$estimates$Person == "NA", ]
      oracle <- portable_gpcm_oracle(a, f$rows[f$rows$Person == "NA", ],
        unlist(target[c("Lower", "Upper")]), if (is.null(prior)) list(mean = 0, sd = 1) else prior)
      expect_equal(as.numeric(target[1, c("Estimate", "SD")]), unname(oracle[1:2]), tolerance = 1e-6)
      expect_equal(unname(oracle["IntervalMass"]), .95, tolerance = 1e-7)
    }
    expect_identical(serialize(f$fit, NULL), before)
  }
})

test_that("GPCM JML joint derivatives and local curvature match an independent likelihood", {
  f <- portable_gpcm_jml_fixture()$fit
  cfg <- f$config; sizes <- mfrmr:::build_param_sizes(cfg)
  idx <- mfrmr:::build_indices(f$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
  evaluator <- mfrmr:::make_mfrm_direct_evaluator("JML",
    mfrmr:::make_param_cache(sizes, cfg, idx, is_mml = FALSE), idx, cfg, sizes)
  d <- f$prep$data
  oracle <- function(par) {
    p <- mfrmr:::expand_params(par, sizes, cfg)
    location <- p$theta[as.integer(d$Person)]
    for (facet in cfg$facet_names) location <- location + cfg$facet_signs[[facet]] * p$facets[[facet]][as.integer(d[[facet]])]
    owner <- as.integer(d[[cfg$step_facet]])
    -sum(vapply(seq_len(nrow(d)), function(i) {
      z <- c(0, cumsum(p$slopes[owner[i]] * (location[i] - p$steps_mat[owner[i], ])))
      z[d$score_k[i] + 1L] - max(z) - log(sum(exp(z - max(z))))
    }, numeric(1)))
  }
  par <- f$opt$par + .02 * sin(seq_along(f$opt$par))
  expect_equal(evaluator$value(par), oracle(par), tolerance = 1e-10)
  gradient <- vapply(seq_along(par), function(i) {
    lo <- hi <- par; lo[i] <- lo[i] - 1e-5; hi[i] <- hi[i] + 1e-5
    (oracle(hi) - oracle(lo)) / 2e-5
  }, numeric(1))
  expect_equal(unname(evaluator$gradient(par)), gradient, tolerance = 1e-6)
  H <- stats::optimHess(f$opt$par, evaluator$value, evaluator$gradient)
  expect_gt(min(eigen(H, symmetric = TRUE, only.values = TRUE)$values), 0)
  for (i in 1:4) {
    v <- sin(seq_along(par) * i); v <- v / sqrt(sum(v^2)); h <- 1e-3
    curvature <- (oracle(f$opt$par + h*v) - 2*oracle(f$opt$par) + oracle(f$opt$par - h*v)) / h^2
    expect_equal(drop(crossprod(v, H %*% v)), curvature, tolerance = 1e-5)
  }
})

test_that("GPCM JML local checks refuse unresolved sources and singular curvature", {
  f <- portable_gpcm_jml_fixture()$fit
  bad <- f; bad$opt$par[1] <- bad$opt$par[1] + 1
  expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  for (field in c("EstimabilityState", "BoundaryState", "NumericalState")) {
    bad <- f; bad$readiness$fit[[field]] <- "failed"
    expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  }
  bad <- f; bad$config$slope_facet <- "Rater"
  expect_error(extract_mfrm_calibration(bad), "same slope and step owner")
  bad <- f; bad$config$theta_spec$anchors[1] <- 0
  expect_error(extract_mfrm_calibration(bad), "MODEL_STRUCTURE_UNSUPPORTED")
  bad <- f; bad$config$boundary_audit$joint_additive <- NULL
  expect_error(extract_mfrm_calibration(bad), "SOURCE_READINESS_INELIGIBLE")
  a <- extract_mfrm_calibration(f)
  for (change in c("curvature", "boundary", "prior", "owner", "estimator")) {
    bad <- a
    if (change == "curvature") bad$eligibility$source_scoring_evidence$local_calibration_review$checks$positive_definite <- FALSE
    if (change == "boundary") bad$eligibility$source_scoring_evidence$local_calibration_review$checks$known_boundary_absent <- FALSE
    if (change == "prior") bad$scoring_basis$type <- "frozen_estimated_normal"
    if (change == "owner") bad$model$slope_owner <- "Rater"
    if (change == "estimator") bad$model$estimator <- "MML"
    bad$integrity$semantic_components <- mfrmr:::mfrmr_calibration_semantic_components(bad)
    expect_gt(nrow(review_mfrm_calibration(bad)), 0)
    expect_error(validate_mfrm_calibration(bad), class = "mfrm_calibration_error")
  }
  local_mocked_bindings(optimHess = function(par, ...) diag(0, length(par)), .package = "mfrmr")
  expect_false(mfrmr:::prediction_source_scoring_readiness(f)$ready)
  expect_error(extract_mfrm_calibration(f), "SOURCE_READINESS_INELIGIBLE")
})

test_that("GPCM JML source boundary is not qualified by finite optimizer traces", {
  skip_if_not_installed("lpSolve")
  d <- load_mfrmr_data("example_core"); d <- d[d$Person %in% unique(d$Person)[1:14], ]
  d$Score[d$Person == unique(d$Person)[1]] <- max(d$Score)
  f <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    method = "JML", model = "GPCM", step_facet = "Criterion", maxit = 80))
  expect_true(any(is.infinite(f$facets$person$PrimaryEstimate)))
  expect_false(mfrmr:::prediction_source_scoring_readiness(f)$ready)
  expect_error(extract_mfrm_calibration(f), "SOURCE_READINESS_INELIGIBLE")
})

test_that("GPCM JML artifacts and saved scores replay without fitting", {
  skip_if_not_installed("callr")
  root <- normalizePath(find.package("mfrmr"))
  for (owner in c("Criterion", "Rater")) {
    f <- portable_gpcm_jml_fixture(owner)
    a <- freeze_mfrm_calibration(validate_mfrm_calibration(extract_mfrm_calibration(f$fit, scoring_quad_points = 141)))
    path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
    save_mfrm_calibration(a, path)
    scores <- score_mfrm_calibration(a, f$rows, missing_response = "omit")
    worker <- function(root, path, rows) {
      if (file.exists(file.path(root, "R", "api-calibration.R"))) pkgload::load_all(root, quiet=TRUE, compile=FALSE) else library("mfrmr",lib.loc=dirname(root))
      testthat::local_mocked_bindings(fit_mfrm=function(...) stop("unexpected fit"),
        predict_mfrm_units=function(...) stop("unexpected native scorer"),
        prediction_jml_calibration_review=function(...) stop("unexpected source evaluation"), .package="mfrmr")
      s <- mfrmr::score_mfrm_calibration(mfrmr::load_mfrm_calibration(path), rows, missing_response="omit")
      save <- tempfile(); on.exit(unlink(save), add = TRUE); saveRDS(s,save); s<-readRDS(save)
      p <- plot(s,draw=FALSE,main="New-cohort EAP")
      if(requireNamespace("ggplot2",quietly=TRUE)) ggplot2::ggplot_build(mfrmr::as_ggplot(s))
      list(scores=s,summary=summary(s),plot=p$data)
    }
    environment(worker)<-baseenv()
    replay <- callr::r(worker,args=list(root,path,f$rows),libpath=.libPaths())
    expect_identical(replay$scores, scores)
    expect_identical(replay$summary, summary(scores))
    expect_identical(replay$plot, plot(scores,draw=FALSE,main="New-cohort EAP")$data)
  }
})
