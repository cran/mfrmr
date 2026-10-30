mml_starts_fixture <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- load_mfrmr_data("example_operational")
      f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
        mml_starts = c(0, .5, -.5), reltol = 1e-10)
      saved <<- list(data = d, fit = f)
    }
    saved
  }
})

test_that("ordinary starts retain same-objective outcomes and the default result", {
  fx <- mml_starts_fixture(); f <- fx$fit; review <- f$opt$start_review
  original <- fit_mfrm(fx$data, "Person", c("Rater", "Criterion"), "Score", reltol = 1e-10)
  expect_identical(review$attempts$default$par, original$opt$par)
  expect_identical(review$attempts$default$value, original$opt$value)
  expect_equal(sum(review$table$Selected), 1L)
  expect_identical(review$table$Start, c("default", "offset_1", "offset_2"))
  expect_identical(review$parameter_map$OptimizerIndex, seq_along(f$opt$par))
  expect_equal(unname(review$starts$offset_1 - review$starts$default),
    rep(c(.5, -.5), length.out = length(f$opt$par)))
  expect_identical(f$opt$value, min(review$table$NegativeLogLik))
  expect_identical(f$summary$LogLik, -f$opt$value)
  expect_true(all(review$table$LogLikLoss >= 0))
  expect_identical(summary(f)$starting_values, review$table)
  expect_match(paste(capture.output(print(summary(f))), collapse = " "), "Starting-value review")
  expect_null(original$opt$start_review)
  expect_null(original$config$replay_inputs$mml_starts)
})

test_that("named starts, weights and anchors retain the constrained objective", {
  # Full-suite scenario: ordinary checks retain live starts and engine-specific tests.
  skip_on_cran()
  d <- mml_starts_fixture()$data
  d$Weight <- rep(c(.5, 1, 2), length.out = nrow(d))
  anchor <- data.frame(facet = "Rater", level = as.character(d$Rater[1]), anchor = .2)
  f <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    model = "PCM", step_facet = "Criterion", weight = "Weight", anchors = anchor,
    mml_starts = c(0, .5, -.5), reltol = 1e-10))
  r <- f$opt$start_review
  for (attempt in r$attempts) {
    p <- expand_params(attempt$par, build_param_sizes(f$config), f$config)
    expect_equal(unname(p$facets$Rater[as.character(d$Rater[1])]), .2)
    # Independent raw-label PCM likelihood: sum weighted category log
    # probabilities within Person, then integrate once for each Person.
    q <- gauss_hermite_normal(31)
    dat <- f$prep$data
    row_logprob <- vapply(seq_len(nrow(dat)), function(j) {
      base <- sum(vapply(f$config$facet_names, function(facet)
        f$config$facet_signs[[facet]] * p$facets[[facet]][as.character(dat[[facet]][j])], numeric(1)))
      steps <- p$steps_mat[match(as.character(dat$Criterion[j]), f$config$facet_levels$Criterion), ]
      logits <- sweep(outer(q$nodes + base, 0:length(steps)), 2, c(0,cumsum(steps)), "-")
      maximum <- apply(logits, 1, max)
      logits[, dat$Score[j] - f$config$rating_min + 1L] - maximum -
        log(rowSums(exp(logits - maximum)))
    }, numeric(length(q$nodes)))
    marginal <- vapply(split(seq_len(nrow(dat)), dat$Person), function(rows) {
      ll <- drop(row_logprob[,rows,drop=FALSE] %*% dat$Weight[rows]) + log(q$weights)
      max(ll) + log(sum(exp(ll-max(ll))))
    }, numeric(1))
    expect_lt(abs(-sum(marginal) - attempt$value), 1e-8)
  }
  key <- r$parameter_map$Coordinate[1]
  custom <- stats::setNames(.75, key)
  starts <- mfrmr_mml_start_vectors(unname(r$starts$default), list(custom = custom), r$parameter_map)
  expect_equal(starts$custom[[key]], .75)
  expect_identical(starts$custom[-1], starts$default[-1])
  expect_true(all(r$table$LogLikLoss >= 0))
  expect_identical(f$config$replay_inputs$weight, "Weight")
  custom_fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    model = "PCM", step_facet = "Criterion", weight = "Weight", anchors = anchor,
    mml_starts = list(custom = custom), reltol = 1e-10))
  expect_identical(custom_fit$opt$start_review$starts$custom, starts$custom)
  expect_error(resolve_dff_refit_controls(custom_fit), "subgroup linking changes")
  expect_identical(resolve_dff_refit_controls(mml_starts_fixture()$fit)$mml_starts, c(0,.5,-.5))
})

test_that("start review uses actual engines and preserves raw EM review", {
  # Full-suite scenario: ordinary checks retain live starts and engine-specific tests.
  skip_on_cran()
  d <- mml_starts_fixture()$data
  em <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    mml_engine = "em", mml_starts = c(0, .75, -.75), maxit = 1000, reltol = 1e-10))
  expect_true(all(em$opt$start_review$table$EngineUsed == "em"))
  expect_true(all(em$opt$start_review$table$NumericalStatus != "pass"))
  expect_false(em$summary$InferenceReady)
  joint <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    facet_interactions = "Rater:Criterion", mml_engine = "em", mml_starts = 0))
  expect_identical(joint$opt$start_review$table$EngineRequested, "em")
  expect_identical(joint$opt$start_review$table$EngineUsed, "em")
  expect_false(joint$opt$mml_engine$Fallback)
})

test_that("GPCM adaptive review preserves the integration basis and slope constraint", {
  # Full-suite scenario: ordinary checks retain live starts and engine-specific tests.
  skip_on_cran()
  d <- mml_starts_fixture()$data
  f <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
    mml_integration = "adaptive", mml_starts = c(0, .5, -.5), reltol = 1e-10))
  r <- f$opt$start_review; sizes <- build_param_sizes(f$config)
  idx <- build_indices(f$prep, f$config$step_facet, f$config$slope_facet, f$config$interaction_specs)
  evaluate <- mfrmr_make_adaptive_mml_evaluator(idx, f$config, sizes, 31L)
  expect_identical(r$integration, list(mode = "adaptive", points = 31))
  expect_true("log_sigma2" %in% r$parameter_map$Block)
  for (i in seq_along(r$attempts)) {
    opt <- r$attempts[[i]]; expected <- evaluate(opt$par)
    expect_equal(r$table$NegativeLogLik[i], expected$value, tolerance = 1e-10)
    expect_equal(r$table$TerminalGradientSupNorm[i], max(abs(expected$gradient)), tolerance = 1e-10)
    p <- expand_params(opt$par, sizes, f$config)
    expect_length(p$slopes, length(unique(d$Criterion)))
    expect_equal(prod(p$slopes), 1, tolerance = 1e-10)
  }
  expect_identical(mfrmr_gqs_refit_arguments(f, d, 41)$mml_starts, c(0, .5, -.5))
  expect_warning(bad <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    model = "GPCM", step_facet = "Criterion", slope_facet = "Criterion",
    mml_integration = "adaptive",
    mml_starts = list(overflow = c("log_sigma2:residual_variance" = 10000)), reltol = 1e-10),
    "non-passing attempts")
  expect_identical(bad$opt$start_review$table$NumericalStatus, c("pass", "error"))
  expect_identical(bad$opt$start_review$table$Selected, c(TRUE, FALSE))
  expect_true(nzchar(bad$opt$start_review$table$Error[2]))
})

test_that("failed and better unfinished attempts remain visible and cannot lose to worse passes", {
  count <- 0L
  local_mocked_bindings(
    mfrmr_mml_optimizer_parameter_map = function(...) data.frame(Coordinate = "location", OptimizerIndex = 1L),
    make_param_cache = function(...) NULL,
    make_mfrm_direct_evaluator = function(...) list(value = function(p) p, gradient = function(p) if (p == 3) 0 else 1),
    run_mfrm_optimization = function(...) {
      count <<- count + 1L
      if (count == 2L) { warning("record this warning"); stop("failed trial") }
      list(par = if (count == 1L) 3 else 1, value = if (count == 1L) 3 else 1,
        convergence = if (count == 1L) 0L else 1L,
        mml_engine = list(Used = "direct"),
        optimizer_diagnostics = list(ConvergenceSeverity = if (count == 1L) "pass" else "review",
          GradientReviewTolerance = 1e-4))
    }, .package = "mfrmr")
  config <- list(estimation_control = list(mml_engine_requested = "direct", mml_integration = "fixed"))
  expect_warning(out <- mfrmr_run_mml_starts(0, c(0, 1, -1), NULL, NULL, config,
    NULL, 31, 20, 1e-9, "auto"), "non-passing attempts")
  t <- out$start_review$table
  expect_equal(out$par, 1)
  expect_identical(t$Selected, c(FALSE, FALSE, TRUE))
  expect_identical(t$NumericalStatus, c("pass", "error", "review"))
  expect_match(t$Warning[2], "record this warning")
  expect_match(t$Error[2], "failed trial")
  expect_true(is.na(t$NegativeLogLik[2]))
  expect_length(out$start_review$attempts, 3)
  expect_null(out$start_review$attempts[[2]])
  expect_equal(t$LogLikLoss[c(1,3)], c(2,0))
  expect_equal(t$MaxAbsCoordinateChange[c(1,3)], c(2,0))
  local_mocked_bindings(run_mfrm_optimization = function(...) stop("all failed"), .package = "mfrmr")
  err <- tryCatch(mfrmr_run_mml_starts(0, c(0,1), NULL, NULL, config, NULL,
    31, 20, 1e-9, "auto"), error = identity)
  expect_s3_class(err, "mfrm_start_error")
  expect_equal(nrow(err$start_review$table), 2)
  expect_true(all(err$start_review$table$NumericalStatus == "error"))
  expect_false(any(err$start_review$table$Selected))
  local_mocked_bindings(run_mfrm_optimization = function(...) list(par=1,value=1,
    convergence=0L,mml_engine=list(Used="direct"),optimizer_diagnostics=list(
      ConvergenceSeverity="pass",GradientReviewTolerance=1e-4,TerminalGradientSupNorm=0)),
    .package="mfrmr")
  expect_warning(checked <- mfrmr_run_mml_starts(0, 0, NULL, NULL, config, NULL,
    31,20,1e-9,"auto"), "non-passing attempts")
  expect_identical(checked$optimizer_diagnostics$ConvergenceSeverity,"review")
  expect_equal(checked$optimizer_diagnostics$TerminalGradientSupNorm,1)
  expect_equal(checked$optimizer_diagnostics$TerminalGradientRMS,1)
})

test_that("start review rejects invalid controls and unsupported estimators", {
  map <- data.frame(Coordinate = c("a", "b"))
  for (bad in list(numeric(), NA_real_, Inf, 1i, matrix(1,1), "random", list(a=Inf),
                   list(default=c(a=1)), list(x=c(c=1)), list(x=c(a=1,a=2))))
    expect_error(mfrmr_mml_start_vectors(c(0,0), bad, map), "mml_starts|named start")
  d <- mml_starts_fixture()$data
  args <- list(data=d, person="Person", facets=c("Rater","Criterion"), score="Score", mml_starts=0)
  for (extra in list(list(method="JML"), list(method="CML"), list(mml_integration="auto"),
                    list(checkpoint=list(path=tempfile())),
                    list(model="GPCM", slope_facet=c("Rater","Criterion"))))
    expect_error(do.call(fit_mfrm,c(args,extra)), "mml_starts")
  # NULL is a no-op even for strict alternative-estimator argument dispatch.
  local_mocked_bindings(mfrmr_fit_cml = function(arguments, supplied) {
    expect_false("mml_starts" %in% supplied)
    "default dispatch"
  }, .package="mfrmr")
  args$mml_starts <- NULL
  expect_identical(do.call(fit_mfrm,c(args,list(method="CML",mml_starts=NULL))), "default dispatch")
})

test_that("starting-value outcomes reach reports, long tables and executable replay", {
  fx <- mml_starts_fixture(); f <- fx$fit
  res <- mfrm_results(f, include="fit", compute="never")
  report <- mfrm_report(res)
  expect_identical(report$tables$starting_values, f$opt$start_review$table)
  expect_identical(res$tables$starting_values, f$opt$start_review$table)
  coords <- report$tables$starting_value_coordinates
  expect_equal(nrow(coords), 3*length(f$opt$par))
  expect_identical(coords$Initial[coords$Start == "default"], unname(f$opt$start_review$starts$default))
  expect_match(report$markdown, "Starting-value review")
  expect_match(report$markdown, "global maximum")
  path <- tempfile(fileext=".rds"); on.exit(unlink(path), add=TRUE)
  saveRDS(f,path); restored <- readRDS(path)
  expect_identical(summary(restored)$starting_values, f$opt$start_review$table)
  expect_no_error(plot(restored, draw=FALSE))
  expect_equal(coef(restored), coef(f))
  lines <- build_replay_fit_mfrm_lines(f$config$replay_inputs,f$population,NULL,
    f$config$source_columns,f$config)
  env <- new.env(parent=environment()); env$data <- fx$data
  env$anchors <- env$group_anchors <- NULL
  eval(parse(text=lines),env)
  expect_identical(env$fit$config$replay_inputs$mml_starts, c(0,.5,-.5))
  expect_equal(env$fit$opt$par,f$opt$par,tolerance=1e-10)
  expect_identical(env$fit$opt$start_review$starts,f$opt$start_review$starts)
})
