gmfrm_public_args <- function() {
  d <- expand.grid(Person = paste0("p", 1:12), Task = c("same:1", "t2", "t3"),
    Rater = c("same:1", "r2"))
  d <- d[!(d$Task == "t3" & d$Rater == "r2"), ]
  d$Score <- rep(c(0, 2, 1, 1, 2, 0, 2, 1, 0, 1), length.out = nrow(d))
  names(d) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  list(data = d, person = "Candidate", facets = c("Judge-ID", "観点 名"),
    score = "Rating", model = "GPCM", method = "MML",
    slope_facet = c("観点 名", "Judge-ID"), step_facet = "Judge-ID",
    noncenter_facet = "Judge-ID", gpcm_mml_identification = "fixed_standard_normal",
    mml_engine = "em", quad_points = 7L, maxit = 1L)
}

test_that("public two-family dispatch preserves its explicit likelihood and replay", {
  args <- gmfrm_public_args()
  fit <- do.call(fit_mfrm, args)
  expect_s3_class(fit, "mfrm_fit")
  expect_identical(fit$config$slope_facet, args$slope_facet)
  expect_identical(fit$summary$MMLEngineRequested, "em")
  expect_identical(fit$summary$MMLEngineUsed, "em")
  expect_identical(fit$summary$ConvergenceBasis, "marginal_score_per_person")
  expect_false(fit$summary$Converged)
  expect_false(fit$summary$InferenceReady)
  expect_false(fit$summary$ICEligible)
  problem <- do.call(mfrm_gmfrm_problem, fit$gmfrm$specification)
  reference <- mfrm_gmfrm_em(problem, maxit = 1L)
  expect_identical(fit$opt$par, reference$par)
  expect_equal(fit$summary$LogLik, reference$logLik)
  expect_identical(summary(fit)$slope_overview$Slopes, c(3L, 2L))
  replay <- fit$config$replay_inputs
  replay$package_version <- NULL
  replay$data <- args$data
  replayed <- do.call(fit_mfrm, replay)
  expect_identical(replayed$opt$par, fit$opt$par)
  expect_identical(replayed$prep$score_map, fit$prep$score_map)
  expect_identical(replayed$summary, fit$summary)
})

test_that("two-family input never silently changes requested model options", {
  args <- gmfrm_public_args()
  bad <- function(name, value) { a <- args; a[name] <- list(value); do.call(fit_mfrm, a) }
  expect_error(bad("mml_engine", "direct"), "Two slope families require")
  expect_error(bad("method", "JML"), "method = 'MML'", fixed = TRUE)
  expect_error(bad("gpcm_mml_identification", "free_population"), "fixed_standard_normal")
  expect_error(bad("step_facet", "観点 名"), "second slope facet")
  expect_error(bad("noncenter_facet", "Person"), "second slope facet")
  expect_error(bad("slope_facet", c("観点 名", "Judge-ID", "Third")), "exactly two")
  for (option in c("anchors", "weight", "population_formula", "reltol", "checkpoint",
                   "attach_diagnostics", "missing_codes", "facet_shrinkage")) {
    value <- switch(option, reltol=1e-7, attach_diagnostics=TRUE,
      facet_shrinkage="laplace", population_formula=~1, "unsupported")
    expect_error(bad(option, value), "not supported by the two-family route")
  }
  expect_error(bad("optimizer", "L-BFGS-B"), "optimizer")
  expect_error(bad("em_score_tol", 0), "finite positive")
  expect_error(bad("rating_min", 1), "rating_min = 0")
  expect_error(bad("rating_min", NA), "rating_min = 0")
  expect_error(bad("rating_min", "zero"), "rating_min = 0")
  expect_error(bad("rating_max", c(1,2)), "positive integer")
  incomplete <- args$data; incomplete$Rating[1] <- NA
  expect_error(bad("data", incomplete), "observed finite numeric scores")
  collapsed <- args$data; collapsed$Rating[collapsed$Rating == 1] <- 2
  expect_error(bad("data", collapsed), "does not collapse score categories")
  preserved <- args; preserved$data <- collapsed
  preserved$category_policy <- "preserve"; preserved$rating_max <- 2
  # Preserve the declared gap and reject its unsupported contrast before fitting.
  expect_error(do.call(fit_mfrm, preserved), class="mfrmr_category_readiness_error")
  duplicate <- rbind(args$data, args$data[1,])
  expect_error(bad("data", duplicate), "Repeated person-facet")
  a <- args; a$slope_facet <- "Judge-ID"; a$em_score_tol <- 1e-6
  expect_error(do.call(fit_mfrm, a), "only to the two-family")
})

test_that("public conditional curves match the product equation and retain unavailable intervals", {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  # The fixture keeps Japanese facet names; the default PDF font cannot encode them.
  grDevices::pdf(NULL, family = "Japan1")
  on.exit(grDevices::dev.off(), add = TRUE)
  args <- gmfrm_public_args(); fit <- do.call(fit_mfrm, args)
  grid <- data.frame(Theta = rep(c(-1,0,1), 2))
  grid[[args$slope_facet[1]]] <- rep(c("same:1", "t3"), each=3)
  grid[[args$slope_facet[2]]] <- rep(c("same:1", "r2"), each=3)
  curves <- mfrm_curve_intervals(fit, grid)
  info <- mfrm_curve_intervals(fit, grid, type="information")
  problem <- do.call(mfrm_gmfrm_problem, fit$gmfrm$specification)
  pars <- problem$unpack(fit$opt$par)
  expected <- t(vapply(seq_len(nrow(grid)), function(j) {
    first <- grid[[args$slope_facet[1]]][j]; second <- grid[[args$slope_facet[2]]][j]
    slope <- pars$slopes[[1]][first] * pars$slopes[[2]][second]
    eta <- grid$Theta[j] - pars$locations[[1]][first] - pars$locations[[2]][second]
    z <- c(0, cumsum(slope * (eta - pars$steps[second,])))
    exp(z-max(z))/sum(exp(z-max(z)))
  }, numeric(3)))
  expect_equal(curves$table$Estimate, as.vector(t(expected)), tolerance=1e-12)
  variance <- rowSums(expected * (matrix(0:2,nrow(grid),3,byrow=TRUE) - drop(expected %*% (0:2)))^2)
  slopes <- pars$slopes[[1]][grid[[args$slope_facet[1]]]] * pars$slopes[[2]][grid[[args$slope_facet[2]]]]
  expect_equal(info$table$Estimate, unname(slopes^2 * variance), tolerance=1e-12)
  expect_identical(curves$contexts$ObservedContext, rep(c(TRUE,FALSE),each=3))
  expect_true(all(!curves$table$CIEligible))
  expect_true(all(is.na(curves$table$SE) & is.na(curves$table$Lower) & is.na(curves$table$Upper)))
  expect_true(curves$settings$point_only)
  expect_error(mfrm_curve_intervals(fit,grid,method="sandwich"), "do not yet support")
  expect_error(mfrm_curve_intervals(fit,grid,simultaneous="bonferroni"), "not supported")
  p <- plot(curves, draw=FALSE)
  expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "intervals are unavailable")
  expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "unavailable for two-family GPCM")
  expect_false(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1))))
  contexts <- unique(plot_data(p, "table")$Context)
  labels <- p$facet$params$labeller(data.frame(Context = contexts))[[1L]]
  expect_true(all(grepl("\n", labels, fixed = TRUE)))
  expect_true(all(grepl("Judge-ID = ", labels, fixed = TRUE)))
  expect_true(all(grepl("観点 名 = ", labels, fixed = TRUE)))
  expect_false(grepl("95|Approximate", paste(plot_data(p, "notes")$Text, collapse = "\n")))
  expect_null(plot(curves, draw=FALSE, title=NULL,subtitle=NULL,caption=NULL)$labels$title)
  expect_equal(plot_data(p, component="table")$Estimate, curves$table$Estimate)
  single <- mfrm_curve_intervals(fit, grid[1, , drop = FALSE])
  one_point <- plot(single, draw = FALSE)
  built <- ggplot2::ggplot_build(one_point)
  expect_equal(nrow(built$data[[2L]]), 0L)
  expect_equal(built$data[[3L]]$y, single$table$Estimate)
  expect_true(all(built$data[[3L]]$shape == 19))
  expect_no_warning(ggplot2::ggplotGrob(one_point))
  ci <- confint(fit, parm = "slopes")
  local_mocked_bindings(fit_mfrm=function(...) stop("unexpected fitting"),
    diagnose_mfrm=function(...) stop("unexpected diagnostics"),
    mfrm_curve_intervals=function(...) stop("unexpected curves"),
    compute_mml_parameter_covariance=function(...) stop("unexpected information"),
    .package="mfrmr")
  res <- mfrm_results(fit, include=c("fit","plots"), compute="never", intervals=list(curves=curves))
  expect_true(res$plot_map$Available[res$plot_map$Type == "gpcm_curves"])
  expect_false(any(res$plot_map$Available[res$plot_map$Type %in% c("wright","fit_pathway")]))
  expect_true(all(res$tables$gpcm_curves_curves$ConfidenceLevel == "Not available"))
  report <- mfrm_report(res)
  expect_s3_class(report, "mfrm_report")
  expect_match(report$markdown, "Conditional new-Person EAP is separately available", fixed = TRUE)
  expect_match(report$markdown, "Calibration intervals are unavailable", fixed=TRUE)
  expect_false(grepl("## GPCM uncertainty", report$markdown, fixed=TRUE))
  expect_false(grepl("required_wright_map_missing|Rebuild.*diagnostics|diagnose_mfrm\\(", report$markdown))
  expect_s3_class(summary(report), "summary.mfrm_report")
  expect_identical(report$tables$slope_estimates, fit$slopes)
  expect_equal(report$tables$fitted_slopes$FittedSlope, fit$slopes$OptimizerEstimate)
  for (name in c("fitted_slopes", "fitted_locations", "fitted_steps"))
    expect_identical(report$tables[[name]], res$tables[[paste0("gpcm_", name)]])
  expect_identical(report$tables$fitted_locations$FittedLocation, unname(fit$facets$others$Estimate))
  expect_identical(report$tables$fitted_locations[c("Facet", "Level")], as.data.frame(fit$facets$others[c("Facet", "Level")]))
  expect_identical(report$tables$fitted_steps$FittedStep, unname(fit$steps$Estimate))
  expect_identical(report$tables$fitted_steps$StepOwner, rep(fit$config$step_facet, nrow(fit$steps)))
  expect_identical(report$tables$fitted_steps$StepLevel, fit$steps$StepFacet)
  expect_identical(report$tables$fitted_steps$Step, fit$steps$Step)
  expect_true(all(report$tables$fitted_locations$Interval == "See attached location-interval tables, if requested"))
  expect_true(all(report$tables$fitted_steps$Interval == "Not available"))
  expect_identical(report$tables$fitted_locations$Constraint,
    rep(c("Sum to zero across levels", "Uncentered; fixed N(0,1) ability"), c(3,2)))
  expect_match(report$markdown, "## fitted locations", fixed=TRUE)
  expect_match(report$markdown, "## fitted steps", fixed=TRUE)
  expect_match(report$markdown, "not standalone category thresholds", fixed=TRUE)
  legacy <- res
  legacy$tables <- legacy$tables[!startsWith(names(legacy$tables), "gpcm_fitted_")]
  expect_identical(mfrm_report(legacy)$tables$fitted_locations, report$tables$fitted_locations)
  expect_identical(mfrm_report(legacy)$tables$fitted_steps, report$tables$fitted_steps)
  expect_identical(report$tables$model_settings$Value[4:5], fit$config$slope_facet)
  expect_s3_class(plot(res,type="gpcm_curves",draw=FALSE), "ggplot")
  expect_identical(plot_data(plot(res,draw=FALSE)), plot_data(plot(res,type="gpcm_curves",draw=FALSE)))
  expect_equal(plot_data(as_ggplot(res),component="table")$Estimate, curves$table$Estimate)
  slopes <- mfrm_results(fit, intervals=list(slopes=ci))
  expect_identical(plot_data(plot(slopes,draw=FALSE)), plot_data(plot(ci,draw=FALSE)))
  expect_false(any(plot_data(plot(slopes,draw=FALSE),component="table")$CIEligible))
  multiple <- mfrm_results(fit, intervals=list(curves=curves, slopes=ci))
  expect_error(plot(multiple,draw=FALSE), "choose `type` from: gpcm_curves, gpcm_slopes", fixed=TRUE)
  expect_error(plot(multiple,type=NULL,draw=FALSE), "Several saved two-family plots")
  expect_identical(plot_data(plot(multiple,type="gpcm_slopes",draw=FALSE)), plot_data(plot(ci,draw=FALSE)))
  no_plots <- mfrm_results(fit, include="fit", intervals=curves)
  expect_error(plot(no_plots,draw=FALSE), "did not include plots")
  expect_error(plot(no_plots,type="gpcm_inference",draw=FALSE), "not available")
  file <- tempfile(fileext=".rds"); withr::defer(unlink(file))
  saveRDS(res,file)
  expect_identical(readRDS(file)$tables, res$tables)
  expect_identical(plot_data(plot(readRDS(file),draw=FALSE)), plot_data(plot(res,draw=FALSE)))
  export_dir <- tempfile("gmfrm-export-"); withr::defer(unlink(export_dir,recursive=TRUE))
  exported <- export_mfrm_results(res, output_dir=export_dir, preset="starter",
    acknowledge_sensitive=TRUE, plot_width=1200, plot_height=900)
  expect_equal(nrow(exported$plot_errors), 0L)
  for (name in c("fitted_slopes", "fitted_locations", "fitted_steps")) {
    for (prefix in c("table_gpcm_", "report_")) {
      path <- exported$written_files$Path[exported$written_files$Component == paste0(prefix,name)]
      expect_length(path, 1L)
      expect_equal(read.csv(path, encoding="UTF-8", check.names=FALSE), report$tables[[name]], tolerance=1e-12)
    }
  }
  index <- readLines(file.path(export_dir,"index.html"),warn=FALSE)
  expect_true(any(grepl("Curve intervals",index,fixed=TRUE)))
  expect_false(any(grepl("required Wright|Required scale display",index)))
  expect_true(file.exists(file.path(export_dir,"mfrmr_results_plot_gpcm_curves.png")))
  expect_error(mfrm_results(fit, include="standard"), "Two-family results support only")
})

test_that("two-family default results collect and replay saved evidence without computation", {
  fit <- do.call(fit_mfrm, gmfrm_public_args())
  local_mocked_bindings(fit_mfrm=function(...) stop("unexpected fitting"),
    diagnose_mfrm=function(...) stop("unexpected diagnostics"),
    mfrm_curve_intervals=function(...) stop("unexpected curves"),
    compute_mml_parameter_covariance=function(...) stop("unexpected information"),
    .package="mfrmr")
  res <- mfrm_results(fit)
  expect_identical(as.character(res$include), c("fit", "plots"))
  expect_identical(res$fit, fit)
  expect_null(res$diagnostics)
  expect_null(res$gpcm_inference)
  expect_true(all(paste0("gpcm_fitted_", c("slopes", "locations", "steps")) %in% names(res$tables)))
  expect_false(any(paste0("gpcm_fitted_", c("slopes", "locations", "steps")) %in%
    names(mfrm_results(fit, include="plots")$tables)))
  expect_identical(res$status$Status[res$status$Section == "diagnostics"], "not_available")
  expect_false(any(res$plot_map$RequiredArtifact))
  expect_false(any(res$plot_map$Available[res$plot_map$Type != "tables"]))
  triage <- summary(res)$triage
  expect_identical(triage$Signal[triage$Area == "Diagnostics"], "diagnostics_unsupported")
  expect_false(any(grepl("diagnose_mfrm|Rebuild|compute =", triage$Route)))
  expect_false(any(grepl("required_wright_map", triage$Signal)))
  expect_error(plot(res,draw=FALSE), "No saved two-family plot")
  html <- mfrm_results(fit,output="html"); withr::defer(unlink(html$path))
  expect_false(any(grepl("required_wright_map_missing|diagnose_mfrm", html$html)))
  expect_identical(mfrm_results(fit, compute="never")$tables, res$tables)
  for (include in list(NULL, character()))
    expect_identical(mfrm_results(fit, include=include)$tables, res$tables)
  expect_identical(mfrm_results(fit, include="fit", compute="auto")$tables, res$tables)
  expect_identical(mfrm_results(fit, output="tables"), res$tables)
  expect_identical(mfrm_results(fit, output="summary"), summary(res))
  report <- mfrm_report(res)
  expect_match(report$markdown, "EM stopped without meeting")
  expect_identical(report$tables$slope_estimates, fit$slopes)
  for (include in c("diagnostics", "standard"))
    expect_error(mfrm_results(fit, include=include), "Two-family results support only")
  expect_error(mfrm_results(fit, scores=data.frame()), "must be saved")
  folder <- tempfile("gmfrm-default-replay-"); dir.create(folder)
  withr::defer(unlink(folder, recursive=TRUE)); withr::local_dir(folder)
  code <- paste(summary(res)$reproducible_code$Code, collapse="\n")
  expect_match(code, "saveRDS", fixed=TRUE)
  expect_match(code, "readRDS", fixed=TRUE)
  replay <- new.env(); replay$res <- res
  eval(parse(text=code), envir=replay)
  expect_identical(replay$res, res)
  expect_identical(replay$report$markdown, report$markdown)
  exported <- export_mfrm_results(res, output_dir="export", include=c("replay", "report"),
    acknowledge_sensitive=TRUE)
  expect_true("rds" %in% exported$written_files$Format)
  replay_file <- list.files("export", pattern="_replay[.]R$", full.names=FALSE)
  restored <- new.env(parent=globalenv())
  withr::with_dir("export", sys.source(replay_file, envir=restored))
  expect_identical(restored$res, res)
  expect_identical(restored$report$markdown, report$markdown)
})

test_that("unsupported two-family consumers give a specific scope error", {
  fit <- do.call(fit_mfrm, gmfrm_public_args())
  for (fun in list(diagnose_mfrm, plot, compute_information,
      plot_wright_unified, plot_threshold_ladder, bootstrap_mfrm_gpcm,
      extract_mfrm_sim_spec, make_anchor_table, apply_empirical_bayes_shrinkage,
      analyze_facet_equivalence, estimate_bias, estimate_all_bias,
      plot_rater_severity_profile, export_mfrm_tables)) {
    expect_error(fun(fit), class="mfrmr_gpcm_scope_error")
  }
  expect_error(mfrm_response_diagnostics(fit), "Resolve the two-family convergence")
  expect_error(predict_mfrm_units(fit, fit$gmfrm$specification$data), "Resolve the two-family convergence")
  expect_error(extract_mfrm_calibration(fit), "Resolve the two-family convergence")
  expect_error(compare_mfrm(fit,fit), class="mfrmr_gpcm_scope_error")
  expect_error(plot_compare_mfrm(fit,fit), class="mfrmr_gpcm_scope_error")
  expect_error(plot_rater_trajectory(list(a=fit,b=fit)), class="mfrmr_gpcm_scope_error")
  expect_error(summary(fit,profile="facets"), "only the fit summary")
  expect_error(summary(fit,profile="reporting"), "only the fit summary")
  ci <- confint(fit, parm = "slopes")
  expect_true(all(is.na(ci)))
  expect_match(attr(ci,"diagnostics")$InferenceReview[1],"per-Person")
})
