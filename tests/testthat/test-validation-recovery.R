# Repository-only recovery adapters. Retained fits and deterministic accounting
# fixtures are not additional simulated datasets or coverage evidence.
recovery <- new.env(parent = globalenv())
for (name in c("mml-stages-20261001", "output-summary-20261001", "recovery-20261002")) {
  path <- test_path("..", "..", "inst", "validation", paste0("mfrm-wide-map-", name, ".R"))
  skip_if_not(file.exists(path), "Repository-only recovery helper")
  sys.source(path, recovery)
}
recovery_ab <- function(truth = "S-PCM") {
  root <- test_path("..", "..", "validation-results", "mfrm-wide-map-r50-20261001", "ab-inputs")
  skip_if_not(dir.exists(root), "Frozen A/B inputs")
  input <- recovery$wide_recovery_ab_input(root, paste0("AB:paired-L1:N20:SD1:", truth), 1L)
  spec <- list(ConditionId = input$ConditionId, Replicate = 1L, Arm = "PCM-MML",
    args = list(data = input$data, person = "Person", score = "Score", facets = c("Rater", "Criterion"),
      model = "PCM", method = "MML", step_facet = "Criterion", rating_min = 0,
      rating_max = 3, category_policy = "preserve", population_formula = ~1))
  list(input = input, spec = spec)
}
recovery_pilot <- function(method = "JML") {
  root <- test_path("..", "..", "validation-results", "mfrm-facet-structure-pilot-20261001")
  skip_if_not(dir.exists(root), "Retained matched method witnesses")
  x <- readRDS(file.path(root, "inputs", "20-base_two_facet.rds"))
  for (f in c("Rater", "Criterion")) names(x$truth[[f]]) <- as.character(seq_along(x$truth[[f]]))
  rownames(x$truth$steps) <- names(x$truth$log_slopes) <- names(x$truth$Criterion)
  input <- recovery$wide_recovery_input("retained-K3-GPCM", 1L, x$data, x$truth, "retained test witness")
  fit <- readRDS(file.path(root, "fits", paste0("20-base_two_facet-", method, ".rds")))$fit
  args <- list(data = x$data, person = "Person", score = "Score", facets = x$facets,
    model = "GPCM", method = if (method == "MML") "MML" else "JML", step_facet = "Criterion",
    slope_facet = "Criterion", rating_min = 0, rating_max = 2, category_policy = "preserve")
  if (method == "MML") args$population_formula <- ~1
  if (method %in% c("JML2", "JML4")) {
    args$jml_correction_order <- as.integer(sub("JML", "", method)); args$jml_correction_sampling <- "fixed_rosters"
  }
  list(input = input, spec = list(ConditionId = input$ConditionId, Replicate = 1L, Arm = method, args = args),
    fitting = list(value = fit, error = ""))
}
recovery_collect <- function(x) recovery$wide_recovery_collect(x$input, x$spec, x$fitting, "initial:fixture")

test_that("frozen PCM targets preserve true contrasts, fixed slopes and all probability contexts", {
  x <- recovery_ab(); p <- recovery$wide_recovery_plan(x$input, x$spec)
  one <- p[p$Analysis == "initial_finite_return", ]
  expect_equal(one$ReferenceValue[one$Kind == "contrast"], c(-.3, -.5))
  expect_equal(sum(one$Eligibility == "fixed"), 3L)
  expect_true(all(one$ReferenceValue[one$Kind == "log_slope"] == 0))
  expect_equal(sum(one$Kind == "probability"), 180L)
  expect_equal(recovery$wide_recovery_contrasts(x$input)$Rater, matrix(c(1,-1,0), 1,
    dimnames = list("1-2", c("1","2","3"))))
  # Check the literal formula against a separate adjacent-category recurrence.
  for (truth in c("S-PCM", "S-GPCM")) {
    a <- recovery_ab(truth)$input; v <- recovery$wide_recovery_probabilities(a, a$truth)
    g <- a$grid; expected <- t(vapply(seq_len(nrow(g)), function(i) {
      eta <- g$Theta[i] - a$truth$Rater[g$Rater[i]] - a$truth$Criterion[g$Criterion[i]]
      odds <- exp(exp(a$truth$log_slopes[g$Criterion[i]]) * (eta-a$truth$steps[g$Criterion[i], ]))
      w <- c(1, cumprod(odds)); w/sum(w)
    }, numeric(4)))
    expect_equal(v, expected, ignore_attr = TRUE, tolerance = 1e-13)
    expect_equal(rowSums(v), rep(1, nrow(v)), ignore_attr = TRUE, tolerance = 1e-13)
  }
  bad <- x$input; bad$truth$Rater[1] <- 5
  expect_error(recovery$wide_recovery_plan(bad, x$spec))
  z <- recovery$wide_recovery_collect(x$input, x$spec, list(value = NULL, error = "fit failure"), "initial:error")
  expect_equal(nrow(z$results), 2L * (nrow(one)-3L))
  expect_true(all(z$results$Status == "error" & !z$results$Available))
  expect_true(all(recovery$wide_recovery_summary(z$plan, z$results)$Complete))
})

test_that("ordinary initial returns survive scoring refusal and retain numerical status separately", {
  x <- recovery_pilot()
  local_mocked_bindings(fit_mfrm = function(...) stop("No fitting"),
    prediction_source_scoring_readiness = function(...) stop("No source computation"),
    .package = "mfrmr")
  z <- recovery_collect(x)
  expect_true(all(z$results$Available)); expect_true(z$numerical_pass)
  expect_identical(z$diagnostics$readiness, x$fitting$value$readiness)
  x$fitting$value$summary$ConvergenceSeverity <- "review"
  bad <- recovery_collect(x)
  expect_true(all(bad$results$Available[bad$results$Analysis == "initial_finite_return"]))
  expect_false(any(bad$results$Available[bad$results$Analysis == "initial_numerical_pass"]))
  expect_identical(z$results$Estimate, bad$results$Estimate)
  expect_true(all(is.na(z$results$SE) & is.na(z$results$Lower)))
  x$fitting$value$prep$data$Score[1] <- (x$fitting$value$prep$data$Score[1]+1L) %% 3L
  expect_error(recovery_collect(x))
})

test_that("matched saved methods use labels and common native probability coordinates", {
  jobs <- lapply(c("MML","JML","JML2","JML4"), recovery_pilot)
  out <- lapply(jobs, recovery_collect)
  for (i in seq_along(jobs)) {
    x <- jobs[[i]]; fit <- x$fitting$value
    p <- recovery$wide_recovery_parameters(fit, x$input); g <- x$input$grid
    eta <- g$Theta - p$Rater[g$Rater] - p$Criterion[g$Criterion]
    native <- mfrmr:::category_prob_gpcm(eta, t(apply(p$steps, 1, function(z) c(0,cumsum(z)))),
      match(g$Criterion, names(p$Criterion)), exp(p$log_slopes))
    expect_equal(recovery$wide_recovery_probabilities(x$input, p), native, ignore_attr = TRUE, tolerance = 1e-12)
    x$fitting$value$facets$others <- fit$facets$others[rev(seq_len(nrow(fit$facets$others))), ]
    x$fitting$value$steps <- fit$steps[rev(seq_len(nrow(fit$steps))), ]
    x$fitting$value$slopes <- fit$slopes[rev(seq_len(nrow(fit$slopes))), ]
    expect_identical(recovery_collect(x)$results$Estimate, out[[i]]$results$Estimate)
  }
  p <- do.call(rbind, lapply(out, `[[`, "plan")); r <- do.call(rbind, lapply(out, `[[`, "results"))
  pair <- recovery$wide_recovery_paired(p, r, "MML", "JML2")
  expect_true(all(pair$grid$PlannedPairs == 1L & pair$grid$BothAvailable == 1L))
  expect_true(all(is.na(pair$grid$PairedMCSE)))
  p$ScaleReference[p$Arm == "JML2"] <- "wrong scale"
  expect_error(recovery$wide_recovery_paired(p, r, "MML", "JML2"), "scale references")
  x <- jobs[[3]]; x$fitting$value$jml_adjustment$point$status <- "unresolved_starts"
  z <- recovery_collect(x)
  expect_true(all(z$results$Available[z$results$Analysis == "initial_finite_return"]))
  expect_false(any(z$results$Available[z$results$Analysis == "initial_numerical_pass"]))
})

test_that("grid MCSE uses dataset losses and incomplete grids remain unavailable", {
  x <- recovery_ab(); x$spec$args$model <- "GPCM"; x$spec$args$slope_facet <- "Criterion"
  p0 <- recovery$wide_recovery_plan(x$input, x$spec)
  # Artificial accounting values, not fitted probability vectors.
  p <- do.call(rbind, lapply(1:3, function(i) transform(p0, Replicate = i, InputId = paste0("fixture-",i))))
  r <- p[c("ConditionId","Arm","Replicate","Target","Output","Analysis")]
  r$Status <- "returned"; r$Available <- TRUE; r$Estimate <- p$ReferenceValue+.01*r$Replicate
  r$SE <- r$Lower <- r$Upper <- NA_real_; r$SourceStage <- "initial:fixture"; r$Reason <- ""
  loss <- recovery$wide_recovery_grid_losses(p,r); s <- recovery$wide_recovery_grid_summary(loss)
  expect_equal(s$MSE, rep(mean(c(.01,.02,.03)^2),2))
  expect_equal(s$MSEMCSE, rep(sd(c(.01,.02,.03)^2)/sqrt(3),2))
  expect_equal(s$Planned, c(3L,3L))
  rs <- recovery$wide_recovery_summary(p,r)
  expect_lt(max(abs(rs$RMSEMCSE - sd(c(.01,.02,.03)^2)/(2*sqrt(mean(c(.01,.02,.03)^2))*sqrt(3)))), 1e-12)
  j <- which(p$Kind == "probability" & p$Replicate == 3L)[1]
  r$Available[j] <- FALSE; r$Reason[j] <- "Missing target"
  loss <- recovery$wide_recovery_grid_losses(p,r)
  expect_false(loss$Available[loss$Analysis == "initial_finite_return" & loss$Replicate == 3L])
  incomplete <- recovery$wide_recovery_grid_summary(recovery$wide_recovery_grid_losses(p,r[-j,]))
  expect_false(incomplete$Complete[incomplete$Analysis == "initial_finite_return"])
  expect_true(is.na(incomplete$MSE[incomplete$Analysis == "initial_finite_return"]))
  expect_error(recovery$wide_recovery_records(p,rbind(r,r[1,])), "Duplicate selected output")
})

test_that("initial phase reader keeps q31 even when selected source points use q121", {
  root <- test_path("..","..","validation-results","mfrm-wide-map-r50-20261001","ab-inputs")
  dir <- test_path("..","..","validation-results","mfrm-wide-map-multi-output-20261001","native60")
  skip_if_not(file.exists(file.path(dir,"selected.rds")))
  input <- recovery$wide_recovery_ab_input(root,"AB:full-L1:N60:SD1:S-GPCM",1L)
  z <- recovery$wide_recovery_initial_phase(input,dir)
  expect_true(all(grepl("^initial:q31:",z$results$SourceStage)))
  expect_match(readRDS(file.path(dir,"selected.rds"))$payload$point$stage,"^q121:")
})
