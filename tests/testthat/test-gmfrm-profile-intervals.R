gmfrm_profile_fixture <- function() readRDS(test_path("fixtures", "gmfrm-profile-intervals.rds"))

test_that("aggregated profiles preserve total marginal scores and invalid-trial recovery", {
  numeric_gradient <- function(f, x) vapply(seq_along(x), function(j) {
    lo <- hi <- x; lo[j] <- lo[j] - 1e-5; hi[j] <- hi[j] + 1e-5
    (f(hi) - f(lo)) / 2e-5
  }, numeric(1))
  # Cover binary/polytomous responses, absent crossings, arbitrary names,
  # nonstandard factor order, and both declared owner orders.
  for (k in c(1L, 3L)) for (reverse in c(FALSE, TRUE)) {
    d <- expand.grid(Candidate = paste0("p", 1:5), `観点 名` = c("a", "b", "c"),
      Assessor = c("a", "b", "c"), check.names = FALSE)
    d$Rating <- rep(0:k, length.out = nrow(d))
    d <- d[!(d[["観点 名"]] == "c" & d$Assessor == "b"), ]
    d[["観点 名"]] <- factor(d[["観点 名"]], levels = c("c", "a", "b"))
    owners <- c("観点 名", "Assessor"); if (reverse) owners <- rev(owners)
    q <- gauss_hermite_normal(31L)
    p <- mfrm_gmfrm_problem(d[nrow(d):1, ], k, q, owners, "Candidate", "Rating")
    s <- p$common
    direct <- make_mfrm_direct_evaluator("MML",
      make_param_cache(s$sizes, s$config, s$idx, is_mml = TRUE), s$idx, s$config, s$sizes, q)
    grouped <- mfrm_gmfrm_profile_evaluator(p, direct)
    x <- unname(p$start + seq(-.3, .4, length.out = length(p$start)))
    for (v in list(x, x + .01, x)) {
      expect_lte(abs(grouped$value(v) - direct$value(v)), 1e-10)
      expect_lte(max(abs(grouped$gradient(v) - direct$gradient(v))), 1e-10)
      expect_equal(grouped$gradient(v), numeric_gradient(grouped$value, v), tolerance = 1e-6)
    }
    # Rejected optim() trials must not become stale-cache hits or valid fits.
    safe <- make_mfrm_boundary_safe_objective(grouped)
    bad <- x; bad[length(bad)] <- 1000
    for (i in 1:2) {
      expect_error(grouped$value(bad), class = "mfrmr_gpcm_slope_numeric_boundary_error")
      expect_identical(safe$value(bad), safe$penalty)
      expect_equal(safe$gradient(bad), numeric(length(bad)))
      expect_lte(abs(grouped$value(x) - direct$value(x)), 1e-10)
    }
    # Exercise the optimizer's actual mutable callback path, not just R copies.
    fn <- function(v) grouped$value(v) + sum(v^2)
    gr <- function(v) grouped$gradient(v) + 2*v
    solved <- optim(x, fn, gr, method = "BFGS", control = list(maxit = 8))
    expect_lte(abs(fn(solved$par) - direct$value(solved$par) - sum(solved$par^2)), 1e-10)
    expect_lte(max(abs(gr(solved$par) - direct$gradient(solved$par) - 2*solved$par)), 1e-10)
  }
})

# Literal response formula and normal integration, independent of the native
# probability/likelihood evaluator. Coordinates come from the documented fixture.
gmfrm_profile_literal_nll <- function(par, data, order = 241L) {
  rule <- gauss_hermite_normal(order)
  task <- as.integer(factor(data$Task)); rater <- as.integer(factor(data$Rater))
  nr <- max(rater) # Both fixtures have three tasks and scores 0:2.
  location <- c(par[1:2], -sum(par[1:2]))[task] + par[2+seq_len(nr)][rater]
  task_log_slope <- par[2+2*nr+1:2]
  slope <- exp(c(task_log_slope, -sum(task_log_slope))[task] + par[4+2*nr+seq_len(nr)][rater])
  step <- par[2+nr+seq_len(nr)][rater]
  response_logp <- vapply(rule$nodes, function(theta) {
    z <- cbind(0, slope * (theta - location - step), 2 * slope * (theta - location))
    shift <- apply(z, 1, max)
    z[cbind(seq_len(nrow(data)), data$Score + 1L)] - shift - log(rowSums(exp(z - shift)))
  }, numeric(nrow(data)))
  joint <- sweep(rowsum(response_logp, data$Person), 2, log(rule$weights), "+")
  shift <- apply(joint, 1, max)
  -sum(shift + log(rowSums(exp(joint - shift))))
}

test_that("EM M-step refinement meets the requested score while retaining ascent", {
  old <- readRDS(test_path("fixtures", "gmfrm-numerical-repair.rds"))$em
  expect_false(old$summary$Converged)
  expect_gt(tail(old$opt$em_trace$max_score,1),1e-7)
  p <- do.call(mfrm_gmfrm_problem,old$gmfrm$specification)
  result <- mfrm_gmfrm_em(p,maxit=500L,score_tol=1e-7)
  expect_true(result$converged)
  expect_lte(result$max_score,1e-7)
  expect_true(any(result$trace$mstep_polished))
  expect_gte(min(result$trace$q_gain,na.rm=TRUE),-1e-12)
  expect_gte(min(diff(result$trace$logLik)),-1e-10)
  expect_lte(-result$logLik,old$opt$value+1e-9)
  expect_lte(abs(-result$logLik-gmfrm_profile_literal_nll(result$par,p$specification$data)),1e-5)
  native <- mfrm_gmfrm_fit_result(p,result)
  expect_true(native$summary$Converged)
  expect_identical(native$opt$em_trace,result$trace)
  # Failed curvature refinement must preserve the ordinary EM step and reason.
  local_mocked_bindings(mfrm_optimizer_curvature_proposal=function(...) list(par=NULL,error="unusable curvature"),
    .package="mfrmr")
  rejected <- mfrm_gmfrm_em(p,start=old$opt$par,maxit=1L,score_tol=1e-7)
  expect_false(rejected$converged)
  expect_false(any(rejected$trace$mstep_polished))
  expect_identical(tail(rejected$trace$mstep_polish_error,1),"unusable curvature")
  expect_gte(min(diff(rejected$trace$logLik)),-1e-10)
})

test_that("two-family result assembly preserves normal rules across eigensolver roundoff", {
  fit <- gmfrm_profile_fixture()$fit
  rule <- gauss_hermite_normal
  alternative <- new.env(parent = environment(rule))
  alternative$eigen <- function(x, symmetric, only.values)
    base::eigen(x, symmetric = symmetric, only.values = FALSE)
  environment(rule) <- alternative
  spec <- fit$gmfrm$specification
  spec$quadrature <- rule(length(spec$quadrature$nodes))
  problem <- do.call(mfrm_gmfrm_problem, spec)
  result <- mfrm_gmfrm_em(problem, start = fit$opt$par, maxit = 1L, score_tol = 1e-7)
  expect_true(result$converged)
  native <- mfrm_gmfrm_fit_result(problem, result)
  expect_identical(native$config$estimation_control$quadrature, spec$quadrature)
  expect_identical(native$opt$par, result$par)
  # Substantive changes remain invalid, including changes to tiny tail weights.
  for (field in c("nodes", "weights")) {
    changed <- problem
    if (field == "nodes") changed$specification$quadrature$nodes[1] <-
      changed$specification$quadrature$nodes[1] + 1e-4 else
        changed$specification$quadrature$weights[1] <-
          changed$specification$quadrature$weights[1] * 1.001
    changed_result <- result
    changed_result$specification <- changed$specification
    expect_error(mfrm_gmfrm_fit_result(changed, changed_result),
                 "requires standard-normal Gauss-Hermite")
  }
})

test_that("neutral two-family profile starts do not escape under total-likelihood scaling", {
  x <- readRDS(test_path("fixtures", "gmfrm-numerical-repair.rds"))
  fit <- x$profile; spec <- fit$gmfrm$specification
  context <- function(q) {
    spec$quadrature <- gauss_hermite_normal(q)
    p <- do.call(mfrm_gmfrm_problem,spec); s <- p$common
    direct <- make_mfrm_direct_evaluator("MML",make_param_cache(s$sizes,s$config,s$idx,is_mml=TRUE),
      s$idx,s$config,s$sizes,spec$quadrature)
    mfrm_gmfrm_profile_evaluator(p,direct)
  }
  low <- context(121L); high <- context(241L)
  expect_equal(low$optimization_scale,240)
  a <- numeric(length(fit$opt$par)); a[15:16] <- -1
  ans <- lapply(list(fit$opt$par,rep(0,length(a))),function(start)
    mfrm_gpcm_profile_nuisance(low,high,start,a,x$log_slope,400L))
  expect_gt(x$old_checks$neutral$NLL-x$old_checks$retained$NLL,80)
  for(z in ans) {
    expect_true(z$check$Passed)
    expect_lte(abs(z$check$NLL-x$old_checks$retained$NLL),1e-5)
    expect_equal(sum(a*z$parameters),x$log_slope,tolerance=1e-12)
    expect_lte(abs(z$check$NLL-gmfrm_profile_literal_nll(z$parameters,spec$data)),1e-5)
  }
  expect_lte(abs(ans[[1]]$check$NLL-ans[[2]]$check$NLL),1e-5)
})

test_that("curvature coordinates recover a weak nuisance slope without changing the target", {
  x <- readRDS(test_path("fixtures","gmfrm-numerical-repair.rds"))$weak_profile
  fit <- x$fit; spec <- fit$gmfrm$specification
  context <- function(q) {
    spec$quadrature <- gauss_hermite_normal(q)
    p <- do.call(mfrm_gmfrm_problem,spec); s <- p$common
    direct <- make_mfrm_direct_evaluator("MML",make_param_cache(s$sizes,s$config,s$idx,is_mml=TRUE),
      s$idx,s$config,s$sizes,spec$quadrature)
    mfrm_gmfrm_profile_evaluator(p,direct)
  }
  low <- context(121L); high <- context(241L)
  a <- numeric(length(fit$opt$par)); a[15:16] <- -1
  ans <- lapply(list(fit$opt$par,rep(0,length(a))),function(start)
    mfrm_gpcm_profile_nuisance(low,high,start,a,x$log_slope,400L))
  expect_true(all(vapply(x$attempts,function(z) z$check$OptimizerCode==1L,TRUE)))
  for(z in ans) {
    expect_true(z$check$Passed)
    expect_identical(z$check$OptimizerCode,0L)
    expect_true(any(vapply(z$stages,function(s) identical(s$stage,"curvature_scaled"),TRUE)))
    expect_true(any(vapply(z$stages,function(s) identical(s$code,1L),TRUE)))
    expect_equal(sum(a*z$parameters),x$log_slope,tolerance=1e-12)
    expect_lte(abs(z$check$NLL-gmfrm_profile_literal_nll(z$parameters,spec$data)),1e-5)
  }
  expect_lte(abs(ans[[1]]$check$NLL-ans[[2]]$check$NLL),1e-5)
  # The optimizer's coordinate change cannot remove a genuine integration error.
  shifted <- high; shifted$value <- function(par) high$value(par)+.01
  bad <- mfrm_gpcm_profile_nuisance(low,shifted,fit$opt$par,a,x$log_slope,400L)
  expect_false(bad$check$Passed); expect_gt(bad$check$NLLChange,.009)
  local_mocked_bindings(mfrm_optimizer_curvature_scale=function(...) list(transform=NULL,error="unusable scaling"),
    .package="mfrmr")
  failed <- mfrm_gpcm_profile_nuisance(low,high,fit$opt$par,a,x$log_slope,1L)
  expect_false(failed$check$Passed)
  expect_identical(tail(failed$stages,1)[[1]]$error,"unusable scaling")
})

test_that("both families profile the literal marginal likelihood with nuisance reoptimization", {
  x <- gmfrm_profile_fixture(); fit <- x$fit; data <- fit$gmfrm$specification$data
  before <- serialize(fit, NULL)
  baseline <- gmfrm_profile_literal_nll(fit$opt$par, data)
  for (ci in x$profiles) {
    d <- attr(ci, "diagnostics"); pr <- attr(ci, "profile")
    expect_true(d$CIEligible)
    expect_true(all(pr$endpoints$Status == "computed"))
    expect_true(all(pr$profile$Passed))
    expect_lte(abs(pr$baseline - baseline), 1e-5)
    expect_true(d$SlopeOwner %in% c("Task", "Rater"))
    contrast <- numeric(13)
    if (d$SlopeOwner == "Task") contrast[9:10] <- -1 else contrast[13] <- 1
    expect_equal(sum(contrast * fit$opt$par), log(d$Estimate), tolerance = 1e-12)
    for (bound in pr$endpoints$Bound) {
      at <- which.min(abs(pr$profile$LogSlope - log(bound)))
      key <- sprintf("%.17g", pr$profile$LogSlope[at])
      for (attempt in pr$attempts[[key]]) {
        p <- attempt$parameters
        expect_true(attempt$check$Passed)
        expect_equal(sum(contrast * p), log(bound), tolerance = 1e-6)
        # Changes must include nuisance locations/steps, not only the selected slope.
        expect_gt(max(abs(p[1:8] - fit$opt$par[1:8])), 1e-4)
        value <- gmfrm_profile_literal_nll(p, data)
        expect_lte(abs(value - attempt$check$NLL), 1e-5)
        expect_lte(abs(2 * (value - baseline) - qchisq(.95, 1)), 1e-4)
      }
    }
  }
  expect_identical(serialize(fit, NULL), before)
})

test_that("two-family profile targets need unambiguous owner and level identities", {
  x <- gmfrm_profile_fixture(); fit <- x$fit
  for (s in list(NULL, "t3", c(Task = "unknown"), c(unknown = "t3"), c(Task = "t1", Rater = "r1")))
    expect_error(confint(fit, parm = "slopes", method = "profile", slope = s), "identify one owner and level")
  expect_error(confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3"), scale = "relative"), "Two-family intervals currently require")
  expect_error(confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3"), simultaneous = "bonferroni"), "multiplicity")
  expect_error(confint(fit, parm = "slopes", method = "profile", slope = c(Task = "t3"), profile_control = list(maxit = NULL)), "positive finite")
  expect_error(confint(fit, parm = "slopes", method = "sandwich"), "Two-family intervals currently require")
  failed <- fit; failed$summary$Converged <- FALSE
  expect_error(confint(failed, parm = "slopes", method = "profile", slope = c(Task = "t3")), "per-Person")
  # Cached QR factors are numerical workspace, not part of model identity.
  uncached <- failed
  uncached$config$gpcm_spec$log_slope_design@factors <- list()
  expect_error(confint(uncached, parm = "slopes", method = "profile", slope = c(Task = "t3")), "per-Person")
  changed <- uncached
  changed$config$gpcm_spec$log_slope_design@x[1] <-
    changed$config$gpcm_spec$log_slope_design@x[1] + .01
  expect_error(confint(changed, parm = "slopes", method = "profile", slope = c(Task = "t3")),
               "product specification must match")
  # Reuse completed searches; this tests dispatch, checks and result identity,
  # not another costly calculation of the same endpoints.
  local_mocked_bindings(mfrm_gpcm_profile_search = function(evaluator, reference, start, contrast, ...) {
    ci <- if (contrast[13] == 1) x$profiles$Rater else x$profiles$Task
    attr(ci, "profile")
  }, .package = "mfrmr")
  for (owner in c("Task", "Rater")) {
    s <- setNames(if (owner == "Task") "t3" else "r3", owner)
    ci <- suppressWarnings(confint(fit, parm = "slopes", method = "profile", slope = s))
    expect_equal(as.vector(ci), as.vector(x$profiles[[owner]]))
    expect_identical(attr(ci, "settings")$slope, s)
    expect_identical(attr(ci, "diagnostics")$SlopeOwner, owner)
    expect_identical(attr(ci, "diagnostics")$SlopeLevel, unname(s))
    expect_identical(attr(ci, "source"), mfrm_gpcm_inference_source(fit))
    expect_true(all(attr(ci, "checks")$Passed))
  }
  # Names and shared level labels are not statistical roles or identity keys.
  spec <- fit$gmfrm$specification
  names(spec$data) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  spec$data[["観点 名"]] <- sub("t", "level:", spec$data[["観点 名"]], fixed = TRUE)
  spec$data[["Judge-ID"]] <- sub("r", "level:", spec$data[["Judge-ID"]], fixed = TRUE)
  spec$slope_facets <- c("観点 名", "Judge-ID"); spec$person <- "Candidate"; spec$score <- "Rating"
  p <- do.call(mfrm_gmfrm_problem, spec)
  renamed <- mfrm_gmfrm_fit_result(p, mfrm_gmfrm_em(p, start = fit$opt$par, maxit = 1, score_tol = 1e-7))
  for (owner in spec$slope_facets) {
    ci <- suppressWarnings(confint(renamed, parm = "slopes", method = "profile", slope = setNames("level:3", owner)))
    expect_identical(attr(ci, "diagnostics")$SlopeOwner, owner)
    expect_identical(attr(ci, "diagnostics")$SlopeLevel, "level:3")
  }
})

test_that("two-family saved profiles retain failures and meanings through existing outputs", {
  x <- gmfrm_profile_fixture()
  local_mocked_bindings(mfrm_gpcm_profile_search = function(...) stop("must not reprofile"),
    mfrm_gpcm_product_inference = function(...) stop("must not recompute covariance"),
    fit_mfrm = function(...) stop("must not refit"), .package = "mfrmr")
  for (ci in x$profiles) {
    p <- plot(ci, draw = FALSE)
    expect_identical(p$labels$x, "Component discrimination (log scale)")
    expect_match(paste(plot_data(p, "notes")$Text, collapse = "\n"), "Dotted: Wald limits")
    expect_identical(plot_data(p, "slope_metadata"), attr(ci, "diagnostics"))
    expect_identical(plot_data(as_ggplot(attr(p, "mfrmr_plot_data"))), plot_data(p))
    expect_no_error(ggplot2::ggplot_build(p))
  }
  limited <- x$unresolved
  p <- plot(limited, title = NULL, subtitle = NULL, caption = NULL, draw = FALSE)
  expect_null(p$labels$title); expect_null(p$labels$subtitle); expect_null(p$labels$caption)
  expect_true(is.na(plot_data(p, "endpoints")$Bound[2]))
  expect_false(attr(limited, "diagnostics")$CIEligible)
  expect_match(attr(limited, "diagnostics")$InferenceReview, "numerical checks")
  res <- mfrm_results(x$fit, include = c("fit", "plots"), compute = "never", intervals = x$profiles)
  report <- mfrm_report(res)
  expect_match(report$markdown, "Profile likelihood; chi-square")
  expect_identical(report$tables$gpcm_rater_profile_endpoints, attr(x$profiles$Rater, "profile")$endpoints)
  expect_identical(res$tables$gpcm_task_settings$SlopeOwner, "Task")
  folder <- tempfile(); on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  exported <- export_mfrm_results(res, output_dir = folder, preset = "starter", acknowledge_sensitive = TRUE)
  expect_equal(nrow(exported$plot_errors), 0L)
  saved <- readRDS(exported$written_files$Path[exported$written_files$Component == "results_rds"])
  expect_identical(saved$gpcm_inference, res$gpcm_inference)
  replay <- exported$written_files$Path[exported$written_files$Component == "replay_code"]
  withr::with_dir(folder, sys.source(replay, envir = new.env(parent = globalenv())))
})

test_that("a recovered real interval retains inaccurate outward trials in public output", {
  x <- readRDS(test_path("fixtures","gmfrm-numerical-repair.rds"))$recovered
  ci <- x$ci; pr <- attr(ci,"profile")
  expect_true(attr(ci,"diagnostics")$CIEligible)
  expect_true(all(pr$endpoints$Status == "computed"))
  expect_true(any(!pr$profile$Passed))
  expect_length(pr$attempts,nrow(pr$profile))
  local_mocked_bindings(mfrm_gpcm_profile_search=function(...) stop("must not reprofile"),
    mfrm_gpcm_product_inference=function(...) stop("must not recompute covariance"),
    fit_mfrm=function(...) stop("must not refit"),.package="mfrmr")
  p <- plot(ci,draw=FALSE)
  expect_identical(plot_data(p,"table")$Passed,pr$profile$Passed)
  expect_true(all(plot_data(p,"table")$Display[!pr$profile$Passed] == "Unresolved"))
  expect_identical(plot_data(p,"endpoints"),pr$endpoints)
  expect_no_error(ggplot2::ggplot_build(p))
  res <- mfrm_results(x$fit,include=c("fit","plots"),compute="never",intervals=list(recovered=ci))
  expect_identical(mfrm_report(res)$tables$gpcm_recovered_profile,pr$profile)
  path <- tempfile(fileext=".rds"); on.exit(unlink(path),add=TRUE)
  saveRDS(res,path)
  expect_identical(readRDS(path)$gpcm_inference,res$gpcm_inference)
})
