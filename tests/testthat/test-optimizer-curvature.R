test_that("polishing cannot trade resolved likelihood loss for a smaller gradient", {
  # Two local minima have different heights. Adding a constant must not make
  # the worse one preferable; both gradients can satisfy the stopping rule.
  fn <- function(x, offset) offset + (x^2 - 1)^2 - .003 * (x^3/3 - x)
  gr <- function(x) (x^2 - 1) * (4*x - .003)
  stage <- function(x, offset) list(opt = list(par = x, value = fn(x, offset)),
    diagnostics = list(ConvergenceSeverity = "pass", TerminalGradientSupNorm = abs(gr(x))))
  for (offset in c(0, 1e6)) {
    current <- stage(-1 + 1e-6, offset)
    worse <- stage(1, offset)
    expect_gt(worse$opt$value - current$opt$value, .0039)
    expect_false(mfrm_optimizer_stage_is_better(worse, current))
    current$diagnostics$ConvergenceSeverity <- "review"
    expect_false(mfrm_optimizer_stage_is_better(worse, current))
    expect_true(mfrm_optimizer_stage_is_better(stage(-1, offset), current))
  }
  current <- stage(-1 + 1e-6, 1e6)
  rounded <- stage(-1, 1e6)
  rounded$opt$value <- current$opt$value + .Machine$double.eps * 1e6
  expect_true(mfrm_optimizer_stage_is_better(rounded, current))
  rounded$opt$value <- Inf
  expect_false(mfrm_optimizer_stage_is_better(rounded, current))
})

test_that("interaction polishing removes only the stopping constant and restores raw likelihood", {
  fn <- function(p) 1e6 + p^2
  bad <- FALSE
  local_mocked_bindings(
    make_param_cache = function(...) list(),
    make_mfrm_direct_evaluator = function(...) list(value=fn, gradient=function(p) 2*p,
      diagnostics=function() list()),
    mfrm_optimizer_curvature_proposal = function(...) list(par=NULL, error="not needed"),
    optim = function(par, fn, ...) {
      if (abs(fn(par)) < 1) {
        expect_equal(fn(par), 0)
        expect_equal(fn(1)-fn(0), 1)
        par <- if (bad) 100 else 0
      }
      list(par=par,value=fn(par),convergence=0L,counts=c('function'=1L,'gradient'=0L))
    }, .package="mfrmr")
  run <- function() suppressWarnings(run_mfrm_direct_optimization(.01,"MML",list(),
    list(model="RSM",interaction_specs=list(cell=list(n_params=1L))), list(),
    quad_points=3,maxit=20,reltol=1e-10))
  out <- run()
  expect_identical(out$par, 0)
  expect_equal(out$value, fn(0))
  expect_identical(out$optimizer_diagnostics$ConvergenceSeverity, "pass")
  expect_identical(out$optimizer_polish$Stages$StageLabel[out$optimizer_polish$Stages$Selected],
    "centered_objective_polish")
  bad <- TRUE
  rejected <- run()
  expect_equal(rejected$par, .01)
  expect_equal(rejected$value, fn(.01))
  expect_identical(rejected$optimizer_diagnostics$ConvergenceSeverity, "review")
})

test_that("the optimizer retains its better point after rejected polishing stages", {
  fn <- function(x) 1e6 + (x^2 - 1)^2 - .003 * (x^3/3 - x)
  gr <- function(x) (x^2 - 1) * (4*x - .003)
  calls <- 0L
  local_mocked_bindings(
    make_param_cache = function(...) list(),
    make_mfrm_direct_evaluator = function(...) list(value = fn, gradient = gr, diagnostics = function() list()),
    optim = function(par, ...) {
      calls <<- calls + 1L
      x <- if (calls == 1L) par else 1
      list(par = x, value = fn(x), convergence = 0L, counts = c("function" = 1L, "gradient" = 1L))
    }, .package = "mfrmr")
  expect_warning(out <- run_mfrm_direct_optimization(-1 + 1e-4, "JML", list(),
    list(model = "RSM"), list(), quad_points = 3L, maxit = 20L, reltol = 1e-9), "convergence review")
  expect_gt(calls, 1L)
  expect_equal(out$par, -1 + 1e-4)
  expect_equal(out$value, fn(out$par))
  expect_identical(out$optimizer_diagnostics$ConvergenceSeverity, "review")
  expect_identical(out$optimizer_polish$Stages$StageLabel[out$optimizer_polish$Stages$Selected], "initial")
  expect_false(out$optimizer_polish$Succeeded)
})

test_that("curvature restart resolves a large-offset stationary point accurately", {
  h <- diag(c(900, 15000))
  fn <- function(p) 86000 + sum(p * (h %*% p)) / 2
  gr <- function(p) as.vector(h %*% p)
  start <- c(4e-7, -1e-7)
  proposal <- mfrm_optimizer_curvature_proposal(start, fn, gr)
  expect_identical(proposal$error, "")
  expect_lt(max(abs(gr(proposal$par))), 1e-10)
  expect_lte(fn(proposal$par), fn(start))
  expect_true(all(proposal$counts > 0))
  expect_equal(proposal$par, c(0, 0), tolerance = 1e-15)
})

test_that("curvature scaling changes search coordinates without changing the objective", {
  h <- matrix(c(5, 2, 2, -1), 2)
  fn <- function(p) drop(crossprod(p, h %*% p)) / 2
  gr <- function(p) drop(h %*% p)
  p <- c(.2, -.3)
  scale <- mfrm_optimizer_curvature_scale(p, fn, gr)
  expect_true(scale$negative)
  expect_identical(scale$error, "")
  expect_true(all(is.finite(scale$transform)))
  z <- c(.1, -.1)
  mapped <- p + drop(scale$transform %*% z)
  numerical <- vapply(seq_along(z), function(j) {
    up <- down <- z; up[j] <- up[j] + 1e-6; down[j] <- down[j] - 1e-6
    (fn(p + drop(scale$transform %*% up)) -
       fn(p + drop(scale$transform %*% down))) / 2e-6
  }, numeric(1))
  expect_equal(drop(crossprod(scale$transform, gr(mapped))), numerical, tolerance = 1e-8)
  expect_equal(drop(solve(scale$transform, mapped - p)), z)
  expect_gt(scale$counts[['gradient']], 0)
  for (case in list(list(fn=function(p) NA_real_, gr=function(p) rep(NA_real_,2)),
                   list(fn=function(p) 0, gr=function(p) c(0,0)))) {
    bad <- mfrm_optimizer_curvature_scale(p, case$fn, case$gr)
    expect_null(bad$transform)
    expect_true(nzchar(bad$error))
  }
  large <- mfrm_optimizer_curvature_scale(rep(0,65), function(p) stop("evaluated"), identity)
  expect_null(large$transform)
  expect_equal(unname(large$counts),c(0L,0L))
})

test_that("a small raw gradient cannot hide recoverable negative curvature", {
  fn <- function(p) (1e-6 * p[1]^2 - 1)^2 + p[2]^2
  gr <- function(p) c(4e-6 * p[1] * (1e-6 * p[1]^2 - 1), 2*p[2])
  testthat::local_mocked_bindings(
    make_param_cache = function(...) list(ensure = function(p) NULL),
    make_mfrm_direct_evaluator = function(...) list(value=fn, gradient=gr, diagnostics=function() list()),
    .package = "mfrmr")
  fit <- run_mfrm_direct_optimization(c(1,0), "MML", list(), list(model="GPCM"),
    list(), quad_points=3, maxit=400, reltol=1e-9, optimizer="BFGS")
  stages <- fit$optimizer_polish$Stages
  expect_lt(stages$TerminalGradientSupNorm[1],1e-4)
  expect_lt(stages$SmallestCurvature[1],0)
  expect_lt(fit$value,1e-10)
  expect_equal(abs(fit$par[1]),1000,tolerance=1e-5)
  expect_identical(fit$optimizer_diagnostics$ConvergenceSeverity,"pass")
  expect_identical(stages$StageLabel[stages$Selected],"curvature_rescale")
  expect_gte(stages$SmallestCurvature[stages$Selected],0)
  expect_lte(sum(stages$StageLabel=="curvature_rescale"),3)
  expect_equal(fit$optimizer_diagnostics$GradientEvaluations, sum(stages$GradientEvaluations))
})

test_that("unsuccessful negative-curvature restarts retain the point and warn", {
  fn <- function(p) 1 - 1e-6*p[1]^2 + p[2]^2
  gr <- function(p) c(-2e-6*p[1],2*p[2])
  testthat::local_mocked_bindings(
    make_param_cache = function(...) list(ensure = function(p) NULL),
    make_mfrm_direct_evaluator = function(...) list(value=fn, gradient=gr, diagnostics=function() list()),
    optim = function(par, fn, gr, ...) list(par=par, value=fn(par), convergence=0L,
      counts=c("function"=1L,"gradient"=1L)), .package = "mfrmr")
  run <- function(model="GPCM",method="MML",reltol=1e-9) {
    run_mfrm_direct_optimization(c(1,0),method,list(),list(model=model),list(),
      quad_points=3,maxit=20,reltol=reltol)
  }
  expect_warning(fit <- run(),"convergence review")
  expect_equal(fit$par,c(1,0))
  expect_identical(fit$optimizer_diagnostics$ConvergenceReason,"negative_curvature_review")
  expect_identical(fit$optimizer_diagnostics$ConvergenceSeverity,"review")
  expect_false(fit$optimizer_polish$Succeeded)
  expect_equal(sum(fit$optimizer_polish$Stages$StageLabel=="curvature_rescale"),1)
  for (other in list(run(model="RSM"),run(model="PCM"),run(method="JML"),run(reltol=1e-6))) {
    expect_identical(other$optimizer_diagnostics$ConvergenceSeverity,"pass")
    expect_false("curvature_rescale" %in% other$optimizer_polish$Stages$StageLabel)
  }
})

test_that("unsafe curvature proposals do not replace the retained point", {
  cases <- list(
    list(fn = function(p) -sum(p^2), gr = function(p) -2*p),
    list(fn = function(p) p[1]^2, gr = function(p) c(2*p[1],0)),
    list(fn = function(p) NA_real_, gr = function(p) p),
    list(fn = function(p) sum(p^2), gr = function(p) c(NA_real_,p[2])),
    list(fn = function(p) if (sum(p^2)<1) stop("outside valid domain") else sum(p^2),
         gr = function(p) 2*p),
    list(fn = function(p) -sum(p^2), gr = function(p) 2*p)
  )
  for (case in cases) {
    proposal <- mfrm_optimizer_curvature_proposal(c(1,1),case$fn,case$gr)
    expect_null(proposal$par)
    expect_true(nzchar(proposal$error))
  }
  too_large <- mfrm_optimizer_curvature_proposal(rep(1,65),function(p) stop("evaluated"),identity)
  expect_equal(unname(too_large$counts),c(0L,0L))
  expect_match(too_large$error,"1-64",fixed=TRUE)
})

test_that("a stalled optimizer restarts only from a safe curvature proposal", {
  run <- function(indefinite = FALSE, code = 0L, method = "MML",
                  model = "RSM", integration = "fixed", population = FALSE,
                  reltol = 1e-9) {
    h <- diag(c(900, if (indefinite) -15000 else 15000))
    fn <- function(p) 86000 + sum(p * (h %*% p))/2
    gr <- function(p) as.vector(h %*% p)
    testthat::local_mocked_bindings(
      make_param_cache = function(...) list(ensure = function(p) NULL),
      make_mfrm_direct_evaluator = function(...) list(value=fn, gradient=gr, diagnostics=function() list()),
      optim = function(par, fn, gr, ...) list(par=par, value=fn(par), convergence=code,
                                             counts=c("function"=1L,"gradient"=1L)),
      .package = "mfrmr"
    )
    config <- list(model = model, population_spec = list(active = population),
      estimation_control = list(mml_integration = integration))
    run_mfrm_direct_optimization(c(4e-7,-1e-7), method, list(), config,
      list(), quad_points=3, maxit=20, reltol=reltol, suppress_convergence_warning=TRUE)
  }
  good <- run()
  expect_identical(good$optimizer_diagnostics$ConvergenceSeverity,"pass")
  expect_lt(good$optimizer_diagnostics$TerminalGradientSupNorm,1e-10)
  stages <- good$optimizer_polish$Stages
  expect_identical(stages$StageLabel[stages$Selected],"curvature_restart")
  expect_gt(tail(stages$GradientEvaluations,1),1L)
  expect_gt(tail(stages$MaxParameterChange,1),0)
  bad <- run(indefinite=TRUE)
  expect_identical(bad$optimizer_diagnostics$ConvergenceSeverity,"review")
  expect_true(nzchar(tail(bad$optimizer_polish$Stages$Error,1)))
  expect_false(tail(bad$optimizer_polish$Stages$Selected,1))
  adaptive <- run(model = "GPCM", integration = "adaptive")
  expect_identical(adaptive$optimizer_diagnostics$ConvergenceSeverity, "pass")
  expect_lt(adaptive$optimizer_diagnostics$TerminalGradientSupNorm, 1e-10)
  expect_identical(tail(adaptive$optimizer_polish$Stages$StageLabel, 1), "curvature_restart")
  expect_lte(adaptive$value, good$value)
  unsafe <- run(indefinite = TRUE, model = "GPCM", integration = "adaptive")
  expect_identical(unsafe$optimizer_diagnostics$ConvergenceSeverity, "review")
  expect_true(nzchar(tail(unsafe$optimizer_polish$Stages$Error, 1)))
  expect_false(tail(unsafe$optimizer_polish$Stages$Selected, 1))
  for (model in c("RSM", "PCM", "GPCM")) {
    fixed <- run(model = model, population = TRUE)
    expect_identical(fixed$optimizer_diagnostics$ConvergenceSeverity, "pass")
    expect_lt(fixed$optimizer_diagnostics$TerminalGradientSupNorm, 1e-10)
    expect_identical(fixed$optimizer_polish$Stages$StageLabel[
      fixed$optimizer_polish$Stages$Selected], "curvature_restart")
    rejected <- run(model = model, population = TRUE, indefinite = TRUE)
    expect_identical(rejected$optimizer_diagnostics$ConvergenceSeverity, "review")
    restart <- rejected$optimizer_polish$Stages$StageLabel == "curvature_restart"
    expect_true(any(restart))
    expect_true(nzchar(rejected$optimizer_polish$Stages$Error[restart]))
    expect_false(any(rejected$optimizer_polish$Stages$Selected[restart]))
  }
  fixed_gpcm <- run(model = "GPCM")
  expect_identical(fixed_gpcm$optimizer_diagnostics$ConvergenceSeverity, "pass")
  expect_identical(fixed_gpcm$optimizer_polish$Stages$StageLabel[
    fixed_gpcm$optimizer_polish$Stages$Selected], "curvature_restart")
  for (x in list(run(code=1L),run(method="JML"),
                run(model="GPCM", integration="adaptive", population=TRUE),
                run(model="GPCM", integration="adaptive", reltol=1e-6),
                run(integration="adaptive"))) {
    expect_false("curvature_restart" %in% x$optimizer_polish$Stages$StageLabel)
  }
})
