test_that("additive JML information matches independent native-gradient differences", {
  d <- load_mfrmr_data("example_core")
  for (model in c("RSM", "PCM")) for (anchored in c(FALSE, TRUE)) {
    a <- if (anchored) data.frame(Facet="Rater", Level=as.character(d$Rater[1]), Anchor=.2) else NULL
    f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method="JML",
      model=model, step_facet=if (model == "PCM") "Criterion" else NULL,
      anchors=a, optimizer="BFGS", reltol=1e-12, maxit=500)
    expect_equal(nrow(f$config$replay_inputs$anchors), as.integer(anchored))
    cfg <- f$config
    idx <- build_indices(f$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
    # Exercise the weighted likelihood too, away from an optimum.
    idx$weight <- rep(c(.5, 1, 2), length.out=length(idx$score_k))
    sizes <- build_param_sizes(cfg); par <- f$opt$par + .003
    fn <- function(p) mfrm_loglik_jml(p, idx, cfg, sizes)
    gr <- function(p) mfrm_grad_jml(p, idx, cfg, sizes)
    H <- mfrm_jml_polish_information(par, idx, cfg, sizes)
    reference <- optimHess(par, fn, gr, control=list(ndeps=rep(1e-5, length(par))))
    expect_equal(unname(as.matrix(H)), unname(reference), tolerance=1e-7)
  }
})

test_that("bounded information finishing preserves the likelihood and rejects unsafe proposals", {
  H <- Matrix::Diagonal(x=c(900, 15000))
  fn <- function(p) 86000 + sum(p * (H %*% p)) / 2
  gr <- function(p) as.numeric(H %*% p)
  start <- c(4e-7, -1e-7)
  result <- mfrm_jml_information_polish(start, fn, gr, function(p) H)
  expect_identical(result$convergence, 0L)
  expect_identical(result$steps, 1L)
  expect_lt(max(abs(result$gradient)), 1e-10)
  expect_lte(result$value, fn(start))
  expect_lte(result$residual, 1e-8)
  bad <- list(function(p) Matrix::Diagonal(x=c(900, -15000)),
    function(p) Matrix::Diagonal(x=c(900, 0)),
    function(p) Matrix::Diagonal(x=c(1e-9, 1e-9)))
  for (information in bad) {
    rejected <- mfrm_jml_information_polish(start, fn, gr, information)
    expect_identical(rejected$par, start)
    expect_identical(rejected$convergence, 1L)
    expect_true(nzchar(rejected$error))
  }
  loss_start <- c(.001, .001)
  rejected <- mfrm_jml_information_polish(loss_start, function(p) -fn(p), gr, function(p) H)
  expect_identical(rejected$par, loss_start)
  expect_match(rejected$error, "likelihood loss")
  invalid <- mfrm_jml_information_polish(start, function(p) NA_real_,
    function(p) c(0, 0), function(p) H)
  expect_identical(invalid$convergence, 1L)
  expect_match(invalid$error, "finite objective")
})

test_that("the public JML fit records finishing without changing its equation", {
  d <- load_mfrmr_data("example_core")
  args <- list(data=d, person="Person", facets=c("Rater", "Criterion"), score="Score",
    method="JML", optimizer="BFGS", maxit=500, reltol=1e-12)
  reference <- do.call(fit_mfrm, args)
  stalled <- reference$opt$par + 1e-4
  # Simulate optim's premature objective stop; all proposed points still use
  # the actual native likelihood/gradient and exact model information.
  local_mocked_bindings(optim=function(par, fn, gr, ...) list(par=stalled,
    value=fn(stalled), counts=c("function"=1L, "gradient"=1L), convergence=0L),
    .package="mfrmr")
  result <- do.call(fit_mfrm, args)
  stages <- result$opt$optimizer_polish$Stages
  expect_gt(stages$TerminalGradientSupNorm[1], 1e-4)
  expect_identical(stages$StageLabel[stages$Selected], "jml_information_polish")
  expect_identical(stages$Method[stages$Selected], "Newton")
  expect_lte(stages$InformationSteps[stages$Selected], 3L)
  expect_lt(mfrm_fair_refit_check(result, 1e-4)$max_gradient, 1e-4)
  expect_equal(result$opt$par, reference$opt$par, tolerance=1e-4)
  expect_lte(result$opt$value, reference$opt$value + 1e-9)
  expect_equal(result$config$replay_inputs$optimizer, "BFGS")
  expect_equal(result$opt$optimizer_diagnostics$FunctionEvaluations, sum(stages$FunctionEvaluations))
  expect_equal(result$opt$optimizer_diagnostics$GradientEvaluations, sum(stages$GradientEvaluations))
})

test_that("finishing preserves iteration limits and cannot hide a frozen gradient", {
  d <- load_mfrmr_data("example_core")
  f <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method="JML")
  cfg <- f$config; sizes <- build_param_sizes(cfg)
  idx <- build_indices(f$prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs)
  code <- 0L
  local_mocked_bindings(optim=function(par, fn, gr, ...) list(par=par, value=fn(par),
    counts=c("function"=1L, "gradient"=1L), convergence=code), .package="mfrmr")
  run <- function() run_mfrm_direct_optimization(f$opt$par + .01, "JML", idx,
    cfg, sizes, quad_points=3, maxit=20, reltol=1e-12, suppress_convergence_warning=TRUE)
  code <- 1L
  limited <- run()
  expect_false(any(limited$optimizer_polish$Stages$StageLabel == "jml_information_polish"))
  code <- 0L; idx$score_k[idx$person == 1L] <- 0L
  extreme <- run(); stages <- extreme$optimizer_polish$Stages
  expect_false(extreme$optimizer_polish$Succeeded)
  expect_gt(max(abs(extreme$optimizer_diagnostics$TerminalGradientSupNorm)), 1e-4)
  expect_false(tail(stages$Selected, 1))
})

test_that("block Newton finishing holds inactive coordinates and checks the full gradient", {
  H <- Matrix::Diagonal(x=c(0,900,15000))
  fn <- function(p) 86000+sum(p*(H%*%p))/2
  gr <- function(p) as.numeric(H%*%p)
  start <- c(28,4e-7,-1e-7)
  answer <- mfrm_jml_information_polish(start,fn,gr,function(p)H,active=2:3)
  expect_identical(answer$par[1],start[1])
  expect_identical(answer$fixed_coordinates,1L)
  expect_identical(answer$convergence,0L)
  expect_lt(max(abs(answer$gradient)),1e-10)
  expect_lte(answer$value,fn(start))
  # An unresolved gradient in a held coordinate must still block convergence.
  bad <- mfrm_jml_information_polish(start,function(p)fn(p)+p[1],
    function(p)gr(p)+c(1,0,0),function(p)H,active=2:3)
  expect_identical(bad$convergence,1L)
  expect_identical(bad$par[1],start[1])
  expect_error(mfrm_jml_information_polish(start,fn,gr,function(p)H,active=c(1,1)),"distinct")
})

test_that("public block finishing preserves extreme limits and fixed Person coordinates", {
  d <- load_mfrmr_data("example_core")
  d$Score[d$Person=="P001"] <- 1
  d$Score[d$Person=="P002"] <- 4
  args <- list(data=d,person="Person",facets=c("Rater","Criterion"),score="Score",
    method="JML",optimizer="BFGS",rating_min=1,rating_max=4,maxit=500,reltol=1e-12)
  reference <- do.call(fit_mfrm,args)
  extreme <- match(c("P001","P002"),reference$config$person_levels)
  theta <- build_param_slices(build_param_sizes(reference$config))$theta
  fixed <- theta[extreme]
  stalled <- reference$opt$par
  stalled[-fixed] <- stalled[-fixed]+1e-4
  local_mocked_bindings(optim=function(par,fn,gr,...)list(par=stalled,value=fn(stalled),
    counts=c("function"=1L,"gradient"=1L),convergence=0L),.package="mfrmr")
  fit <- do.call(fit_mfrm,args)
  stages <- fit$opt$optimizer_polish$Stages
  expect_identical(stages$StageLabel[stages$Selected],"jml_information_polish")
  expect_identical(fit$opt$par[fixed],stalled[fixed])
  expect_identical(fit$facets$person$PrimaryEstimate[extreme],c(-Inf,Inf))
  expect_lt(mfrm_fair_refit_check(fit,1e-4)$max_gradient,1e-4)
  expect_lte(fit$opt$value,reference$opt$value+1e-8)
})
