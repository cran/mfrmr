test_that("GPCM slope numeric-boundary conditions retain audit payloads", {
  condition <- mfrmr:::new_gpcm_slope_numeric_boundary_error(
    c(800, -800, 0)
  )

  expect_s3_class(condition, "mfrmr_gpcm_slope_numeric_boundary_error")
  expect_match(
    conditionMessage(condition),
    "finite and strictly positive"
  )
  expect_equal(condition$expanded_log_slopes, c(800, -800, 0))
  expect_true(any(!is.finite(condition$expanded_slopes)))
  expect_true(any(condition$expanded_slopes <= 0))
})

test_that("direct objectives reject only typed non-representable slope trials", {
  evaluator <- list(value = function(par) {
    if (abs(par[1]) > 0.1) {
      stop(mfrmr:::new_gpcm_slope_numeric_boundary_error(par[1]))
    }
    (par[1] - 1)^2
  })
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)

  expect_equal(safe$value(0), 1)
  expect_equal(safe$value(2), 1e100)
  expect_identical(safe$rejections(), 1L)
  expect_true(is.finite(safe$value(2)))
  expect_identical(safe$rejections(), 2L)

  unrelated <- mfrmr:::make_mfrm_boundary_safe_objective(list(
    value = function(par) stop("unrelated evaluator failure"),
    gradient = function(par) stop("unrelated gradient failure")
  ))
  expect_error(unrelated$value(0), "unrelated evaluator failure")
  expect_error(unrelated$gradient(0), "unrelated gradient failure")
  expect_identical(unrelated$rejections(), 0L)
})

test_that("boundary-safe optimizers retain valid iterates and expose stalled gradients", {
  for (method in c("BFGS", "L-BFGS-B")) {
    evaluator <- list(
      value = function(par) {
        if (abs(par[1]) > 0.1) {
          stop(mfrmr:::new_gpcm_slope_numeric_boundary_error(par[1]))
        }
        (par[1] - 0.08)^2
      },
      gradient = function(par) {
        if (abs(par[1]) > 0.1) {
          stop(mfrmr:::new_gpcm_slope_numeric_boundary_error(par[1]))
        }
        2 * (par - 0.08)
      }
    )
    safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)
    opt <- stats::optim(
      par = 0, fn = safe$value, gr = safe$gradient, method = method,
      control = mfrmr:::build_mfrm_optim_control(method, 50L, 1e-12)
    )

    expect_true(is.finite(opt$value))
    expect_lte(abs(opt$par[1]), 0.1)
    gradient <- evaluator$gradient(opt$par)
    diagnostics <- mfrmr:::build_optimizer_diagnostics(
      opt, gradient, reltol = 1e-12, maxit = 50L, optimizer_method = method
    )
    if (abs(opt$par[1] - 0.08) > 1e-6) {
      # A dominating penalty can stall a line search. Its zero derivative
      # must never replace the actual terminal gradient in readiness review.
      expect_false(identical(diagnostics$ConvergenceSeverity, "pass"))
      expect_gt(diagnostics$TerminalGradientSupNorm, 1e-6)
    }
    expect_gt(safe$rejections(), 0L)
    expect_identical(safe$population_variance_rejections(), 0L)
  }
})

test_that("population variance boundaries cannot poison an MML evaluation cache", {
  fit <- make_toy_fit(method = "MML")
  config <- fit$config
  n_persons <- length(fit$prep$levels$Person)
  config$population_spec <- list(
    active = TRUE,
    design_matrix = matrix(1, nrow = n_persons, ncol = 1L),
    design_columns = "(Intercept)",
    person_lookup = seq_len(n_persons)
  )
  sizes <- mfrmr:::build_param_sizes(config)
  idx <- mfrmr:::build_indices(fit$prep, config$step_facet)
  par <- c(fit$opt$par, 0, 0)
  sigma_index <- mfrmr:::build_param_slices(sizes)$log_sigma2
  cache <- mfrmr:::make_param_cache(sizes, config, idx, is_mml = TRUE)
  evaluator <- mfrmr:::make_mfrm_direct_evaluator(
    "MML", cache, idx, config, sizes, mfrmr:::gauss_hermite_normal(7L)
  )
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)
  value <- safe$value(par)
  gradient <- evaluator$gradient(par)
  params <- cache$params()

  for (log_sigma2 in c(800, -800)) {
    trial <- par
    trial[sigma_index] <- log_sigma2
    error <- tryCatch(cache$ensure(trial), error = identity)
    expect_s3_class(error, "mfrmr_population_variance_numeric_boundary_error")
    expect_identical(error$log_sigma2, log_sigma2)
    expect_identical(error$sigma2, exp(log_sigma2))
    expect_identical(cache$params(), params)
    # Repeated rejected points must never be mistaken for successful cache hits.
    expect_equal(safe$value(trial), safe$penalty)
    expect_equal(safe$value(trial), safe$penalty)
    expect_identical(safe$gradient(trial), numeric(length(trial)))
    expect_error(evaluator$gradient(trial),
                 class = "mfrmr_population_variance_numeric_boundary_error")
    expect_identical(safe$value(par), value)
    expect_identical(evaluator$gradient(par), gradient)
  }
  expect_identical(safe$population_variance_rejections(), 4L)
  expect_identical(safe$rejections(), 0L)
  expect_error(mfrmr:::run_mfrm_direct_optimization(
    trial, "MML", idx, config, sizes, quad_points = 7L,
    maxit = 2L, reltol = 1e-6
  ), class = "mfrmr_population_variance_numeric_boundary_error")
  for (log_sigma2 in c(-10, 10)) {
    trial <- par
    trial[sigma_index] <- log_sigma2
    cache$ensure(trial)
    expect_identical(cache$params()$population$sigma2, exp(log_sigma2))
  }
})

test_that("finite GPCM slopes can produce explicitly rejected adaptive moments", {
  kernel <- mfrmr:::mfrmr_adaptive_person_kernel(
    list(person=1L,score_k=1L,weight=1,step_idx=1L,slope_idx=1L),
    list(model="GPCM",n_cat=2L),
    list(slopes=1e200,steps_mat=matrix(0,1,1)),0,0,1)
  expect_error(kernel(0), class="mfrmr_adaptive_quadrature_numeric_error")
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(list(
    value=function(p)kernel(p)$information,gradient=function(p)kernel(p)$score))
  expect_identical(safe$value(0),safe$penalty)
  expect_identical(safe$gradient(0),0)
  expect_identical(safe$adaptive_quadrature_rejections(),1L)
  expect_identical(safe$rejections(),0L)
})

test_that("failed adaptive trials do not cache old likelihoods or qualify invalid starts", {
  fit <- suppressWarnings(fit_mfrm(load_mfrmr_data("example_core"),
    "Person",c("Rater","Criterion"),"Score",model="GPCM",
    step_facet="Criterion",slope_facet="Criterion",quad_points=5L,maxit=2L))
  cfg <- fit$config
  cfg$estimation_control$mml_integration <- "adaptive"
  sizes <- mfrmr:::build_param_sizes(cfg)
  idx <- mfrmr:::build_indices(fit$prep,cfg$step_facet,cfg$slope_facet)
  evaluator <- mfrmr:::make_mfrm_direct_evaluator("MML",
    mfrmr:::make_param_cache(sizes,cfg,idx,is_mml=TRUE),idx,cfg,sizes,
    mfrmr:::gauss_hermite_normal(5L))
  safe <- mfrmr:::make_mfrm_boundary_safe_objective(evaluator)
  valid <- rep(0,length(fit$opt$par)); invalid <- valid
  invalid[mfrmr:::build_param_slices(sizes)$log_slopes[1]] <- 580
  expect_true(all(is.finite(mfrmr:::expand_params(invalid,sizes,cfg)$slopes)))
  expected <- evaluator$value(valid); gradient <- evaluator$gradient(valid)
  for(i in 1:2) {
    expect_identical(safe$value(invalid),safe$penalty)
    expect_identical(safe$gradient(invalid),numeric(length(invalid)))
    expect_error(evaluator$value(invalid),class="mfrmr_adaptive_quadrature_numeric_error")
  }
  expect_identical(safe$adaptive_quadrature_rejections(),2L)
  expect_identical(evaluator$value(valid),expected)
  expect_identical(evaluator$gradient(valid),gradient)
  expect_error(mfrmr:::run_mfrm_direct_optimization(invalid,"MML",idx,cfg,sizes,
    quad_points=5L,maxit=2L,reltol=1e-6),class="mfrmr_adaptive_quadrature_numeric_error")
})
