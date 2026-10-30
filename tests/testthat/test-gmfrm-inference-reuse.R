gmfrm_auto_reuse_fixture <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- readRDS(test_path("fixtures", "gmfrm-joint-information.rds"))$data
      fit <- suppressWarnings(fit_mfrm(d, "Person", c("Task", "Rater"), "Score",
        model = "GPCM", method = "MML", slope_facet = c("Task", "Rater"),
        step_facet = "Rater", noncenter_facet = "Rater", mml_engine = "em",
        gpcm_mml_identification = "fixed_standard_normal", mml_integration = "auto",
        mml_accuracy_control = list(max_points = 241L), maxit = 500L,
        rating_min = 0, rating_max = 2, category_policy = "preserve"))
      saved <<- list(fit = fit, inference = mfrm_gpcm_product_inference(fit))
    }
    saved
  }
})

test_that("automatic qualification serves slope and location intervals without recalculation", {
  calls <- integer()
  original <- compute_mml_parameter_covariance
  local_mocked_bindings(compute_mml_parameter_covariance = function(res, ...) {
    calls <<- c(calls, res$config$estimation_control$quad_points)
    original(res, ...)
  }, .package = "mfrmr")
  withr::defer(mfrmr_gmfrm_inference_reuse(NULL))
  x <- gmfrm_auto_reuse_fixture(); fit <- x$fit
  checked <- calls
  expect_gt(length(checked), 0L)
  expect_true(x$inference$check$eligible)
  expect_warning(slopes <- confint(fit, parm = "slopes"), "Experimental")
  expect_warning(locations <- mfrm_facet_intervals(fit, "Rater"), "Experimental")
  expect_identical(calls, checked)
  expect_identical(fit, unserialize(serialize(fit, NULL)))
  expect_false(any(vapply(fit, is.environment, FALSE)))
  mfrmr_gmfrm_inference_reuse(NULL)
  fresh <- mfrm_gpcm_product_inference(fit)
  expect_identical(fresh, x$inference)
  expect_identical(tail(calls, 2L), c(fit$config$estimation_control$quad_points,
    2L * fit$config$estimation_control$quad_points - 1L))
  local_mocked_bindings(mfrm_gpcm_product_inference = function(...) fresh,
    .package = "mfrmr")
  expect_identical(suppressWarnings(confint(fit, parm = "slopes")), slopes)
  expect_identical(suppressWarnings(mfrm_facet_intervals(fit, "Rater")), locations)
  expect_identical(fit, x$fit)
})

test_that("reuse requires the entire unchanged fit and current memory budget", {
  x <- gmfrm_auto_reuse_fixture()
  withr::defer(mfrmr_gmfrm_inference_reuse(NULL))
  for (field in names(x$fit)) {
    mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
    expect_identical(mfrmr_gmfrm_inference_reuse(x$fit), x$inference)
    changed <- x$fit; changed[field] <- list(NULL)
    expect_null(mfrmr_gmfrm_inference_reuse(changed), info = field)
    expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
  }
  mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
  changed <- x$fit; changed$opt$par[1] <- changed$opt$par[1] + .01
  expect_null(mfrmr_gmfrm_inference_reuse(changed))
  changed <- x$fit; changed$summary$Converged <- FALSE
  mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
  expect_false(mfrm_gpcm_product_inference(changed)$check$eligible)
  backend <- getOption("mfrmr.use_cpp11_backend", TRUE)
  mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
  withr::local_options(mfrmr.use_cpp11_backend = !isTRUE(backend))
  expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
  withr::local_options(mfrmr.use_cpp11_backend = backend)
  mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
  withr::local_options(mfrmr.max_information_bytes = 1)
  expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
  mfrmr_gmfrm_inference_reuse(x$fit, x$inference)
  expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
  limited <- mfrm_gpcm_product_inference(x$fit)
  expect_false(limited$check$eligible)
  expect_match(limited$check$review, "workspace")
  withr::local_options(mfrmr.max_information_bytes = NA_real_)
  expect_error(mfrm_gpcm_product_inference(x$fit), "finite positive")
})

test_that("warnings and refused qualification never create reusable automatic evidence", {
  x <- gmfrm_auto_reuse_fixture()
  withr::defer(mfrmr_gmfrm_inference_reuse(NULL))
  args <- list(mml_engine = "em", quad_points = 31L,
    mml_accuracy_control = mfrmr_gmfrm_accuracy_control())
  local_mocked_bindings(mfrm_fit_product_slopes = function(...) x$fit,
    mfrmr_gmfrm_accuracy_check = function(...) {
      warning("qualification warning")
      list(passed = TRUE, retry = FALSE, inference = x$inference, per_person = NULL,
        metrics = c(MaxLogMarginalDifference = 0, MaxEAPDifference = 0,
          MaxPosteriorSDDifference = 0), detail = "checked")
    }, .package = "mfrmr")
  expect_warning(fit <- mfrm_fit_product_slopes_auto(args, character()), "qualification warning")
  expect_null(mfrmr_gmfrm_inference_reuse(fit))
  refused <- x$inference; refused$check$eligible <- FALSE
  mfrmr_gmfrm_inference_reuse(x$fit, refused)
  expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
  local_mocked_bindings(mfrmr_gmfrm_accuracy_check = function(...) stop("check failed"),
    .package = "mfrmr")
  expect_error(mfrm_fit_product_slopes_auto(args, character()), "check failed")
  expect_null(mfrmr_gmfrm_inference_reuse(x$fit))
})
