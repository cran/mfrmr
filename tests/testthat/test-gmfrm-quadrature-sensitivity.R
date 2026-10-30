gmfrm_quadrature_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      d <- readRDS(test_path("fixtures","gmfrm-joint-information.rds"))$data
      # Both families deliberately use the same levels and arbitrary owner names.
      d$Task <- sub("^.*?([0-9]+)$","\\1",as.character(d$Task))
      d$Rater <- sub("^.*?([0-9]+)$","\\1",as.character(d$Rater))
      names(d)[match(c("Task","Rater"),names(d))] <- c("Task label","Assessor label")
      fit <- fit_mfrm(d,person="Person",facets=c("Task label","Assessor label"),score="Score",
        model="GPCM",method="MML",slope_facet=c("Task label","Assessor label"),
        step_facet="Assessor label",noncenter_facet="Assessor label",
        mml_engine="em",gpcm_mml_identification="fixed_standard_normal",
        quad_points=5L,maxit=500L,em_score_tol=1e-6,
        rating_min=0,rating_max=2,category_policy="preserve")
      review <- suppressWarnings(mml_quadrature_sensitivity(fit,d,
        quad_points=c(5L,31L),theta_points=21L,adaptive_quad_points=c(15L,31L)))
      cached <<- list(data=d,fit=fit,review=review)
    }
    cached
  }
})

test_that("two-family quadrature refits preserve the model and interval evidence", {
  x <- gmfrm_quadrature_fixture(); review <- x$review
  expect_s3_class(review,"mfrm_quadrature_sensitivity")
  expect_identical(review$fits$q5,x$fit)
  expect_identical(review$fits$q31$gmfrm$controls,x$fit$gmfrm$controls)
  expect_identical(review$fits$q31$config$replay_inputs$rating_max,2)
  expect_true(mfrmr_gqs_same_prepared_data(x$fit,review$fits$q31))
  expect_true(all(review$runs$EstimationConverged))
  expect_false(any(review$runs$InferenceReady))
  expect_true(all(is.finite(review$summary$ProbabilityMaxAbsChange)))
  expect_true(all(is.finite(review$summary$RawSlopeSEMaxAbsChange)))
  expect_true(all(is.na(review$summary$EAPMaxAbsChange)))
  expect_true(all(is.na(review$summary$PosteriorSDMaxAbsChange)))
  expect_identical(review$settings$person_score_comparison,"unavailable_two_family")
  expect_true(all(is.na(review$intervals$q5)))
  expect_true(all(is.finite(review$intervals$q31)))
  expect_identical(review$summary$SlopeIntervalEligibilityChanged,c(0L,6L))
  fine <- attr(review$intervals$q31,"diagnostics")
  expect_equal(review$slopes$SlopeOwner,rep(fine$SlopeOwner,2L))
  expect_equal(review$slopes$ScaleReference,rep(fine$ScaleReference,2L))
  expect_false(anyDuplicated(review$slopes[c("Nodes","SlopeOwner","SlopeFacet")])>0)
  expect_true(any(grepl("Experimental two-family",review$conditions$Message)))
  expect_output(print(review),"does not compare two-family Person scores")
  expect_output(print(summary(review)),"does not compare two-family Person scores")
  expect_s3_class(apa_table(review,digits=5),"apa_table")
  expect_s3_class(plot(review$intervals$q31,draw=FALSE),"ggplot")
  expect_equal(as.data.frame(review),review$summary)
  expect_equal(summary(review)$intervals,review$intervals)
  path <- tempfile(fileext=".rds"); saveRDS(review,path)
  withr::defer(unlink(path))
  expect_identical(readRDS(path)$intervals,review$intervals)
  expect_identical(readRDS(path)$quadrature_review,review$quadrature_review)
})

test_that("nonconverged two-family refits keep missing intervals and their reason", {
  x <- gmfrm_quadrature_fixture()
  args <- mfrmr_gqs_refit_arguments(x$fit,x$data,5L); args$maxit <- 1L
  fit <- do.call(fit_mfrm,args)
  review <- suppressWarnings(gpcm_mml_quadrature_sensitivity(fit,x$data,
    quad_points=c(5L,7L),theta_points=21L))
  expect_false(any(review$runs$EstimationConverged))
  expect_false(any(review$slopes$CIEligible))
  expect_true(all(is.na(review$slopes$CI_Lower)))
  expect_match(review$slopes$InferenceReview[1],"EM must meet")
  expect_true(all(is.finite(review$summary$ProbabilityMaxAbsChange)))
  expect_identical(review$fits$q7$gmfrm$controls$maxit,1L)
})

test_that("quadrature probabilities and diagnostic SEs use both slope owners", {
  fit <- gmfrm_quadrature_fixture()$review$fits$q31
  got <- mfrmr_gqs_design_probabilities(fit,c(-4,4),21L)
  grid <- got$cells[rep(seq_len(nrow(got$cells)),each=21L),,drop=FALSE]
  grid$Theta <- rep(got$theta,nrow(got$cells))
  expected <- mfrm_curve_intervals(fit,grid)$table$Estimate
  expect_equal(got$values,expected,tolerance=1e-12)
  info <- compute_mml_parameter_covariance(fit)
  # Independent map: first family log slopes sum to zero, second is free.
  J <- matrix(0,6,13); J[1:3,9:10] <- rbind(c(1,0),c(0,1),c(-1,-1))
  J[4:6,11:13] <- diag(3)
  expected_se <- fit$slopes$Estimate*sqrt(diag(J%*%info$cov%*%t(J)))
  tab <- gmfrm_quadrature_fixture()$review$slopes
  expect_equal(tab$RawObservedInformationSE[tab$Nodes==31],expected_se,tolerance=1e-10)
  reordered <- fit; reordered$slopes <- fit$slopes[6:1,]
  expect_equal(mfrmr_gqs_max_named_difference(mfrmr_gqs_measurement_parameters(fit),
    mfrmr_gqs_measurement_parameters(reordered),"measurement"),0)
  expect_identical(mfrmr_gqs_slope_keys(reordered),rev(mfrmr_gqs_slope_keys(fit)))
})

test_that("two-family sensitivity shares fresh information without reusing eligibility", {
  fit <- gmfrm_quadrature_fixture()$review$fits$q31
  original <- compute_mml_parameter_covariance
  orders <- integer()
  local_mocked_bindings(compute_mml_parameter_covariance = function(x, ...) {
    orders <<- c(orders, x$config$estimation_control$quad_points)
    original(x, ...)
  }, .package = "mfrmr")
  expected <- suppressWarnings(confint(fit, parm = "slopes"))
  orders <- integer()
  actual <- suppressWarnings(mfrmr_gqs_raw_information(fit))
  expect_identical(orders, c(31L, 61L))
  expect_identical(actual$intervals, expected)
  expect_true(all(is.finite(actual$slope_se)))
  expect_true(all(attr(actual$intervals, "diagnostics")$CIEligible))

  # Refusal before information must preserve raw SEs without borrowing the
  # preceding fit's permission to report intervals.
  changed <- fit
  changed$summary$Converged <- FALSE
  expected <- confint(changed, parm = "slopes")
  orders <- integer()
  refused <- mfrmr_gqs_raw_information(changed)
  expect_identical(orders, 31L)
  expect_identical(refused$intervals, expected)
  expect_identical(refused$slope_se, actual$slope_se)
  expect_false(any(attr(refused$intervals, "diagnostics")$CIEligible))
  expect_match(attr(refused$intervals, "diagnostics")$InferenceReview[1], "EM must meet")
  expect_true(fit$summary$Converged)
  imported <- changed
  class(imported) <- c("mfrm_imported_fit", class(imported))
  expect_error(mfrmr_gqs_raw_information(imported), "require a native GPCM fit")
})

test_that("two-family replay controls and invalid integration requests stay explicit", {
  x <- gmfrm_quadrature_fixture()
  args <- mfrmr_gqs_refit_arguments(x$fit,x$data,31L)
  expect_identical(args$em_score_tol,1e-6)
  expect_null(args$reltol)
  expect_null(args$package_version)
  changed <- x$fit; changed$config$replay_inputs$em_score_tol <- 2e-6
  expect_identical(mfrmr_gqs_refit_arguments(changed,x$data,31L)$em_score_tol,2e-6)
  for (fun in list(mml_quadrature_sensitivity,gpcm_mml_quadrature_sensitivity)) {
    expect_error(fun(x$fit,x$data,quad_points=c(5L,31L),adaptive_quad_points=c(15L,15L)),
      "adaptive_quad_points")
  }
  args$mml_integration <- "adaptive"
  expect_error(do.call(fit_mfrm,args),"requires MML with the direct engine",fixed=TRUE)
  changed <- x$review$fits$q31
  changed$prep$data$score_k[1] <- changed$prep$data$score_k[1]+1L
  expect_false(mfrmr_gqs_same_prepared_data(x$fit,changed))
})

test_that("two-family adaptive review retains model identity and diagnostic-only outputs", {
  x <- gmfrm_quadrature_fixture(); review <- x$review; q <- review$quadrature_review
  expect_equal(nrow(q), 4L * length(x$fit$prep$levels$Person))
  expect_true(all(q$Status == "computed"))
  expect_identical(summary(review)$quadrature_review,q)
  expect_identical(summary(review)$overview$ReadinessEffect,"none_diagnostic_only")
  expect_identical(review$settings$adaptive_quad_points,c(15L,31L))
  expect_true(all(is.na(review$summary$EAPMaxAbsChange)))
  for (nodes in c(5L,31L)) {
    fit <- review$fits[[paste0("q",nodes)]]
    rows <- q[q$FixedNodes==nodes & q$AdaptiveNodes==31L,]
    expect_identical(rows$Person,fit$prep$levels$Person)
    expect_equal(sum(rows$FixedLogMarginal),-fit$opt$value,tolerance=1e-9)
    expect_identical(fit$config$estimation_control$mml_integration,"fixed")
  }
  expect_gt(max(abs(q$LogMarginalChange[q$FixedNodes==5L])),1e-3)
  # Replaying summaries never performs inference or changes the calibration.
  expect_identical(review$fits$q5,x$fit)
  expect_error(predict_mfrm_units(x$fit,x$data),"source calibration did not pass")
})

# Independent literal probability equation and continuous posterior integrals.
# This does not call the package probability kernel or GH integration rule.
gmfrm_quadrature_reference <- function(fit, person, mode, scale) {
  spec <- fit$gmfrm$specification; owners <- spec$slope_facets
  data <- spec$data[as.character(spec$data[[spec$person]])==person,]
  log_joint <- function(theta) {
    sum(vapply(seq_len(nrow(data)),function(i) {
      location <- sum(vapply(owners,function(owner) {
        tab <- fit$facets$others
        tab$Estimate[tab$Facet==owner & tab$Level==data[[owner]][i]]
      },0))
      slope <- prod(vapply(owners,function(owner) {
        tab <- fit$slopes
        tab$Estimate[tab$SlopeOwner==owner & tab$SlopeFacet==data[[owner]][i]]
      },0))
      tab <- fit$steps[fit$steps$StepFacet==data[[owners[2]]][i],]
      steps <- tab$Estimate[order(tab$Step)]
      eta <- c(0,cumsum(slope*(theta-location-steps)))
      eta[data[[spec$score]][i]+1L]-max(eta)-log(sum(exp(eta-max(eta))))
    },0))+dnorm(theta,log=TRUE)
  }
  height <- log_joint(mode)
  moments <- vapply(0:2,function(power) integrate(function(z) vapply(z,function(x)
    x^power*exp(log_joint(mode+scale*x)-height)*scale,0),-Inf,Inf,
    rel.tol=1e-10,abs.tol=1e-12)$value,0)
  c(log_marginal=height+log(moments[1]),eap=mode+scale*moments[2]/moments[1],
    sd=scale*sqrt(moments[3]/moments[1]-(moments[2]/moments[1])^2))
}

test_that("both slope families match continuous integrals at unchanged parameters", {
  review <- gmfrm_quadrature_fixture()$review
  for (nodes in c(5L,31L)) {
    fit <- review$fits[[paste0("q",nodes)]]
    q <- review$quadrature_review
    q <- q[q$FixedNodes==nodes & q$AdaptiveNodes==31L,]
    selected <- unique(c(1L,which.min(q$LocalPosteriorSD),which.max(abs(q$PosteriorMode))))
    for (i in selected) {
      ref <- gmfrm_quadrature_reference(fit,q$Person[i],q$PosteriorMode[i],q$LocalPosteriorSD[i])
      got <- unlist(q[i,c("AdaptiveLogMarginal","AdaptiveEAP","AdaptivePosteriorSD")],use.names=FALSE)
      expect_lt(max(abs(got-ref)),1e-7)
    }
  }
})

test_that("adaptive review uses observed crossings in an incomplete two-family design", {
  x <- gmfrm_quadrature_fixture(); fit <- x$review$fits$q31
  owners <- fit$config$slope_facet; data <- x$data
  data <- data[!(data[[owners[1]]]=="1" & data[[owners[2]]]=="3"),]
  problem <- mfrm_gmfrm_problem(data,2L,gauss_hermite_normal(31L),slope_facets=owners)
  fit <- mfrm_gmfrm_fit_result(problem,mfrm_gmfrm_em(problem,start=fit$opt$par,maxit=1L))
  config <- fit$config
  expect_lt(nrow(config$gpcm_spec$cells),prod(lengths(config$gpcm_spec$component_levels)))
  idx <- build_indices(fit$prep,config$step_facet,config$slope_facet,
    config$interaction_specs,gpcm_spec=config$gpcm_spec)
  q <- mfrmr_adaptive_quadrature_review(idx,config,
    expand_params(fit$opt$par,build_param_sizes(config),config),gauss_hermite_normal(31L),
    fit$prep$levels$Person,c(15L,31L))
  row <- q[q$AdaptiveNodes==31L,][1,]
  ref <- gmfrm_quadrature_reference(fit,row$Person,row$PosteriorMode,row$LocalPosteriorSD)
  expect_lt(max(abs(unlist(row[c("AdaptiveLogMarginal","AdaptiveEAP","AdaptivePosteriorSD")],
    use.names=FALSE)-ref)),1e-7)
  expect_equal(sum(q$FixedLogMarginal[q$AdaptiveNodes==31L]),-fit$opt$value,tolerance=1e-9)
})
