mml_two_binary_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      patterns <- rbind(c(0,0),c(0,1),c(1,0),c(1,1))
      responses <- patterns[rep(1:4, c(35,15,20,30)), ]
      data <- data.frame(Person=rep(paste0("P",1:100),each=2),
        Item=rep(c("A","B"),100),Score=as.vector(t(responses)))
      cached <<- fit_mfrm(data,"Person","Item","Score",model="GPCM",method="MML",
        step_facet="Item",slope_facet="Item",rating_min=0,rating_max=1,
        category_policy="preserve",quad_points=61,maxit=400,reltol=1e-10)
    }
    cached
  }
})

test_that("two binary responses admit distinct calibrations with the same distribution", {
  # Independent continuous-normal calculation: symmetry fixes each marginal
  # at 1/2, leaving just the joint success probability to match.
  joint <- function(log_slope, sigma) integrate(function(z)
    plogis(exp(log_slope)*sigma*z)*plogis(exp(-log_slope)*sigma*z)*dnorm(z),
    -Inf,Inf,rel.tol=1e-11)$value
  first <- joint(0,1)
  other_sd <- uniroot(function(s) joint(.7,s)-first,c(.1,10),tol=1e-10)$root
  second <- joint(.7,other_sd)
  expect_gt(abs(other_sd-1),.01)
  expect_equal(c(first,.5-first,.5-first,first),
    c(second,.5-second,.5-second,second),tolerance=1e-10)
  expect_equal(exp(.7)*exp(-.7),1)
})

test_that("response-distribution dimension overrides positive numerical curvature", {
  fit <- mml_two_binary_fixture()
  check <- mfrm_mml_inference_identification(fit)
  expect_true(check$excluded)
  expect_equal(check$free_dimension,4)
  expect_equal(check$probability_dimension_bound,3)
  expect_match(check$detail,"more Persons")
  expect_true(all(is.finite(coef(fit))))
  before <- fit
  # Refuse before numerical differentiation, even if a caller supplies a
  # positive old information result or an old fit lacks the optional audit.
  local_mocked_bindings(mfrm_refine_mml_information=function(...) stop("No refinement"),
    .package="mfrmr")
  info <- compute_mml_parameter_covariance(fit)
  expect_null(info$cov)
  expect_identical(info$status,"unavailable")
  expect_match(info$detail,"at most 3 independent probabilities for 4 free parameters")
  expect_false(mfrm_gpcm_bootstrap_check(fit,list(information=info))$eligible)
  positive <- list(status="evaluated_diagnostic_only",free_dimension=4,
    evaluation_summary=data.frame(ReevaluatedObjective=fit$opt$value,GradientMaxAbs=0),
    eigenvalue_summary=data.frame(Smallest=1,AbsoluteScale=2))
  decision <- mfrm_ic_fit_check(fit,mfrm_extract_fit_ic_contract(fit),positive)
  expect_false(decision$eligible)
  expect_match(decision$review,"cannot identify")
  fit$config$estimability_audit <- NULL
  expect_true(mfrm_mml_inference_identification(fit)$excluded)
  fit$config$estimation_control$mml_integration <- "adaptive"
  expect_true(mfrm_mml_inference_identification(fit)$excluded)
  expect_identical(mml_two_binary_fixture(),before)
})

test_that("unidentified GPCM outputs retain points, unavailable rows and reasons", {
  fit <- mml_two_binary_fixture()
  contrast <- matrix(c(.001,-.001),1,dimnames=list("scaled log ratio",c("A","B")))
  for (args in list(list(),list(contrasts=contrast),list(scale="standardized"),
      list(method="sandwich"),list(contrasts=contrast,contrast_scale="difference"))) {
    ci <- do.call(confint,c(list(object=fit,parm="slopes"),args))
    tab <- attr(ci,"diagnostics")
    expect_true(all(is.finite(tab$Estimate)))
    expect_true(all(is.na(ci)))
    expect_false(any(tab$CIEligible))
    expect_match(tab$InferenceReview[1],"cannot identify")
  }
  grid <- data.frame(Theta=c(0,.5),Item=c("A","B"))
  for (type in c("probability","information")) {
    curves <- mfrm_curve_intervals(fit,grid,type=type)
    expect_true(all(is.finite(curves$table$Estimate)))
    expect_false(any(curves$table$CIEligible))
    expect_match(curves$table$InferenceReview[1],"cannot identify")
  }
  result <- mfrm_results(fit,intervals=list(slopes=ci),compute="never")
  expect_identical(mfrm_report(result)$tables$gpcm_slopes_intervals$InferenceReview,
    attr(ci,"diagnostics")$InferenceReview)
  expect_identical(unserialize(serialize(result,NULL)),result)
  expect_error(confint(fit,parm="slopes",method="profile",slope="A"),"cannot identify")
})

test_that("the dimension bound distinguishes designs and is not a sample-size gate", {
  fit <- mml_two_binary_fixture()
  fit$config$estimability_audit <- NULL
  # Known N(0,1) removes the mean/variance coordinates; this bound alone no
  # longer excludes the two-item model. It does not certify identification.
  fixed <- fit
  fixed$config$population_spec$active <- FALSE
  fixed$opt$par <- fixed$opt$par[seq_len(sum(unlist(build_param_sizes(fixed$config))))]
  expect_false(mfrm_mml_inference_identification(fixed)$excluded)
  # Distinct population covariates give distinct distributions on the same
  # observed roster and must not be pooled as one probability simplex.
  covariates <- fit
  covariates$config$population_spec$design_matrix[1,1] <- 2
  expect_false(mfrm_mml_inference_identification(covariates)$excluded)
  # Repeated independent rating events expand each Person's response support.
  repeated <- fit
  repeated$prep$data <- rbind(fit$prep$data,fit$prep$data)
  expect_false(mfrm_mml_inference_identification(repeated)$excluded)
  # Order of Persons/rows and missing unused factor labels cannot affect it.
  reordered <- fit
  reordered$prep$data <- fit$prep$data[nrow(fit$prep$data):1,]
  expect_identical(mfrm_mml_inference_identification(reordered),mfrm_mml_inference_identification(fit))
  subset <- fit
  subset$prep$data <- fit$prep$data[fit$prep$data$Person %in% levels(fit$prep$data$Person)[2:4],]
  expect_equal(mfrm_mml_inference_identification(subset)$probability_dimension_bound,3)
  # A rank-deficient observed score subset is not this dimension argument.
  repeated$config$estimability_audit$nonlinear_local_estimability <- list(
    local_first_order_rank_deficient=FALSE)
  expect_false(mfrm_mml_inference_identification(repeated)$excluded)
  repeated$config$estimability_audit$nonlinear_local_estimability$local_first_order_rank_deficient <- TRUE
  expect_match(mfrm_mml_inference_identification(repeated)$detail,"first-order rank deficient")
})
