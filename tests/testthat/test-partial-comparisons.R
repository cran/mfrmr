partial_comparison_cases <- function() readRDS(test_path("fixtures", "partial-comparisons.rds"))
partial_comparison_args <- function(x) c(list(data=x$data, fit_scope="comparisons",
  comparison_review=x$review, interaction_policy="silent"), x$review$comparison_review$fit_arguments)

test_that("comparison-only fits preserve unidentified rows and dedicated result methods", {
  for (x in partial_comparison_cases()) {
    f <- x$fit; out <- interaction_effect_table(f)
    expect_s3_class(f, "mfrm_comparison_fit")
    expect_false(inherits(f, "mfrm_fit"))
    expect_identical(out$ConditionalState, c("target_unidentified", "target_identified_only"))
    expect_true(is.na(out$Estimate[1]) && is.finite(out$Estimate[2]))
    expect_identical(out$EstimateStatus, c("unavailable", "descriptive_only"))
    expect_true(all(is.na(out$SE) & is.na(out$Lower) & is.na(out$Upper)))
    expect_false(any(out$FormalInferenceEligible))
    expect_identical(summary(f), out)
    expect_identical(interaction_effect_table(unserialize(serialize(f,NULL))), out)
    expect_true(all(is.na(attr(out,"comparison_review")$table$Estimate)))
    expect_output(print(f), "Descriptive estimates only")
    expect_error(predict(f), "do not supply primary Person estimates")
    expect_error(interaction_effect_table(f, x$review), "omit comparison_review")
    if (f$config$method == "MML") {
      expect_identical(out$State, c("target_unidentified", "not_evaluated"))
      expect_true(all(out$QuadratureState == "stable_at_retained_point"))
      expect_equal(nrow(f$checks$starts), 3L)
      expect_true(all(out$BoundaryState == "finite_mle_not_certified"))
    } else expect_true(all(out$BoundaryState == "no_recession_detected"))
  }
})

test_that("a different likelihood representative preserves only invariant targets", {
  for (x in partial_comparison_cases()) {
    f <- x$fit; cfg <- f$config; sizes <- build_param_sizes(cfg)
    idx <- build_indices(f$prep,cfg$step_facet,interaction_specs=cfg$interaction_specs)
    adjacent <- mfrmr_estimability_adjacent_design(f$prep,idx,cfg,sizes)
    A <- as.matrix(adjacent$design); par <- f$optimization$representative
    map <- adjacent$map$OptimizerIndex
    # QR selects independent columns, independently of the fitter's scaled SVD.
    qrA <- qr(A,tol=1e-10,LAPACK=FALSE)
    kept <- qrA$pivot[seq_len(qrA$rank)]
    other <- par
    other[map] <- 0
    other[map[kept]] <- qr.solve(A[,kept,drop=FALSE],drop(A%*%par[map]))
    expect_equal(drop(A%*%other[map]),drop(A%*%par[map]),tolerance=1e-11)
    C <- f$comparison_review$target_map
    expect_equal(drop(C[2,,drop=FALSE]%*%other[map]),f$comparisons$Estimate[2],tolerance=1e-11)
    # An exact same-law direction moves the explicitly unavailable comparison.
    delta <- other[map]-par[map]
    expect_gt(abs(drop(C[1,,drop=FALSE]%*%delta)),1e-4)
    mml <- cfg$method=="MML"
    ev <- make_mfrm_direct_evaluator(cfg$method,make_param_cache(sizes,cfg,idx,is_mml=mml),
      idx,cfg,sizes,if(mml)gauss_hermite_normal(63) else NULL)
    expect_equal(ev$value(other),ev$value(par),tolerance=1e-10)
    expect_lt(max(abs(ev$gradient(other))),1e-4)
  }
})

test_that("the live partial route estimates the requested contrast without weakening full fits", {
  skip_if_not_installed("lpSolve")
  x <- partial_comparison_cases()[["RSM-JML"]]
  a <- partial_comparison_args(x)
  f <- do.call(fit_mfrm,a)
  expect_equal(f$comparisons$Estimate,x$fit$comparisons$Estimate,tolerance=1e-6)
  expect_lt(f$checks$numerical[["GradientMaxAbs"]],1e-4)
  expect_gt(f$representation$nullity,0)
  a$fit_scope <- "full"; a$comparison_review <- NULL
  expect_error(do.call(fit_mfrm,a),class="mfrmr_estimability_error")
  a$comparison_review <- x$review
  expect_error(do.call(fit_mfrm,a),"Use fit_scope")
  a <- partial_comparison_args(x)
  a$data <- a$data[rev(seq_len(nrow(a$data))),]
  a$comparison_review$comparison_review$specification$groups <- c("B","A")
  # Saved rank/target fields are not trusted: the declarative request is rebuilt.
  a$comparison_review$comparison_review$target_map[,] <- 0
  reverse <- do.call(fit_mfrm,a)
  expect_equal(reverse$comparisons$Estimate,-f$comparisons$Estimate,tolerance=1e-5)
})

test_that("source, model and unsupported controls cannot silently change the likelihood", {
  x <- partial_comparison_cases()[["RSM-JML"]]
  a <- partial_comparison_args(x)
  a$data$Score[1] <- if(a$data$Score[1]==1)2 else 1
  expect_error(do.call(fit_mfrm,a),"retained responses")
  a <- partial_comparison_args(x); a$dummy_facets <- NULL
  expect_error(do.call(fit_mfrm,a),"model, constraints")
  a <- partial_comparison_args(x); a$comparison_review <- NULL
  expect_error(do.call(fit_mfrm,a),"Supply the comparison_review")
  for (option in c("anchors","checkpoint","attach_diagnostics","jml_correction_order","mml_engine")) {
    a <- partial_comparison_args(x)
    a[[option]] <- switch(option,anchors=data.frame(),checkpoint="unused",attach_diagnostics=TRUE,
      jml_correction_order=1,mml_engine="em")
    expect_error(do.call(fit_mfrm,a),"Comparison fitting does not accept")
  }
})

test_that("recession and numerical failures retain missing estimates with reasons", {
  skip_if_not_installed("lpSolve")
  x <- partial_comparison_cases()[["RSM-JML"]]
  a <- partial_comparison_args(x); a$maxit <- 1L
  failed <- do.call(fit_mfrm,a)
  expect_true(all(is.na(failed$comparisons$Estimate)))
  expect_true(all(failed$comparisons$NumericalState=="review"))
  expect_match(failed$comparisons$EstimateReason[2],"stationary point")
  # One extreme Person causes a strict recession in the unchanged JML likelihood.
  d <- x$data; d$Score[d$Person==d$Person[1]] <- 4L
  r <- do.call(describe_mfrm_data,c(list(data=d,person="Person",facets=c("Rater","Criterion","Group"),
    score="Score",rating_min=1,rating_max=4,category_policy="preserve",include_agreement=FALSE,
    comparisons=x$review$comparison_review$specification)))
  a <- partial_comparison_args(x); a$data<-d; a$comparison_review<-r
  boundary <- do.call(fit_mfrm,a)
  expect_true(all(is.na(boundary$comparisons$Estimate)))
  expect_true(all(boundary$comparisons$BoundaryState=="recession_detected"))
  expect_true(boundary$checks$boundary$certified)
  expect_null(boundary$optimization)
  m <- partial_comparison_cases()[["RSM-MML"]]
  a <- partial_comparison_args(m); a$quad_points <- 7L
  coarse <- do.call(fit_mfrm,a)
  expect_true(all(is.na(coarse$comparisons$Estimate)))
  expect_true(all(coarse$comparisons$NumericalState=="review"))
  expect_true(all(coarse$comparisons$QuadratureState=="review"))
})
