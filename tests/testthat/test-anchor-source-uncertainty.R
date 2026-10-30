anchor_source_fixture <- local({
  saved <- NULL
  function() {
    if (!is.null(saved)) return(saved)
    data <- load_mfrmr_data("example_core")
    fit <- suppressWarnings(fit_mfrm(data, "Person", c("Rater", "Criterion"), "Score",
      model = "PCM", step_facet = "Criterion", noncenter_facet = "Rater",
      anchors = data.frame(Facet = "Rater", Level = c("R01", "R02"), Anchor = c(.3, -.2)),
      population_formula = ~1, person_data = unique(data["Person"]), population_sd = 1.4,
      quad_points = 121L, maxit = 400L, reltol = 1e-11, anchor_policy = "silent"))
    V <- matrix(c(.04, .015, .015, .09), 2L,
      dimnames = list(c("R01", "R02"), c("R01", "R02")))
    conditional <- suppressWarnings(mfrm_facet_intervals(fit, "Rater"))
    selected <- suppressWarnings(mfrm_facet_intervals(fit, "Rater",
      anchor_covariance = V, independent_anchors = TRUE))
    saved <<- list(fit = fit, covariance = V, conditional = conditional, selected = selected)
    saved
  }
})

test_that("source covariance includes direct-anchor uncertainty and cross-covariance", {
  x <- anchor_source_fixture(); ci <- x$selected; a <- ci$anchor_uncertainty
  expect_true(a$eligible)
  expect_true(all(a$checks$Passed))
  expect_equal(ci$table$Estimate, x$conditional$table$Estimate)
  expect_equal(ci$table$ConditionalSE, x$conditional$table$SE)
  expect_equal(ci$model_covariance, x$conditional$covariance)
  expect_equal(ci$table$ModelLower, x$conditional$table$Lower)
  anchors <- match(a$anchors$Level, ci$table$Target)
  expect_equal(vcov(ci)[anchors, anchors], x$covariance, tolerance = 1e-12)
  expect_true(all(ci$table$CIEligible))
  expect_true(all(ci$table$SE >= ci$table$ConditionalSE))
  expect_true(all(ci$table$ConditionalStatus[anchors] == "fixed"))
  expect_false(any(ci$table$CoverageQualified))
  V <- x$conditional$parameter_covariance
  J <- ci$target_jacobian; D <- J %*% a$sensitivity
  D[anchors, ] <- diag(length(anchors))
  expect_equal(vcov(ci), J %*% V %*% t(J) + D %*% x$covariance %*% t(D),
    tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(ci$parameter_covariance, V + a$sensitivity %*% x$covariance %*% t(a$sensitivity),
    tolerance = 1e-10, ignore_attr = TRUE)
  expect_false(any(grepl("source calibration is not propagated", ci$cautions, fixed = TRUE)))
  other <- suppressWarnings(mfrm_facet_intervals(x$fit, "Criterion",
    anchor_covariance = x$covariance, independent_anchors = TRUE))
  expect_equal(other$anchor_uncertainty$target_sensitivity,
    other$target_jacobian %*% other$anchor_uncertainty$sensitivity)
  expect_lt(max(abs(colSums(other$anchor_uncertainty$target_sensitivity))), 1e-10)
})

test_that("named correlations, singular covariance and nonunit standardization are respected", {
  x <- anchor_source_fixture(); V <- x$covariance
  reversed <- suppressWarnings(mfrm_facet_intervals(x$fit, "Rater",
    anchor_covariance = V[2:1, 2:1], independent_anchors = TRUE))
  expect_equal(vcov(reversed), vcov(x$selected), tolerance = 1e-10)
  C <- matrix(c(1,0,-1,0),1L,dimnames=list("anchor minus new",x$selected$table$Target))
  native <- suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",C,
    anchor_covariance=V,independent_anchors=TRUE))
  standardized <- suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",C,scale="standardized",
    anchor_covariance=V,independent_anchors=TRUE))
  expect_equal(standardized$table$Estimate,native$table$Estimate/1.4)
  expect_equal(standardized$table$SE,native$table$SE/1.4,tolerance=1e-10)
  zero <- suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",
    anchor_covariance=V*0,independent_anchors=TRUE))
  expect_equal(zero$table$SE,x$conditional$table$SE,tolerance=1e-10)
  expect_identical(zero$table$Status,x$conditional$table$Status)
  singular <- V; singular[,]<-.04
  rank1 <- suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",
    anchor_covariance=singular,independent_anchors=TRUE))
  expect_true(rank1$anchor_uncertainty$eligible)
  C <- rbind(Difference = c(1,0,-1,0), Overflow = c(1e200,0,-1e200,0))
  colnames(C) <- x$selected$table$Target
  overflow <- suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",C,
    anchor_covariance=V,independent_anchors=TRUE))
  expect_true(overflow$table$CIEligible[1])
  expect_equal(overflow$table$SE[1],native$table$SE)
  expect_false(overflow$table$CIEligible[2])
  expect_true(is.finite(overflow$table$Estimate[2]))
  expect_true(all(is.na(overflow$table[2,c("SE","Lower","Upper","AnchorSourceSE")])))
})

test_that("population reuse applies source covariance once and survives saved replay", {
  x <- anchor_source_fixture()
  direct <- suppressWarnings(mfrm_population_intervals(x$fit,
    anchor_covariance=x$covariance,independent_anchors=TRUE))
  path <- tempfile(fileext=".rds");on.exit(unlink(path),add=TRUE)
  saveRDS(x$selected,path);restored<-readRDS(path)
  local_mocked_bindings(mfrmr_anchor_sensitivity=function(...) stop("No sensitivity recalculation"),
    mfrm_native_location_inference=function(...) stop("No information recalculation"),.package="mfrmr")
  reused <- suppressWarnings(mfrm_population_intervals(restored))
  expect_equal(reused$table[setdiff(names(reused$table), "InferenceCaution")],
    direct$table[setdiff(names(direct$table), "InferenceCaution")],tolerance=1e-10)
  expect_equal(vcov(reused),vcov(direct),tolerance=1e-10)
  expect_identical(reused$table$Status[2],"fixed")
  expect_true(is.na(reused$table$SE[2]))
  C<-matrix(2,1L,dimnames=list("double mean","(Intercept)"))
  contrast<-suppressWarnings(mfrm_population_intervals(restored,contrasts=C))
  expect_equal(contrast$table$SE[2],2*reused$table$SE[1])
  expect_equal(contrast$contrast_covariance,4*vcov(reused),ignore_attr=TRUE)
  bad<-restored;bad$anchor_uncertainty$covariance[1,1]<-.5
  expect_error(mfrm_population_intervals(bad),"unchanged fit, covariance")
  expect_error(mfrm_population_intervals(restored,anchor_covariance=x$covariance),"already contain")
  expect_error(mfrm_population_intervals(restored,independent_anchors=TRUE),"omit both")
})

test_that("unsupported independence, scope and covariance are refused", {
  x<-anchor_source_fixture();f<-x$fit;V<-x$covariance
  expect_error(mfrm_facet_intervals(f,"Rater",anchor_covariance=V),"independent_anchors = TRUE",fixed=TRUE)
  expect_error(mfrm_facet_intervals(f,"Rater",independent_anchors=TRUE),"Supply.*anchor_covariance")
  for (bad in list(V[,1,drop=FALSE],unname(V),V*NA,replace(V,2,.3),
      matrix(c(.01,.2,.2,.01),2L,dimnames=dimnames(V))))
    expect_error(mfrm_facet_intervals(f,"Rater",anchor_covariance=bad,independent_anchors=TRUE),"anchor_covariance")
  bad<-f;bad$config$estimation_control$mml_integration<-"adaptive"
  expect_error(mfrm_facet_intervals(bad,"Rater",anchor_covariance=V,independent_anchors=TRUE),"fixed-grid")
  bad<-f;bad$config$model<-"GPCM"
  expect_error(mfrm_facet_intervals(bad,"Rater",anchor_covariance=V,independent_anchors=TRUE),"RSM/PCM")
})

test_that("failed sensitivity retains conditional evidence and is not retried on reuse", {
  x<-anchor_source_fixture()
  local_mocked_bindings(mfrmr_anchor_sensitivity=function(...) list(sensitivity=NULL,
    checks=data.frame(Check="Injected derivative failure",Passed=FALSE),eligible=FALSE),.package="mfrmr")
  failed<-suppressWarnings(mfrm_facet_intervals(x$fit,"Rater",
    anchor_covariance=x$covariance,independent_anchors=TRUE))
  expect_equal(failed$table$Estimate,x$conditional$table$Estimate)
  expect_equal(failed$table$ModelSE,x$conditional$table$SE)
  expect_equal(failed$table$ConditionalLower,x$conditional$table$Lower)
  expect_true(all(is.na(failed$table$SE)))
  expect_true(all(!failed$table$CIEligible))
  local_mocked_bindings(mfrmr_anchor_sensitivity=function(...)stop("No retry"),.package="mfrmr")
  population<-suppressWarnings(mfrm_population_intervals(failed))
  expect_true(is.finite(population$table$Estimate[1]))
  expect_false(population$table$CIEligible[1])
  expect_true(is.finite(population$table$ConditionalSE[1]))
})

test_that("failed conditional information never becomes a source-qualified interval", {
  x <- anchor_source_fixture()
  local_mocked_bindings(mfrm_native_location_inference = function(...) list(
    covariance = NULL, checks = data.frame(Check = "Quadrature comparison",
      Passed = FALSE, Detail = "Injected finer-grid failure"),
    check = list(eligible = FALSE, review = "Injected finer-grid failure")),
    mfrmr_anchor_sensitivity = function(...) stop("Must not differentiate"), .package = "mfrmr")
  failed <- suppressWarnings(mfrm_facet_intervals(x$fit, "Rater",
    anchor_covariance = x$covariance, independent_anchors = TRUE))
  expect_equal(failed$table$Estimate, x$conditional$table$Estimate)
  expect_false(failed$anchor_uncertainty$eligible)
  expect_true(all(is.na(failed$table$SE)))
  pop <- suppressWarnings(mfrm_population_intervals(failed))
  expect_true(is.finite(pop$table$Estimate[1]))
  expect_true(is.na(pop$table$ConditionalSE[1]))
  expect_false(pop$table$CIEligible[1])
})

test_that("saved reports and plots identify propagated versus conditional uncertainty", {
  x<-anchor_source_fixture()
  local_mocked_bindings(mfrmr_anchor_sensitivity=function(...)stop("No recalculation"),
    compute_mml_parameter_covariance=function(...)stop("No covariance"),.package="mfrmr")
  p<-plot(x$selected,draw=FALSE)
  expect_true(plot_data(p)$comparison)
  expect_identical(plot_data(p)$table,x$selected$table)
  expect_match(plot_data(p)$subtitle,"Information + anchor source",fixed=TRUE)
  expect_false(grepl("conditional on supplied anchors", plot_data(p)$subtitle, fixed=TRUE))
  rows<-ggplot2::ggplot_build(as_ggplot(p))$data
  expect_true(length(rows)>0L)
  population<-suppressWarnings(mfrm_population_intervals(x$selected))
  pp <- plot(population, draw=FALSE)
  expect_match(plot_data(pp)$subtitle,"Independent anchor-source uncertainty included",fixed=TRUE)
  expect_identical(plot_data(pp)$anchor_uncertainty,population$anchor_uncertainty)
  results<-mfrm_results(x$fit,intervals=list(rater=x$selected,population=population),
    include="fit",compute="never")
  expect_true(any(grepl("anchor_source_covariance",names(results$tables))))
  expect_true(any(grepl("anchor_sensitivity_checks",names(results$tables))))
  expect_identical(results$facet_intervals$rater$table,x$selected$table)
})
