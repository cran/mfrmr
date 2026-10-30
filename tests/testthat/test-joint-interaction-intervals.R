joint_cell_fixture <- local({
  saved <- list()
  function(model = "PCM") {
    if (is.null(saved[[model]])) {
      d <- load_mfrmr_data("example_core")
      f <- suppressWarnings(fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",
        model=model,step_facet=if(model=="PCM")"Criterion"else NULL,
        noncenter_facet="Criterion",facet_interactions="Rater:Criterion",
        quad_points=61L,maxit=400L,reltol=1e-11))
      ci <- suppressWarnings(mfrm_facet_intervals(f,interaction="Rater:Criterion"))
      saved[[model]] <<- list(fit=f,ci=ci)
    }
    saved[[model]]
  }
})

test_that("joint cells expand full nuisance-adjusted information with zero margins", {
  for (model in c("RSM","PCM")) {
    x<-joint_cell_fixture(model);ci<-x$ci;f<-x$fit
    expect_true(all(ci$checks$Passed))
    expect_true(all(ci$table$CIEligible))
    expect_false(any(ci$table$CoverageQualified))
    expect_equal(ci$table$Estimate,f$interactions$effects$Estimate)
    spec<-f$config$interaction_specs[[1]]
    # Independent Kronecker expansion in the stored column-major cell order.
    A<-rbind(diag(spec$n_a-1L),-1);B<-rbind(diag(spec$n_b-1L),-1)
    G<-kronecker(B,A);at<-build_param_slices(build_param_sizes(f$config))$interactions
    expect_equal(ci$target_jacobian[,at],G,ignore_attr=TRUE)
    V<-ci$parameter_covariance
    expected<-G %*% V[at,at] %*% t(G)
    expect_equal(vcov(ci),expected,tolerance=1e-11,ignore_attr=TRUE)
    # Inverting the interaction-only Hessian would omit nuisance uncertainty.
    conditional<-G %*% solve(solve(V)[at,at]) %*% t(G)
    expect_gt(max(abs(diag(expected)-diag(conditional))),1e-5)
    estimate<-matrix(ci$table$Estimate,spec$n_a,spec$n_b)
    expect_equal(rowSums(estimate),rep(0,spec$n_a),tolerance=1e-12)
    expect_equal(colSums(estimate),rep(0,spec$n_b),tolerance=1e-12)
    expect_identical(qr(vcov(ci))$rank,as.integer(spec$n_params))
    expect_identical(ci$interaction_cells$N,as.data.frame(f$interactions$effects)$N)
    expect_equal(ci$table$Lower,ci$table$Estimate-qnorm(.975)*ci$table$SE)
  }
})

test_that("named cell contrasts include correlations and fixed sums stay fixed", {
  x<-joint_cell_fixture();ci<-x$ci;spec<-x$fit$config$interaction_specs[[1]]
  C<-matrix(0,3,nrow(ci$table),dimnames=list(c("difference","double","column sum"),ci$table$Target))
  C[1,1:2]<-c(1,-1);C[2,]<-2*C[1,];C[3,seq_len(spec$n_a)]<-1
  out<-suppressWarnings(mfrm_facet_intervals(x$fit,interaction="Rater:Criterion",
    contrasts=C[,ncol(C):1],level=.9))
  expect_equal(out$table$Estimate,drop(C %*% ci$table$Estimate),ignore_attr=TRUE)
  expect_equal(vcov(out),C %*% vcov(ci) %*% t(C),tolerance=1e-10)
  expect_equal(out$table$SE[2],2*out$table$SE[1])
  expect_equal(out$table$SE[1]^2,vcov(ci)[1,1]+vcov(ci)[2,2]-2*vcov(ci)[1,2])
  expect_identical(out$table$Status[3],"fixed")
  expect_equal(out$table$Estimate[3],0)
  expect_equal(out$table$SE[3],0)
  expect_true(is.na(out$table$Lower[3]))
  expect_false(out$table$CIEligible[3])
  expect_identical(out$interaction_cells,ci$interaction_cells)
  expect_error(mfrm_facet_intervals(x$fit,interaction="Rater:Criterion",contrasts=unname(C)),"every cell Target")
})

test_that("joint-cell source validation does not open unsupported models", {
  x<-joint_cell_fixture();f<-x$fit
  expect_error(mfrm_facet_intervals(f,"Rater",interaction="Rater:Criterion"),"omit facet")
  expect_error(mfrm_facet_intervals(f,interaction="Rater:Criterion",method="sandwich"),"method")
  expect_error(mfrm_facet_intervals(f,interaction="Rater:Criterion",scale="standardized"),"scale")
  expect_error(mfrm_facet_intervals(f,interaction="bad"),"one fitted two-facet")
  for(change in list(list(model="GPCM"),list(method="JML"),list(method="CML"),list(population_spec=list(active=TRUE)),
      list(dummy_facets="Criterion"),list(facet_shrinkage="empirical_bayes"))) {
    bad<-f;bad$config[names(change)]<-change
    expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"))
  }
  bad<-f;bad$config$estimation_control$mml_integration<-"adaptive"
  expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"),"fixed-grid")
  bad<-f;class(bad)<-c("mfrm_comparison_fit",class(bad))
  expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"),"full native")
  bad<-f;bad$config$model<-"GPCM"
  with_mocked_bindings(
    expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"),"full native"),
    mfrm_has_product_slopes=function(...)stop("Unsupported route was entered"),.package="mfrmr")
  bad<-f;bad$config$interaction_specs[[1]]$cell_n[1,1]<-999
  expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"),"coding, constraints and counts")
  bad<-f;bad$interactions$effects$Estimate[1]<-bad$interactions$effects$Estimate[1]+.1
  expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"),"Saved interactions")
  bad<-f;bad$prep$data$Weight[1]<-2
  expect_error(mfrm_facet_intervals(bad,interaction="Rater:Criterion"))
})

test_that("failed joint information retains points and explicit unavailable intervals", {
  x<-joint_cell_fixture()
  local_mocked_bindings(mfrm_native_location_inference=function(...)list(covariance=NULL,
    checks=data.frame(Check="Joint information",Passed=FALSE,Detail="Injected singular information"),
    check=list(eligible=FALSE,review="Injected singular information")),.package="mfrmr")
  out<-suppressWarnings(mfrm_facet_intervals(x$fit,interaction="Rater:Criterion"))
  expect_identical(out$table$Estimate,x$ci$table$Estimate)
  expect_true(all(out$table$Status=="numerical_review_failed"))
  expect_true(all(is.na(out$table$SE)&is.na(out$table$Lower)))
  expect_false(any(out$table$CIEligible))
  expect_true(all(grepl("Injected singular",out$table$InferenceReview)))
  p<-plot(out,draw=FALSE)
  expect_identical(plot_data(p)$table,out$table)
})

test_that("saved joint cells reach plots and reports without numerical recalculation", {
  x<-joint_cell_fixture();ci<-unserialize(serialize(x$ci,NULL))
  local_mocked_bindings(fit_mfrm=function(...)stop("No refit"),
    compute_mml_parameter_covariance=function(...)stop("No covariance"),.package="mfrmr")
  p<-plot(ci,draw=FALSE)
  expect_match(plot_data(p)$subtitle,"Joint interaction effects",fixed=TRUE)
  expect_false(grepl("Weak information",plot_data(p)$subtitle,fixed=TRUE))
  expect_identical(plot_data(p)$interaction_cells,ci$interaction_cells)
  expect_identical(plot_data(p)$target_labels,paste(ci$table$FacetA_Level,ci$table$FacetB_Level,sep=" | "))
  expect_s3_class(as_ggplot(p),"ggplot")
  expect_equal(apa_table(ci,which="interaction_cells",digits=12)$table,ci$interaction_cells,tolerance=1e-10)
  res<-mfrm_results(x$fit,intervals=list(cells=ci),include="fit",compute="never")
  expect_identical(res$facet_intervals$cells$table,ci$table)
  expect_identical(res$tables$facet_cells_interaction_cells,ci$interaction_cells)
  bad<-ci;bad$fit$interactions$effects$Estimate[1]<-bad$fit$interactions$effects$Estimate[1]+.1
  expect_error(mfrm_results(x$fit,intervals=list(cells=bad),include="fit",compute="never"),"match")
})
