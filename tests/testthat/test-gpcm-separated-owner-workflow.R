separate_owner_fixture <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      withr::local_seed(926740)

      d <- expand.grid(Person=paste0('P',1:120),Rater=paste0('R',1:3),Criterion=paste0('C',1:2),stringsAsFactors=FALSE)
      r <- match(d$Rater,unique(d$Rater)); c <- match(d$Criterion,unique(d$Criterion)); p <- match(d$Person,unique(d$Person))
      theta <- rnorm(120,sd=1.2)
      steps <- rbind(c(-.8,.8),c(-.45,.45),c(-1,1)); a <- c(.8,1.25)
      probs <- t(vapply(seq_len(nrow(d)),function(i) {
        eta <- theta[p[i]]-c(-.25,0,.25)[r[i]]-c(-.15,.15)[c[i]]
        z <- c(0,cumsum(a[c[i]]*(eta-steps[r[i],])))
        w <- exp(z-max(z)); w/sum(w)
      },numeric(3)))
      d$Score <- vapply(seq_len(nrow(d)),function(i) sample.int(3,1,prob=probs[i,])-1L,integer(1))
      f <- fit_mfrm(d,'Person',c('Rater','Criterion'),'Score',model='GPCM',method='MML',
        step_facet='Rater',slope_facet='Criterion',
        quad_points=31,maxit=400,reltol=1e-10)

      cached <<- list(data=d,fit=f)
    }
    cached
  }
})

test_that("MML keeps distinct owners through inference, curves and saved output", {
  x <- separate_owner_fixture(); f <- x$fit
  expect_identical(f$config$step_facet, "Rater")
  expect_identical(f$config$slope_facet, "Criterion")
  expect_setequal(f$steps$StepFacet, paste0("R",1:3))
  expect_identical(f$slopes$SlopeFacet, paste0("C",1:2))
  expect_equal(prod(f$slopes$Estimate),1,tolerance=1e-12)
  ci <- confint(f, parm = "slopes")
  expect_true(all(attr(ci,"diagnostics")$CIEligible))
  cv <- mfrm_gpcm_inference(f,method="sandwich")
  expect_true(cv$check$eligible)
  expect_equal(cv$covariance,cv$information$cov %*% crossprod(cv$person_scores) %*%
    cv$information$cov,tolerance=1e-9,ignore_attr=TRUE)
  grid <- expand.grid(Theta=c(-1,0,1),Rater=paste0("R",1:3),Criterion=paste0("C",1:2))
  curves <- mfrm_curve_intervals(f,grid)
  info <- mfrm_curve_intervals(f,grid,type="information")
  expect_identical(curves$contexts$ObservedContext, rep(TRUE, nrow(grid)))
  expect_identical(curves$contexts$InputRow, seq_len(nrow(grid)))
  expect_identical(mfrm_gpcm_inference_tables(curves)$contexts, curves$contexts)
  # Retaining known levels while omitting one design crossing changes the
  # context label, without changing the fitted coefficients used to predict it.
  prep <- f$prep
  prep$data <- prep$data[!(prep$data$Rater == "R1" & prep$data$Criterion == "C1"), ]
  query <- mfrm_gpcm_response_evaluator(prep, f$config, grid)
  expect_identical(query$contexts$ObservedContext,
    !(grid$Rater == "R1" & grid$Criterion == "C1"))
  expect_equal(as.vector(t(query$evaluate(f$opt$par)$probabilities)),
    curves$table$Estimate, tolerance = 1e-12)
  marked <- curves; marked$contexts <- query$contexts
  expect_output(print(marked), "Unobserved combinations")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_match(paste(plot_data(plot(marked, draw = FALSE), "notes")$Text, collapse = "\n"), "not observed together")
    expect_null(plot(marked, caption = NULL, draw = FALSE)$labels$caption)
    unavailable <- marked; unavailable$table$CIEligible[1] <- FALSE
    caption <- paste(plot_data(plot(unavailable, draw = FALSE), "notes")$Text, collapse = "\n")
    expect_match(caption, "Intervals unavailable")
    expect_match(caption, "not observed together")
  }
  pars <- expand_params(f$opt$par,build_param_sizes(f$config),f$config)
  # Independent category recursion, with deliberately different role sizes.
  expected <- t(vapply(seq_len(nrow(grid)),function(i) {
    r <- match(grid$Rater[i],f$prep$levels$Rater)
    c <- match(grid$Criterion[i],f$prep$levels$Criterion)
    eta <- grid$Theta[i]-pars$facets$Rater[r]-pars$facets$Criterion[c]
    z <- c(0,cumsum(pars$slopes[c]*(eta-pars$steps_mat[r,])))
    w <- exp(z-max(z)); w/sum(w)
  },numeric(3)))
  expect_equal(curves$table$Estimate,as.vector(t(expected)),tolerance=1e-12)
  expect_true(all(curves$table$CIEligible))
  variance <- drop(expected %*% (0:2)^2)-drop(expected %*% (0:2))^2
  expect_equal(info$table$Estimate,variance*pars$slopes[match(grid$Criterion,f$prep$levels$Criterion)]^2,
    tolerance=1e-10,ignore_attr=TRUE)
  total <- compute_information(f,theta_range=c(-1,1),theta_points=3)
  expect_equal(total$tif$Information,as.numeric(rowsum(info$table$Estimate,grid$Theta))*120,
    tolerance=1e-9)
  result <- mfrm_results(f,intervals=list(slopes=ci,curves=curves),
    include=c("fit","plots"),compute="never")
  path <- tempfile(fileext=".rds"); withr::defer(unlink(path))
  saveRDS(result,path); reopened <- readRDS(path)
  expect_identical(reopened$gpcm_inference,result$gpcm_inference)
  saveRDS(marked, path)
  expect_identical(readRDS(path)$contexts, query$contexts)
  expect_identical(reopened$fit$config$slope_facet,"Criterion")
  report <- mfrm_report(reopened)
  expect_identical(report$tables$gpcm_slopes_intervals,result$tables$gpcm_slopes_intervals)
  if (requireNamespace("ggplot2",quietly=TRUE)) {
    plot <- as_ggplot(reopened,type="gpcm_curves",title=NULL,subtitle=NULL)
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_equal(plot_data(curves,"table")$Estimate,curves$table$Estimate)
  }
  changed <- f; changed$config$step_facet <- "Criterion"
  expect_false(mfrm_gpcm_slope_inference_check(changed,cv$information)$eligible)
  expect_error(mfrm_results(changed,intervals=ci,compute="never"),"must match")
})

test_that("resampling keeps both owners and specialist routes refuse unsupported structures", {
  x <- separate_owner_fixture(); f <- x$fit
  dat <- mfrm_gpcm_bootstrap_generate(f,926741)
  args <- mfrmr_gqs_refit_arguments(f,dat,31L)
  expect_identical(args$step_facet,"Rater")
  expect_identical(args$slope_facet,"Criterion")
  expect_equal(dat[c("Person","Rater","Criterion")],x$data[c("Person","Rater","Criterion")],ignore_attr=TRUE)
  expect_true(all(dat$Score %in% 0:2))
  expect_error(extract_mfrm_sim_spec(f),"currently requires")
  expect_error(build_mfrm_sim_spec(model="GPCM",step_facet="Rater",slope_facet="Criterion"),
    "slope_facet == step_facet",fixed=TRUE)
  expect_error(simulation_resolve_fit_slope_facet("GPCM","Criterion","Rater"),"Simulation workflows")
  expect_error(fit_mfrm(x$data,"Person",c("Rater","Criterion"),"Score",model="GPCM",method="JML",
    step_facet="Rater",slope_facet="Criterion"),"require MML")
})

test_that("joint-owner reference curves retain every step/slope pair without borrowing fit flags", {
  f <- separate_owner_fixture()$fit
  spec <- build_step_curve_spec(f)
  expect_length(spec$groups,6)
  curves <- build_curve_tables(spec,c(-1,0,1))
  pars <- expand_params(f$opt$par,build_param_sizes(f$config),f$config)
  for (r in 1:3) for (c in 1:2) {
    name <- paste0("Rater = R",r,"; Criterion = C",c)
    expected <- vapply(c(-1,0,1),function(theta) {
      z <- c(0,cumsum(pars$slopes[c]*(theta-pars$steps_mat[r,])))
      p <- exp(z-max(z)); sum((0:2)*p/sum(p))
    },numeric(1))
    expect_equal(curves$expected$ExpectedScore[curves$expected$CurveGroup==name],expected,
      tolerance=1e-12)
  }
  pathway <- build_pathway_map_data(f,theta_points=3)
  expect_false(any(pathway$curve_fit_status$MatchedFitRow))
  expect_true(all(grepl("joint slope/step profile",pathway$curve_fit_status$ReviewReason)))
  expect_gt(nrow(pathway$fit_measures),0)
  expect_setequal(spec$step_points$StepLevel,paste0("R",1:3))
  for (type in c("fit_pathway","pathway","ccc","wright")) {
    expect_no_error(suppressWarnings(plot(f,type=type,draw=FALSE)))
  }
})

test_that("perfectly confounded owners are rejected before estimation", {
  d <- separate_owner_fixture()$data
  d <- d[(d$Rater=="R1" & d$Criterion=="C1") | (d$Rater=="R2" & d$Criterion=="C2"),]
  expect_error(fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",model="GPCM",method="MML",
    step_facet="Rater",slope_facet="Criterion"),"structurally unidentified")
})
