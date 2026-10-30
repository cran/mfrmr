test_that("Jeffreys objective agrees with binomial closed forms, including extremes", {
  n <- 8L
  X <- list(Matrix::Matrix(0,n,1,sparse=TRUE),Matrix::Matrix(1,n,1,sparse=TRUE))
  for (s in 0:n) {
    f <- mfrm_fair_jeffreys_objective(X,matrix(0,n,2),c(rep(1,s),rep(0,n-s)))
    p <- (s+.5)/(n+1); beta <- qlogis(p)
    z <- f$evaluate(beta,TRUE)
    expect_equal(z$gradient,0,tolerance=1e-12)
    expect_equal(z$value,-s*log(p)-(n-s)*log1p(-p)-.5*log(n*p*(1-p)),tolerance=1e-12)
    expect_equal(f$gr(.3),numDeriv::grad(f$fn,.3),tolerance=1e-7)
  }
  # Canonical-logit bias reduction does not mean unbiased expected scores.
  for (p in c(.1,.5,.9)) {
    mass <- dbinom(0:n,n,p)
    expect_equal(sum(mass*(0:n)/n),p,tolerance=1e-14)
    expect_equal(sum(mass*((0:n)+.5)/(n+1))-p,(.5-p)/(n+1),tolerance=1e-14)
  }
  singular <- lapply(X,function(x)cbind(x,x))
  f <- mfrm_fair_jeffreys_objective(singular,matrix(0,n,2),rep(0,n))
  expect_error(f$evaluate(c(0,0),TRUE),"positive definite")
})

fair_jeffreys_case <- local({
  saved <- NULL
  function() {
    if (is.null(saved)) {
      d <- load_mfrmr_data("example_core")
      d <- d[d$Person %in% unique(d$Person)[1:12],]
      saved <<- list(data=d,fit=fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",
        method="JML",rating_min=1,rating_max=4,category_policy="preserve",
        optimizer="BFGS",reltol=1e-12,maxit=1000))
    }
    saved
  }
})

test_that("canonical RSM and anchored PCM penalties have independently checked derivatives", {
  x <- fair_jeffreys_case()
  pcm <- fit_mfrm(x$data,"Person",c("Rater","Criterion"),"Score",method="JML",
    model="PCM",step_facet="Criterion",rating_min=1,rating_max=4,category_policy="preserve",
    anchors=data.frame(Facet="Rater",Level="R01",Anchor=.2),
    optimizer="BFGS",maxit=1000,reltol=1e-12)
  for (fit in list(x$fit,pcm)) {
    d <- mfrm_fair_jeffreys_design(fit)
    f <- mfrm_fair_jeffreys_objective(d$category,d$offset,d$y)
    b <- sin(seq_along(fit$opt$par))*.2
    z <- f$evaluate(b,TRUE)
    idx <- build_indices(fit$prep,fit$config$step_facet,fit$config$slope_facet,fit$config$interaction_specs)
    sizes <- build_param_sizes(fit$config)
    expect_equal(z$nll,mfrm_loglik_jml(b,idx,fit$config,sizes),tolerance=1e-10)
    expect_equal(z$ordinary_gradient,mfrm_grad_jml(b,idx,fit$config,sizes),tolerance=1e-9)
    expect_equal(z$gradient,numDeriv::grad(f$fn,b),tolerance=1e-6)
    expect_lt(z$solve_residual,1e-10)
    # Affine free-coordinate changes add a constant to the log determinant.
    T <- diag(seq(.8,1.2,length.out=length(b))); T[1,2] <- .3
    transformed <- mfrm_fair_jeffreys_objective(lapply(d$category,function(X)X%*%T),d$offset,d$y)
    q <- transformed$evaluate(as.numeric(solve(T,b)),TRUE)
    expect_equal(q$value,z$value-as.numeric(determinant(T,logarithm=TRUE)$modulus),tolerance=1e-10)
    expect_equal(q$gradient,as.numeric(t(T)%*%z$gradient),tolerance=1e-9)
    expect_equal(q$probabilities,z$probabilities,tolerance=1e-12)
  }
})

test_that("Jeffreys generation retains finite probabilities without altering the source", {
  x <- fair_jeffreys_case(); d <- x$data
  d$Score[d$Person==unique(d$Person)[1]] <- 1
  fit <- fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",method="JML",
    rating_min=1,rating_max=4,category_policy="preserve",optimizer="BFGS",maxit=1000,reltol=1e-12)
  before <- serialize(fit,NULL); set.seed(74); rng <- .Random.seed
  g <- mfrm_fair_jeffreys_generator(fit,1e-4)
  expect_identical(serialize(fit,NULL),before)
  expect_identical(.Random.seed,rng)
  expect_true(all(g$probabilities>0 & g$probabilities<1))
  expect_equal(unname(rowSums(g$probabilities)),rep(1,nrow(d)),tolerance=1e-12)
  expect_lt(g$max_gradient,1e-4)
  expect_lt(g$solve_residual,1e-8)
  plan <- mfrm_fair_reference_plan(fit,unique(d[c("Rater","Criterion")]))
  persons <- fit$config$person_levels[1:2]
  target <- mfrm_fair_jeffreys_target(fit,g$par,plan,persons)
  expect_true(all(target>1 & target<4))
  expect_equal(mfrm_fair_reference_scores(fit,plan)$Estimate[1],1)
  # Independent ordinary score calculation at the saved finite generator point.
  pa <- expand_params(g$par,build_param_sizes(fit$config),fit$config)
  probabilities <- category_prob_rsm(pa$theta[1]-rep(pa$facets$Rater,4)-rep(pa$facets$Criterion,each=4),
    c(0,cumsum(pa$steps)))
  expect_equal(unname(target[1]),mean(drop(probabilities%*%1:4)),tolerance=1e-10)
})

test_that("explicit Jeffreys refits propagate same-event anchors and retain ordinary points", {
  x <- fair_jeffreys_case(); d <- x$data; persons <- x$fit$config$person_levels[1:2]
  d$Score[d$Person==persons[1]] <- 1
  calibration <- fit_mfrm(d,"Person",c("Rater","Criterion"),"Score",method="JML",
    rating_min=1,rating_max=4,category_policy="preserve",optimizer="BFGS",maxit=1000,reltol=1e-12)
  a <- subset(calibration$facets$others,Facet=="Rater" & Level=="R01")[,c("Facet","Level","Estimate")]
  names(a)[3] <- "Anchor"
  args <- calibration$config$replay_inputs; args$package_version <- NULL; args$data <- d; args$anchors <- a
  scored <- do.call(fit_mfrm,args)
  plan <- unique(d[c("Rater","Criterion")]); plan$Weight <- seq_len(nrow(plan))
  fair <- fair_average_table(scored,reference_plan=plan)
  maps <- list(calibration=seq_len(nrow(d)),scoring=seq_len(nrow(d)))
  L <- matrix(c(1,-1),1,dimnames=list("difference",persons))
  set.seed(41); rng <- .Random.seed
  ci <- confint(fair,parm=persons,contrasts=L,nsim=3,seed=71,
    anchor_source=calibration,response_map=maps,generator_method="jeffreys")
  expect_identical(.Random.seed,rng)
  expect_equal(ci$Estimate,fair$reference_scores$Estimate[1]-fair$reference_scores$Estimate[2])
  expect_false(ci$CIEligible)
  expect_equal(ci$Resolved+ci$Unresolved,3L)
  expect_identical(attr(ci,"response_map"),maps)
  expect_equal(attr(ci,"settings")$shared_responses,nrow(d))
  expect_identical(attr(ci,"settings")$generator_method,"jeffreys")
  expect_true(all(attr(ci,"generating_model")$probabilities>0))
  expect_match(as.data.frame(ci)$InferenceReview,"Jeffreys generating model")
  expect_match(plot(ci,draw=FALSE)$data$notes$Text,"Jeffreys generating model")
  expect_output(print(ci),"Jeffreys generating model")
  # Reconstruct the first replication independently through the public fitter.
  g <- attr(ci,"generating_model")
  p <- g$probabilities
  set.seed(71); y <- pmin(rowSums(runif(nrow(p))>t(apply(p,1,cumsum))),ncol(p)-1)+1
  cal_args <- calibration$config$replay_inputs; cal_args$package_version <- NULL
  cal_args$data <- as.data.frame(calibration$prep$data[c("Person","Rater","Criterion")])
  cal_args$data$Score <- y
  cal <- do.call(fit_mfrm,cal_args)
  anchors <- a; anchors$Anchor <- subset(cal$facets$others,Facet=="Rater" & Level=="R01")$Estimate
  args$data <- cal_args$data; args$anchors <- anchors
  scoring <- do.call(fit_mfrm,args)
  expected <- fair_average_table(scoring,reference_plan=plan)$reference_scores$Estimate[1:2]
  expect_equal(unname(attr(ci,"person_draws")[1,]),expected,tolerance=1e-10)
  expect_equal(unname(attr(ci,"anchor_draws")[1,]),unname(anchors$Anchor),tolerance=1e-10)
  root <- attr(ci,"draws")[,1]-attr(ci,"generating_target")
  expect_equal(ci$CI_Lower,unname(ci$Estimate-quantile(root,.975,type=1)),tolerance=1e-12)
  path <- tempfile(fileext=".rds"); withr::defer(unlink(path))
  saveRDS(ci,path); expect_identical(readRDS(path),ci)
  expect_error(confint(fair,parm=persons,contrasts=L,nsim=2,seed=71,
    anchor_source=calibration,response_map=maps),"finite primary estimates")
  expect_error(confint(fair,method="delta",generator_method="jeffreys"),"only")
  expect_error(confint(fair,method="bounded",generator_method="ordinary"),"only")
  # Failed calibration remains in the planned denominator and the outer tails.
  actual_refit <- mfrm_fair_refit; calls <- 0L
  local_mocked_bindings(mfrm_fair_refit=function(...) {
    calls <<- calls+1L
    if (calls==1L) stop("Controlled calibration failure")
    actual_refit(...)
  },.package="mfrmr")
  failed <- confint(fair,parm=persons,contrasts=L,nsim=2,seed=71,
    anchor_source=calibration,response_map=maps,generator_method="jeffreys")
  expect_equal(failed$Unresolved,1L)
  expect_equal(failed$Planned,2L)
  expect_false(failed$Available)
  expect_match(attr(failed,"trials")$Reason[1],"Controlled calibration failure")
  # A stopped generator is not replaced by an ordinary or clipped generator.
  local_mocked_bindings(mfrm_fair_jeffreys_generator=function(...)stop("Controlled generator failure"),.package="mfrmr")
  expect_error(confint(fair,parm=persons,nsim=2,seed=71,anchor_source=calibration,
    response_map=maps,generator_method="jeffreys"),"Controlled generator failure")
})
