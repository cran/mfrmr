population_testlet_fixture <- function() {
  counts <- matrix(c(8,3,3,3,4,3,3,3,8),3,3)
  pattern <- list(c(0,0),c(0,1),c(1,1))
  scores <- do.call(rbind,lapply(1:3,function(a) do.call(rbind,lapply(1:3,function(b)
    matrix(rep(c(pattern[[a]],pattern[[b]]),counts[a,b]),ncol=4,byrow=TRUE)))))
  data.frame(Person=rep(seq_len(nrow(scores)),each=4),Block=rep(c('a','a','b','b'),nrow(scores)),
    Score=as.vector(t(scores)))
}

test_that("estimated ability variance flows through scoring, reporting and saved output", {
  d <- population_testlet_fixture()
  fit <- fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,quad_points=61)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$InformationPositive)
  expect_gt(fit$calibration$person_variance,0)
  expect_gt(fit$calibration$variance,0)
  expect_null(fit$settings$fixed_person_sd)
  expect_equal(summary(fit)$person_variance,fit$calibration$person_sd^2)
  # A new Person with no observed ratings receives this population, not N(0,1).
  new <- d[d$Person==1,]; new$Person <- "new"; new$Score <- NA_real_
  prior <- predict(fit,new,missing="omit")
  expect_identical(prior$table$Status,"prior_only")
  expect_equal(prior$table$ConditionalSD,fit$calibration$person_sd)
  expect_equal(prior$table$Lower,qnorm(.025)*fit$calibration$person_sd)
  scores <- predict(fit,d[d$Person==1,])
  expect_identical(scores$table$Status,"available_conditional")
  expect_equal(as.numeric(scores$table[c("Estimate","ConditionalSD")]),
    as.numeric(fit$moments[match("1",fit$input$persons),]),tolerance=1e-6)
  # Independent continuous 2D integration for an all-zero response pattern.
  sd <- fit$calibration$person_sd; local_sd <- sqrt(fit$calibration$variance)
  step <- fit$calibration$steps
  density <- function(theta) vapply(theta,function(t) {
    block <- integrate(function(g) dnorm(g)*plogis(step-t-local_sd*g)^2,
      -Inf,Inf,rel.tol=1e-9)$value
    dnorm(t,sd=sd)*block^2
  },numeric(1))
  mass <- integrate(density,-Inf,Inf,rel.tol=1e-8)$value
  expect_equal(integrate(density,-Inf,scores$table$Lower,rel.tol=1e-8)$value/mass,.025,tolerance=1e-6)
  expect_equal(integrate(density,scores$table$Upper,Inf,rel.tol=1e-8)$value/mass,.025,tolerance=1e-6)
  path <- tempfile(); saveRDS(fit,path)
  expect_identical(predict(readRDS(path),new,missing="omit"),prior)
  res <- mfrm_results(fit,predictions=scores)
  row <- subset(res$tables$variance,Effect=="Person ability")
  expect_equal(row$Variance,fit$calibration$person_variance)
  expect_true(row$Estimated)
  expect_false(row$EstimatedBoundary)
  changed <- fit; changed$calibration$person_variance <- 1
  expect_error(mfrm_results(changed,predictions=scores),"matching source")
})

test_that("known nonunit populations and old fixed-normal fits retain their meaning", {
  d <- population_testlet_fixture()
  fit <- fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,
    testlet_variance=0,person_sd=.7,quad_points=61)
  expect_equal(fit$calibration$person_variance,.49)
  expect_false(fit$checks$EstimatedPersonVarianceBoundary)
  expect_false(subset(mfrm_results(fit)$tables$variance,Effect=="Person ability")$Estimated)
  # The old parameter layout is interpreted as N(0,1), never the new default.
  old <- fit; old$parameters <- head(old$parameters,-1)
  old$calibration$person_variance <- old$calibration$person_sd <- NULL
  old$settings$fixed_person_sd <- old$settings$person_variance_max <- NULL
  old$settings$person_distribution <- "N(0,1)"
  old$checks$EstimatedPersonVarianceBoundary <- old$checks$PersonVarianceScore <- NULL
  reference <- old; reference$parameters <- c(old$parameters,1)
  reference$calibration$person_variance <- reference$calibration$person_sd <- 1
  reference$settings$fixed_person_sd <- 1
  new <- d[d$Person==1,]
  expect_equal(predict(old,new)$table,predict(reference,new)$table)
  new$Score <- NA_real_
  expect_equal(predict(old,new,missing="omit")$table$ConditionalSD,1)
  expect_false(subset(mfrm_results(old)$tables$variance,Effect=="Person ability")$Estimated)
})

test_that("both variance boundaries and invalid or capped populations remain explicit", {
  d <- expand.grid(Item=1:2,Block=c("a","b"),Person=1:3)
  d$Score <- rep(c(0,1),6)
  for (sd in list(0,-1,Inf,NaN,1e-300,1e200,c(1,2)))
    expect_error(fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,person_sd=sd),"person_sd")
  expect_error(fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,person_variance_max=0),"variance bounds")
  fit <- fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,quad_points=31)
  expect_true(fit$checks$NumericalReady)
  expect_true(fit$checks$EstimatedPersonVarianceBoundary)
  expect_equal(fit$calibration$person_variance,0)
  expect_lte(fit$checks$PersonVarianceScore,0)
  expect_true(all(is.na(fit$calibration_table$SE)))
  scores <- predict(fit)
  expect_true(all(scores$table$Status=="unavailable"))
  expect_true(all(grepl("ability variance is zero",scores$table$Reason)))
  expect_true(all(is.na(scores$table$Estimate)))
  expect_true(subset(mfrm_results(fit,predictions=scores)$tables$variance,Effect=="Person ability")$EstimatedBoundary)
  # A search cap is not evidence of a fitted population variance.
  cap <- suppressWarnings(fit_mfrm_testlet(population_testlet_fixture(),"Person","Score","Block",
    score_levels=0:1,person_variance_max=.01,quad_points=31))
  expect_true(cap$checks$SearchBoundary)
  expect_false(cap$checks$NumericalReady)
  expect_error(predict(cap),"checks")
})

test_that("zero local variance reduces to an independently fitted normal Rasch population", {
  d <- population_testlet_fixture()
  fit <- fit_mfrm_testlet(d,"Person","Score","Block",score_levels=0:1,
    testlet_variance=0,quad_points=61)
  expect_true(fit$checks$NumericalReady)
  totals <- table(factor(rowsum(d$Score,d$Person),levels=0:4))
  # Ordered-response probabilities: do not insert binomial coefficients.
  nll <- function(p) -sum(totals*vapply(0:4,function(s) log(integrate(function(z) {
    eta <- sqrt(p[2])*z-p[1]
    exp(s*plogis(eta,log.p=TRUE)+(4-s)*plogis(-eta,log.p=TRUE))*dnorm(z)
  },-Inf,Inf,rel.tol=1e-10)$value),numeric(1)))
  reference <- nlminb(c(0,1),nll,lower=c(-20,0),upper=c(20,16),
    control=list(rel.tol=1e-11,x.tol=1e-9))
  expect_equal(reference$convergence,0)
  expect_equal(fit$parameters[c(1,3)],reference$par,tolerance=2e-5)
  expect_equal(fit$loglik,-reference$objective,tolerance=1e-8)
  reference_covariance <- solve(optimHess(reference$par,nll))
  expect_equal(fit$covariance[c(1,3),c(1,3)],reference_covariance,tolerance=1e-4)
})
