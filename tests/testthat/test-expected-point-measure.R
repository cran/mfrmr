# Native references are literal FACETS 4.5.0 output from fixed-calibration PCM
# cases. Probability moments below are computed independently of mfrmr kernels.
pcm_correlation_moments <- function(z) {
  eta <- z$PersonMeasure-z$TaskMeasure
  logs <- cbind(0,eta-z$Step1,2*eta-z$Step1-z$Step2,
    3*eta-z$Step1-z$Step2-z$Step3)
  prob <- exp(logs-apply(logs,1,max)); prob <- prob/rowSums(prob)
  e <- drop(prob %*% (0:3)); centered <- matrix(0:3,nrow(z),4,byrow=TRUE)-e
  list(eta=eta,e=e,v=rowSums(prob*centered^2),fourth=rowSums(prob*centered^4))
}

test_that("expected correlations and small PCM df match native weighted cases", {
  input <- read.csv(test_path("fixtures","facets-pcm-correlation-input.csv"))
  ref <- read.csv(test_path("fixtures","facets-pcm-correlation-reference.csv"))
  for (i in seq_len(nrow(ref))) {
    rr <- ref[i,]; z <- input[input$Case==rr$Case & input[[rr$Facet]]==rr$Level,]
    mom <- pcm_correlation_moments(z)
    r <- mfrmr:::.mfrm_expected_correlation_moments(mom$eta,mom$e,mom$v,z$Weight)
    expect_equal(r$PTMEAExpectedN,rr$Count)
    if (rr$Count>=2) {
      expect_identical(r$PTMEAExpectedStatus,"available_conditional")
      expect_lt(abs(r$PTMEAExpected-rr$Expected),.00501)
    } else {
      expect_true(is.na(r$PTMEAExpected))
      expect_identical(r$PTMEAExpectedStatus,"insufficient_frequency")
    }
    # The unadjusted ratio is independently the correlation in the conditional
    # joint distribution obtained by selecting a row with probability w/sum(w).
    w <- z$Weight/sum(z$Weight); mx <- sum(w*mom$eta); me <- sum(w*mom$e)
    expected <- sum(w*(mom$eta-mx)*(mom$e-me))/
      sqrt(sum(w*(mom$eta-mx)^2)*sum(w*((mom$e-me)^2+mom$v)))
    expect_equal(r$PTMEAExpectedAsymptotic,expected,tolerance=1e-13)
    obs <- data.frame(Observed=z$Score,Expected=mom$e,Var=mom$v,Weight=z$Weight,
      StdSq=(z$Score-mom$e)^2/mom$v,FourthCentralMoment=mom$fourth)
    fit <- mfrmr:::mfrm_summarize_response_fit(mfrmr:::mfrm_fit_observation_terms(obs),
      whexact=FALSE,fit_df_method="facets",facets_zstd_cap=9)
    expect_lt(abs(fit$DF_Infit-rr$InfitDF),.00501)
    expect_lt(abs(fit$DF_Outfit-rr$OutfitDF),.00501)
  }
})

test_that("frequency convention, invariance and unavailable correlation scopes are explicit", {
  fun <- mfrmr:::.mfrm_expected_correlation_moments
  x <- c(-1,0,1); e <- plogis(x); v <- e*(1-e); w <- c(1,2,1)
  a <- fun(x,e,v,w)
  repeated <- rep(seq_along(x),w)
  b <- fun(x[repeated],e[repeated],v[repeated],rep(1,length(repeated)))
  expect_equal(a,b)
  shift <- fun(x+1e6,e+100,v,w)
  expect_equal(a$PTMEAExpected,shift$PTMEAExpected,tolerance=1e-11)
  reverse <- fun(-x,e,v,w)
  expect_equal(reverse$PTMEAExpected,-a$PTMEAExpected)
  scaled <- fun(x,e,v,2*w)
  expect_equal(scaled$PTMEAExpectedAsymptotic,a$PTMEAExpectedAsymptotic)
  expect_lt(scaled$PTMEAExpected,a$PTMEAExpected)
  expect_equal(fun(c(-1,1),plogis(c(-1,1)),rep(.2,2),c(1,1))$PTMEAExpected,1)
  # Exact two-Bernoulli enumeration: the sample correlation is undefined on
  # equal responses. Its mean conditional on being defined is .60/.68, not 1.
  # This prevents reinterpreting the finite-count reporting rule as E[r].
  two <- fun(c(-1,1),c(.2,.8),c(.16,.16),c(1,1))
  patterns <- expand.grid(a=0:1,b=0:1)
  mass <- dbinom(patterns$a,1,.2)*dbinom(patterns$b,1,.8)
  defined <- patterns$a != patterns$b
  exact <- sum(mass[defined]*sign(patterns$b[defined]-patterns$a[defined]))/sum(mass[defined])
  expect_equal(sum(mass[defined]),.68)
  expect_equal(exact,.60/.68)
  expect_gt(two$PTMEAExpected-exact,.1)
  expect_equal(fun(c(-1,1),c(.8,.2),c(.16,.16),c(1,1))$PTMEAExpected,-1)
  expect_equal(fun(c(x,NA),c(e,NA),c(v,NA),c(w,0)),a)
  bad <- fun(x,e,replace(v,2,NA),w)
  expect_true(is.na(bad$PTMEAExpected)); expect_true(is.na(bad$PTMEAExpectedAsymptotic))
  expect_identical(bad$PTMEAExpectedStatus,"unavailable")
  expect_match(fun(x,e,v,w,c(TRUE,FALSE,TRUE))$PTMEAExpectedReason,"unanchored extreme")
  expect_match(fun(rep(1,3),e,v,w)$PTMEAExpectedReason,"constant")
  expect_true(is.na(fun(x,rep(1,3),rep(0,3),w)$PTMEAExpected))
  expect_equal(fun(x,rep(1,3),rep(.2,3),w)$PTMEAExpected,0)
})

test_that("public diagnostics and fair tables retain expected references and save/reload", {
  z <- read.csv(test_path("fixtures","facets-pcm-correlation-input.csv"))
  z <- z[z$Case=="unit",]
  anchor <- rbind(unique(data.frame(Facet="Person",Level=z$Person,Anchor=z$PersonMeasure)),
    unique(data.frame(Facet="Task",Level=z$Task,Anchor=z$TaskMeasure)))
  expect_warning(fit <- fit_mfrm(z,"Person","Task","Score",method="JML",model="PCM",step_facet="Task",
    anchors=anchor,min_obs_per_element=1,min_obs_per_category=1,maxit=200),"Category support")
  before <- serialize(fit,NULL)
  d <- diagnose_mfrm(fit,residual_pca="none",diagnostic_mode="legacy")
  expect_true(all(d$measures$PTMEAExpectedStatus=="available_conditional"))
  for (i in seq_len(nrow(d$measures))) {
    row <- d$measures[i,]; obs <- d$obs[as.character(d$obs[[row$Facet]])==row$Level,]
    # Independent covariance calculation at this public fit's own calibration.
    w <- obs$Weight/sum(obs$Weight); x <- obs$MeasureSum; e <- obs$Expected
    cv <- sum(w*(x-sum(w*x))*(e-sum(w*e)))
    target <- cv/sqrt(sum(w*(x-sum(w*x))^2)*
      (sum(w*(e-sum(w*e))^2)+(sum(obs$Weight)-2)/sum(obs$Weight)*sum(w*obs$ConditionalVariance)))
    expect_equal(row$PTMEAExpected,target,tolerance=1e-12)
  }
  b <- fair_average_table(fit,d)
  for (f in names(b$raw_by_facet)) {
    tab <- b$raw_by_facet[[f]]
    src <- d$measures[d$measures$Facet==f,]
    expect_equal(tab$PTMEAExpected,src$PTMEAExpected[match(tab$Level,src$Level)])
  }
  saved <- unserialize(serialize(list(fit=fit,d=d,b=b),NULL))
  expect_identical(summary(saved$b),summary(b))
  expect_no_error(plot(saved$fit,type="fair_average",diagnostics=saved$d,draw=FALSE))
  expect_identical(serialize(fit,NULL),before)
  old <- d; old$measures <- old$measures[,!startsWith(names(old$measures),"PTMEAExpected")]
  expect_false("PTMEAExpected" %in% names(fair_average_table(fit,old)$stacked))
  # An invalid row may not disappear from the expected reference silently.
  bad <- d$obs; bad$Expected[1] <- NA_real_
  refs <- mfrmr:::.mfrm_expected_ptmea(fit,bad,c("Person","Task"))
  expect_true(all(is.na(refs$PTMEAExpected[refs$Level=="1"])))
})
