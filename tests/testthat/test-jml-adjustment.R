test_that("general JML adjustments preserve the frozen exact response calculation", {
  saved <- readRDS(test_path("fixtures","jml-adjustment-reference.rds"))
  for (x in saved) {
    p <- mfrm_jml_adjustment_problem(x$data,"Person",c("Rater","Criterion"),"Score",x$owner,2L)
    for (j in seq_along(x$orders)) {
      z <- p$evaluate(x$beta,x$orders[j])
      expect_equal(z$value,x$scores[[j]],tolerance=1e-8,ignore_attr=TRUE)
      expect_lt(z$root_residual,1e-9)
    }
    expect_equal(p$evaluate(x$beta,0L)$q,unname(x$q),tolerance=1e-9)
    theta <- p$evaluate(x$beta,0L)$theta
    expect_true(any(theta == -Inf) && any(theta == Inf))
    expect_true(all(p$evaluate(x$beta,4L)$value[!is.finite(theta),] == 0))
  }
})

# Literal three-judge, two-criterion, four-category model, independently
# enumerating ordered response sequences rather than owner-total convolutions.
jml_literal_reference <- function(owner, cells, beta) {
  ng <- if (owner == "Rater") 3L else 2L
  y <- as.matrix(expand.grid(rep(list(0:3),nrow(cells))))
  n <- nrow(y); extreme <- rowSums(y) %in% c(0,3*ncol(y))
  evaluate <- function(b,theta=NULL) {
    rater <- c(b[1:2],-sum(b[1:2])); criterion <- c(b[3],-b[3])
    step <- matrix(b[3+seq_len(ng*2L)],ng,2L,byrow=TRUE)
    step <- cbind(step,-rowSums(step))
    a <- exp(c(tail(b,ng-1L),-sum(tail(b,ng-1L))))[cells[[owner]]]
    logp <- function(t) {
      out <- t(vapply(seq_len(nrow(cells)),function(j) {
        z <- a[j]*((t-rater[cells$Rater[j]]-criterion[cells$Criterion[j]])*(0:3)-
          c(0,cumsum(step[cells[[owner]][j],])))
        z-max(z)-log(sum(exp(z-max(z))))
      },numeric(4)))
      out
    }
    profiled <- is.null(theta)
    if (profiled) theta <- vapply(seq_len(n),function(i) {
      if (extreme[i]) return(0)
      score <- function(t) sum(a*(drop(exp(logp(t))%*%(0:3))-y[i,]))
      uniroot(score,c(-100,100),tol=1e-12)$root
    },numeric(1))
    q <- vapply(seq_len(n),function(i) -sum(logp(theta[i])[cbind(seq_len(ncol(y)),y[i,]+1L)]),numeric(1))
    if(profiled) q[extreme] <- 0
    list(q=q,theta=theta,logp=logp)
  }
  base <- evaluate(beta)
  U <- vapply(seq_along(beta),function(j) {
    lo <- hi <- beta; lo[j] <- lo[j]-1e-5; hi[j] <- hi[j]+1e-5
    (evaluate(hi)$q-evaluate(lo)$q)/2e-5
  },numeric(n))
  mass <- matrix(0,n,n)
  for(i in which(!extreme)) {
    lp <- base$logp(base$theta[i])
    mass[,i] <- exp(rowSums(vapply(seq_len(ncol(y)),function(j)lp[j,y[,j]+1L],numeric(n))))
  }
  mass[cbind(which(extreme),which(extreme))] <- 1
  stopifnot(max(abs(colSums(mass)-1))<1e-12)
  scores <- list(U)
  for(k in 1:4) scores[[k+1L]] <- scores[[k]]-crossprod(mass,scores[[k]])
  data <- do.call(rbind,lapply(seq_len(n),function(i)
    data.frame(Person=sprintf("p%03d",i),cells,Score=y[i,])))
  list(data=data,scores=scores,q=base$q,beta=beta)
}

test_that("different level/category counts agree with literal full enumeration", {
  cells <- data.frame(Rater=1:3,Criterion=c(1,1,2))
  for(owner in c("Criterion","Rater")) {
    ng <- if(owner=="Rater") 3L else 2L
    beta <- c(.2,-.1,.3,rep(c(-.6,.2),ng),if(ng==2L) .18 else c(.18,-.12))
    x <- jml_literal_reference(owner,cells,beta)
    p <- mfrm_jml_adjustment_problem(x$data,"Person",c("Rater","Criterion"),"Score",owner,3L)
    expect_equal(nrow(p$parameters),length(beta))
    for(k in c(0L,1L,2L,4L)) expect_equal(p$evaluate(beta,k)$value,x$scores[[k+1L]],tolerance=1e-6)
    expect_equal(p$evaluate(beta,0L)$q,x$q,tolerance=1e-9)
    # Reorder judge levels, changing free coordinates but not the response model.
    d <- x$data; d$Rater <- c("z","m","a")[d$Rater]
    if(owner=="Criterion") {
      renamed <- mfrm_jml_adjustment_problem(d,"Person",c("Rater","Criterion"),"Score",owner,3L)
      b <- beta; b[1:2] <- c(-sum(beta[1:2]),beta[2])
      J <- diag(length(beta)); J[1:2,1:2] <- rbind(c(-1,-1),c(0,1))
      expect_equal(renamed$evaluate(b,2L)$value,p$evaluate(beta,2L)$value%*%J,tolerance=1e-9)
    }
  }
})

test_that("sample covariance retains the full equation derivative and actual Persons", {
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  p <- do.call(mfrm_jml_adjustment_problem,x[c("data","person","facets","score","owner","rating_max")])
  expect_equal(unname(p$total_states),c(91L,91L))
  expect_equal(as.integer(table(p$roster)),c(200L,200L))
  for(k in c(2L,4L)) {
    b <- x$results[[as.character(k)]]$beta
    expect_lt(max(abs(p$mean_score(b,k))),1e-7)
    V <- mfrm_jml_adjustment_covariance(p,b,k)
    expect_equal(V$vcov,x$results[[as.character(k)]]$covariance$vcov,tolerance=1e-10)
    expect_gt(max(abs(V$jacobian-t(V$jacobian))),1e-5)
    expect_equal(V$vcov,crossprod(V$influence)/length(p$persons)^2,tolerance=1e-12)
    expect_gt(min(eigen(V$vcov,symmetric=TRUE)$values),0)
    random <- mfrm_jml_adjustment_covariance(p,b,k,sampling="random_rosters")
    global <- colMeans(V$scores)
    between <- Reduce(`+`,lapply(split(seq_along(p$persons),p$roster),function(rows)
      length(rows)/length(p$persons)*tcrossprod(colMeans(V$scores[rows,])-global)))
    expect_lt(max(abs(random$meat-V$meat-between)),1e-12)
    # Finite perturbation of a Person's frequency within its fixed roster.
    rows <- which(p$roster==1L); i <- rows[which.max(rowSums(V$influence[rows,]^2))]
    direction <- numeric(length(p$persons)); direction[rows] <- -1/length(rows); direction[i] <- direction[i]+1
    step <- -drop(solve(V$jacobian,drop(crossprod(direction,V$scores))))
    h <- 1e-5
    corrected <- p$evaluate(b+h*step,k)$value
    shifted <- colMeans(corrected)+h*drop(crossprod(direction,corrected))
    expect_lt(max(abs(shifted-p$mean_score(b,k))),1e-8)
  }
  d <- x$data[nrow(x$data):1L,rev(names(x$data))]
  renamed <- do.call(mfrm_jml_adjustment_problem,c(list(data=d),x[c("person","facets","score","owner","rating_max")]))
  expect_identical(renamed$parameters,p$parameters)
  expect_identical(renamed$evaluate(x$results[[1L]]$beta,2L),p$evaluate(x$results[[1L]]$beta,2L))
  # Align with the existing native JML free coordinates and response kernel.
  expect_warning(prep <- prepare_mfrm_data(x$data,x$person,x$facets,x$score,rating_min=0,
    rating_max=x$rating_max,keep_original=TRUE),"duplicated Person")
  config <- build_estimation_config(prep,"GPCM","JML",x$owner,x$owner,
    weight_col=NULL,facet_signs=setNames(rep(-1,length(x$facets)),x$facets),
    positive_facets=character(),noncenter_facet="Person",dummy_facets=character(),
    anchor_df=NULL,group_anchor_df=NULL)$config
  b <- x$results[[1L]]$beta; raw <- p$evaluate(b,0L)
  person_theta <- raw$theta[match(prep$levels$Person,p$persons)]
  person_theta[!is.finite(person_theta)] <- 0
  params <- expand_params(c(person_theta,b),build_param_sizes(config),config)
  idx <- build_indices(prep,x$owner,x$owner)
  eta <- params$theta[idx$person]+compute_base_eta(idx,params,config)
  cum <- t(apply(params$steps_mat,1L,function(z)c(0,cumsum(z))))
  probabilities <- category_prob_gpcm(eta,cum,idx$step_idx,params$slopes,idx$slope_idx)
  nll <- -log(probabilities[cbind(seq_len(nrow(prep$data)),prep$data$score_k+1L)])
  native <- as.numeric(rowsum(nll,idx$person))
  finite <- is.finite(raw$theta)
  expect_equal(native[finite],raw$q[finite],tolerance=1e-10)
})

test_that("transition reuse respects coordinate changes and a small workspace", {
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  p <- do.call(mfrm_jml_adjustment_problem,x[c("data","person","facets","score","owner","rating_max")])
  beta <- x$results[["4"]]$beta
  shifted <- beta+seq_along(beta)*1e-3
  withr::local_options(mfrmr.max_information_bytes=256*1024^2)
  original <- p$evaluate(beta,4L)
  changed <- p$evaluate(shifted,4L)
  expect_gt(max(abs(original$value-changed$value)),1e-5)
  expect_identical(p$evaluate(beta,4L),original)
  # The existing problem is already constructed. A smaller allowance must
  # disable the optional transition cache without changing its calculation.
  withr::local_options(mfrmr.max_information_bytes=1)
  expect_identical(p$evaluate(beta,4L),original)
  expect_identical(p$evaluate(shifted,4L),changed)
})

test_that("binary categories and absent owners retain singleton roster dimensions", {
  # Three additive facets; each roster omits one owner level. This design is
  # deliberately confounded, but its profile equations still have exact limits.
  patterns <- as.matrix(expand.grid(rep(list(0:1),2L)))
  d <- do.call(rbind,lapply(1:8,function(i) data.frame(
    Person=paste0("p",i),Judge=if(i<=4) "A" else "B",Criterion=c("C","D"),
    Task=if(i<=4) c("X","Y") else c("Y","X"),Score=patterns[1L+(i-1L)%%4L,])))
  args <- list(data=d,person="Person",facets=c("Judge","Criterion","Task"),
    score="Score",owner="Judge",rating_max=1L)
  p <- do.call(mfrm_jml_adjustment_problem,args)
  expect_equal(unname(p$total_states),c(3L,3L))
  expect_false(any(p$parameters$Type=="step"))
  raw <- p$evaluate(rep(0,4L),0L)
  expect_equal(raw$q,rep(c(0,2*log(2),2*log(2),0),2L),tolerance=1e-12)
  expect_equal(dim(p$evaluate(rep(0,4L),4L)$value),c(8L,4L))
  one <- args; one$data <- d[d$Person %in% c("p2","p6"),]
  singleton <- do.call(mfrm_jml_adjustment_problem,one)
  expect_equal(singleton$evaluate(rep(0,4L),4L)$value,
    p$evaluate(rep(0,4L),4L)$value[c(2L,6L),,drop=FALSE])
  local_mocked_bindings(mfrmr_numeric_transformation_jacobian = function(...)
    stop("Numerical derivative requested."), .package = "mfrmr")
  expect_error(mfrm_jml_adjustment_covariance(singleton,rep(0,4L),4L),
    "at least two Persons")
  expect_error(mfrm_jml_adjustment_covariance(singleton,rep(0,4L),4L,
    sampling="random_rosters"), "Numerical derivative requested")
})

test_that("unsupported or uninformative calculations retain explicit failures", {
  x <- readRDS(test_path("fixtures","jml-adjustment-reference.rds"))[[1L]]
  args <- list(data=x$data,person="Person",facets=c("Rater","Criterion"),score="Score",owner=x$owner,rating_max=2L)
  expect_error(do.call(mfrm_jml_adjustment_problem,c(args,list(max_states=2L))),"exceeds max_states")
  bad <- args; bad$data$Score[1L] <- NA_real_
  expect_error(do.call(mfrm_jml_adjustment_problem,bad),"observed integer")
  p <- do.call(mfrm_jml_adjustment_problem,args)
  expect_error(p$evaluate(x$beta,1.5),"correction order")
  expect_error(mfrm_jml_adjustment_covariance(p,x$beta,1L),"adjusted mean equation")
  # All-minimum responses have limiting zero profile scores, not finite abilities.
  zero <- args; zero$data$Score[] <- 0
  z <- do.call(mfrm_jml_adjustment_problem,zero)
  expect_true(all(z$evaluate(x$beta,2L)$theta == -Inf))
  expect_true(all(z$evaluate(x$beta,2L)$value == 0))
  expect_error(mfrm_jml_adjustment_covariance(z,x$beta,2L),"full equation Jacobian")
  withr::local_options(mfrmr.max_information_bytes=1)
  expect_error(do.call(mfrm_jml_adjustment_problem,args),"workspace exceeds")
})

test_that("explicit-order fitting reproduces saved roots and native parameter tables", {
  skip_if_not_installed("nleqslv")
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  p <- do.call(mfrm_jml_adjustment_problem,x[c("data","person","facets","score","owner","rating_max")])
  expect_warning(prep <- prepare_mfrm_data(x$data,x$person,x$facets,x$score,
    rating_min=0,rating_max=x$rating_max,keep_original=TRUE),"duplicated Person")
  config <- build_estimation_config(prep,"GPCM","JML",x$owner,x$owner,
    weight_col=NULL,facet_signs=setNames(rep(-1,length(x$facets)),x$facets),
    positive_facets=character(),noncenter_facet="Person",dummy_facets=character(),
    anchor_df=NULL,group_anchor_df=NULL)$config
  sizes <- build_param_sizes(config)
  native <- function(b) expand_params(c(rep(0,sizes$theta),b),sizes,config)
  expanded <- function(b) {
    q <- native(b)
    c(unlist(q$facets,use.names=FALSE),as.vector(t(q$steps_mat)),q$log_slopes)
  }
  for(k in c(2L,4L)) {
    z <- mfrm_jml_adjustment_fit(p,k)
    expect_true(z$point$available && z$covariance$available)
    expect_identical(z$point$status,"consistent_roots")
    expect_equal(z$point$beta,x$results[[as.character(k)]]$beta,tolerance=1e-7)
    expect_equal(z$covariance$result$vcov,x$results[[as.character(k)]]$covariance$vcov,tolerance=1e-8)
    expect_identical(z$estimator$order,k)
    expect_identical(z$estimator$order_selection,"explicit")
    expect_false(z$estimator$formal_structural_intervals)
    expect_length(z$attempts,2L)
    expect_true(all(vapply(z$attempts,function(a)a$newton_step<1e-5,logical(1))))
    expect_equal(z$point$theta,p$evaluate(z$point$beta,k)$theta)
    tables <- mfrm_jml_adjustment_tables(z)
    actual <- c(tables$facets$others$Estimate,tables$steps$Estimate,tables$slopes$LogEstimate)
    expect_equal(actual,expanded(z$point$beta),tolerance=1e-12)
    J <- mfrmr_numeric_transformation_jacobian(expanded,z$point$beta)$jacobian
    expected_se <- sqrt(diag(J%*%z$covariance$result$vcov%*%t(J)))
    expect_equal(c(tables$facets$others$RootSE,tables$steps$RootSE,tables$slopes$LogRootSE),
      expected_se,tolerance=1e-8)
    expect_equal(tables$slopes$RootSE,tables$slopes$Estimate*tables$slopes$LogRootSE)
    expect_true(all(is.na(tables$facets$person$SE)))
    expect_true(all(!tables$slopes$CIEligible))
    expect_false(inherits(tables,"mfrm_fit"))
    file <- tempfile(fileext=".rds"); on.exit(unlink(file),add=TRUE)
    saveRDS(list(result=z,tables=tables),file)
    expect_identical(readRDS(file),list(result=z,tables=tables))
  }
})

test_that("covariance failure cannot erase a root or trigger another fitting attempt", {
  skip_if_not_installed("nleqslv")
  x <- readRDS(test_path("fixtures","jml-adjustment-general.rds"))
  # Removing one assigned rating creates a roster containing only one Person.
  # A fixed-roster empirical variance is unavailable for that stratum.
  x$data <- x$data[-1L,]
  p <- do.call(mfrm_jml_adjustment_problem,x[c("data","person","facets","score","owner","rating_max")])
  expect_true(any(table(p$roster)==1L))
  z <- mfrm_jml_adjustment_fit(p,2L)
  random <- mfrm_jml_adjustment_fit(p,2L,sampling="random_rosters")
  expect_true(z$point$available)
  expect_false(z$covariance$available)
  expect_true(random$covariance$available)
  expect_identical(z$point,random$point)
  expect_identical(z$attempts,random$attempts)
  expect_match(z$covariance$reason,"at least two Persons")
  tables <- mfrm_jml_adjustment_tables(z)
  expect_true(all(is.finite(tables$slopes$Estimate)))
  expect_true(all(is.na(tables$slopes$RootSE)))
  expect_true(all(grepl("at least two Persons",tables$slopes$CovarianceReason)))
  # Random-roster inference has a different sampling assumption; this comparison
  # checks that covariance policy cannot select different point estimates.
})

test_that("failed and ambiguous equations retain attempts without choosing by likelihood", {
  skip_if_not_installed("nleqslv")
  # Synthetic equations isolate solver decisions; they are not GPCM evidence.
  p <- list(parameters=data.frame(Type="location",Step=NA_integer_),
    specification=list(rating_max=1L),persons=paste0("p",1:4),levels=list(),roster=rep(1L,4))
  p$mean_score <- function(b,order) b^2-1
  p$evaluate <- function(b,order) list(value=matrix(rep(b^2-1,4),4,1),q=rep(b^2,4),theta=rep(0,4))
  roots <- mfrm_jml_adjustment_fit(p,2L,starts=list(left=-.8,right=.8))
  expect_identical(roots$point$status,"multiple_roots")
  expect_false(roots$point$available || roots$covariance$available)
  expect_null(roots$point$beta)
  expect_true(all(vapply(roots$attempts,`[[`,logical(1),"available")))
  expect_error(mfrm_jml_adjustment_tables(roots),"unambiguous")
  p$mean_score <- function(b,order) {
    if(b>10) stop("Starting region cannot be evaluated.")
    b^2-1
  }
  partial <- mfrm_jml_adjustment_fit(p,2L,starts=list(left=.8,right=100))
  expect_identical(partial$point$status,"root_with_unresolved_starts")
  expect_true(partial$point$available)
  expect_false(partial$covariance$available)
  expect_true(is.na(partial$point$root_difference))
  expect_length(partial$attempts,4L)
  expect_match(partial$attempts[[4L]]$error,"Starting region")
  limited <- mfrm_jml_adjustment_fit(p,2L,starts=list(left=-.2,right=.2),maxit=1L)
  expect_false(limited$point$available)
  expect_length(limited$attempts,6L)
  expect_true(all(vapply(limited$attempts,function(a)!is.null(a$solver),logical(1))))
  # A small unscaled residual alone is insufficient when a full Newton step is large.
  p$mean_score <- function(b,order) 1e-12*(b-2)
  small <- mfrm_jml_adjustment_fit(p,2L,starts=list(left=0,right=1))
  expect_false(small$point$available)
  expect_true(all(vapply(small$attempts,function(a)a$reason=="root_step_unresolved",logical(1))))
  p$mean_score <- function(b,order) 0*b
  zero <- mfrm_jml_adjustment_fit(p,2L,starts=list(left=0,right=1))
  expect_false(zero$point$available)
  expect_true(all(vapply(zero$attempts,function(a)a$reason=="equation_rank_deficient",logical(1))))
  expect_error(mfrm_jml_adjustment_fit(p,"auto"),"explicit positive")
  expect_error(mfrm_jml_adjustment_fit(p,0L),"explicit positive")
  expect_error(mfrm_jml_adjustment_fit(p,2L,maxit=0L),"maxit")
  expect_error(mfrm_jml_adjustment_fit(p,2L,starts=list(0,0)),"distinct")
})

test_that("binary table expansion preserves zero steps and unbounded Person profiles", {
  # Only a coordinate-transform test; no estimator is claimed for these values.
  z <- list(point=list(available=TRUE,status="consistent_roots",
    beta=c(.2,-.3,.1,.15),theta=c(-Inf,0,Inf)),
    specification=list(facets=c("Judge","Criterion","Task"),owner="Judge",rating_max=1L),
    levels=list(Judge=c("A","B"),Criterion=c("C","D"),Task=c("X","Y")),
    persons=c("p1","p2","p3"),estimator=list(name="finite_MLE_plugin_profile_score_adjustment",
      order=2L,target="Adjusted-equation root; residual structural bias may remain."),
    covariance=list(available=TRUE,status="local_root_covariance",reason="",result=list(vcov=diag(4L))))
  t <- mfrm_jml_adjustment_tables(z)
  expect_equal(t$facets$others$Estimate,c(.2,-.2,-.3,.3,.1,-.1))
  expect_equal(t$steps$Estimate,c(0,0))
  expect_equal(t$steps$RootSE,c(0,0))
  expect_equal(t$slopes$Estimate,exp(c(.15,-.15)))
  expect_equal(t$slopes$LogRootSE,c(1,1))
  expect_identical(t$facets$person$Extreme,c("low","none","high"))
  expect_equal(t$facets$person$Estimate,c(-Inf,0,Inf))
  withr::local_options(mfrmr.max_information_bytes=1)
  expect_error(mfrm_jml_adjustment_tables(z),"transformation exceeds")
})
