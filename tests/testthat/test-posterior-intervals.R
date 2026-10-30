test_that("continuous posterior bounds recover a normal prior, including small tails", {
  idx <- list(person=c(2L,1L),score_k=c(0L,3L),weight=c(0,0),step_idx=c(1L,1L))
  config <- list(model="RSM",n_cat=4L,n_person=2L)
  params <- list(steps=c(-1,0,1))
  population <- list(active=TRUE,design_matrix=matrix(c(1,2),2L,1L),coefficients=-0.8,
    sigma2=0.16,person_lookup=1:2)
  for (level in c(0.8,0.95,1-1e-10)) {
    actual <- mfrmr:::mfrmr_person_posterior_intervals(idx,config,params,c(0,0),c("A","B"),
      population,level)
    alpha <- (1-level)/2
    expected <- outer(c(-1.6,-0.8),0.4*c(qnorm(alpha),qnorm(alpha,lower.tail=FALSE)),"+")
    expect_equal(unname(actual),unname(expected),tolerance=1e-8)
  }
})

test_that("posterior intervals and decision probabilities ignore constant likelihood offsets", {
  offset <- 0
  local_mocked_bindings(mfrmr_adaptive_person_kernel = function(...)
    function(z, ...) list(log_likelihood = rep(offset, length(z)),
      score = -z, information = rep(1, length(z))), .package = "mfrmr")
  population <- list(active = TRUE, design_matrix = matrix(1, 1L, 1L),
    coefficients = -.8, sigma2 = .16, person_lookup = 1L)
  cuts <- -.8 + .4 * c(-1, 0, 1)
  # An ability-independent likelihood leaves this normal prior unchanged.
  expected <- c(Lower = -.8 + .4 * qnorm(.025), Upper = -.8 + .4 * qnorm(.975),
    setNames(pnorm(c(-1, 0, 1), lower.tail = FALSE), paste0("Above", 1:3)),
    setNames(diff(pnorm(c(-Inf, -1, 0, 1, Inf))), paste0("Band", 1:4)))
  for (constant in c(0, -1e12, -1e16)) {
    offset <- constant
    actual <- mfrmr_person_posterior_intervals(
      list(person = 1L, score_k = 0L, weight = 1),
      list(model = "RSM", n_cat = 2L, n_person = 1L), list(steps = 0),
      0, "P", population, cut_scores = cuts)
    expect_equal(drop(actual), expected, tolerance = 1e-9)
  }
})

test_that("new-Person summaries preserve the prior under negligible ability dependence", {
  # Finite GPCM inputs give an almost ability-independent likelihood exp(-1e16).
  # Its absolute size must not change prior weights, means or SDs.
  idx <- list(person = c(2L, 1L), score_k = c(0L, 0L), weight = c(1, 1),
    step_idx = c(1L, 1L), slope_idx = c(1L, 1L), facets = list(Item = c(1L, 1L)))
  cfg <- list(model = "GPCM", n_cat = 2L, n_person = 2L, facet_names = "Item")
  params <- list(facets = list(Item = -1e32), steps_mat = matrix(0, 1L, 1L), slopes = 1e-16)
  pop <- list(active = TRUE, design_matrix = matrix(c(1, 2), 2L, 1L),
    coefficients = .4, sigma2 = .81, person_lookup = 1:2)
  for (integration in c("fixed", "adaptive")) {
    cfg$estimation_control <- list(mml_integration = integration)
    out <- compute_person_posterior_summary(idx, cfg, params, gauss_hermite_normal(7L),
      c("A", "B"), pop, n_draws = 3L, seed = 421L)
    expect_identical(out$estimates$Person, c("B", "A"))
    expect_equal(unname(out$estimates$Estimate), c(.8, .4), tolerance = 1e-12)
    expect_equal(unname(out$estimates$SD), c(.9, .9), tolerance = 1e-12)
    expect_equal(out$estimates$Lower, c(.8, .4) + .9 * qnorm(.025), tolerance = 1e-9)
    expect_equal(out$estimates$Upper, c(.8, .4) + .9 * qnorm(.975), tolerance = 1e-9)
    expect_identical(out$draws$Person, rep(c("B", "A"), each = 3L))
  }
})

test_that("posterior endpoints carry their stated mass for weighted ordinal likelihoods", {
  for (model in c("RSM","PCM","GPCM")) for (n in c(4L,40L)) {
    owner <- rep(1:2,length.out=n)
    idx <- list(person=rep(1L,n),score_k=if(n==4L) 0:3 else rep(3L,n),
      weight=rep(c(0.25,1,2,0),length.out=n),step_idx=owner,slope_idx=owner)
    config <- list(model=model,n_cat=4L,n_person=1L)
    params <- list(steps=c(-0.8,0,0.8),steps_mat=rbind(c(-1.2,0.3,0.9),c(-0.4,-0.1,0.5)),
      slopes=c(0.6,1.5))
    base <- rep(c(-0.8,0.8),length.out=n)
    actual <- mfrmr:::mfrmr_person_posterior_intervals(idx,config,params,base,"P")
    steps <- if(model=="RSM") matrix(params$steps,n,3L,byrow=TRUE) else params$steps_mat[owner,]
    cumulative <- t(apply(steps,1L,function(x) c(0,cumsum(x))))
    slope <- if(model=="GPCM") params$slopes[owner] else rep(1,n)
    log_density <- function(theta) vapply(theta,function(value) {
      logits <- slope*(outer(value+base,0:3)-cumulative)
      high <- apply(logits,1L,max)
      sum(idx$weight*(logits[cbind(seq_len(n),idx$score_k+1L)]-high-
        log(rowSums(exp(logits-high)))))+dnorm(value,log=TRUE)
    },0)
    peak <- optimize(log_density,c(-12,12),maximum=TRUE)$objective
    density <- function(theta) exp(log_density(theta)-peak)
    mass <- integrate(density,-12,12,rel.tol=1e-11)$value
    cdf <- vapply(actual,function(value) integrate(density,-12,value,rel.tol=1e-11)$value/mass,0)
    expect_equal(unname(cdf),c(0.025,0.975),tolerance=1e-8)
  }
})

test_that("likelihood-only evaluation preserves weighted and shifted-prior kernels", {
  old <- options(mfrmr.use_cpp11_backend=TRUE)
  on.exit(options(old),add=TRUE)
  idx <- list(person=rep(1L,4L),score_k=0:3,weight=c(0,0.25,1,2),
    step_idx=c(2L,1L,2L,1L),slope_idx=c(2L,1L,2L,1L))
  params <- list(steps=c(-0.8,0,0.8),steps_mat=rbind(c(-1.2,0.3,0.9),c(-0.4,-0.1,0.5)),
    slopes=c(0.6,1.5))
  z <- c(-25,-3,0,0.4,4,25)
  for (model in c("RSM","PCM","GPCM")) {
    config <- list(model=model,n_cat=4L,n_person=1L,
      estimation_control=list(mml_integration="adaptive"))
    evaluate <- mfrmr:::mfrmr_adaptive_person_kernel(idx,config,params,c(-1,0,1,0.3),-0.8,0.4)
    expected <- evaluate(z)
    for (compiled in c(TRUE,FALSE)) {
      options(mfrmr.use_cpp11_backend=compiled)
      actual <- evaluate(z,include_moments=FALSE)
      expect_equal(actual$log_likelihood,expected$log_likelihood,tolerance=1e-12)
      expect_identical(evaluate(z),expected)
    }
  }
})
