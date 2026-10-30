# Use accepted fits; extracting and scoring a calibration must not refit them.
anchored_calibration_rows <- function() {
  ratings <- load_mfrmr_data('example_core')
  rows <- ratings[ratings$Person == unique(ratings$Person)[1L], ]
  low <- high <- rows
  rows$Person <- 'Typical'; low$Person <- 'Lower'; high$Person <- 'Upper'
  low$Score <- min(ratings$Score); high$Score <- max(ratings$Score)
  rbind(rows, low, high)
}

# Direct RSM/PCM category probabilities and continuous normal integration,
# independent of the production scoring and posterior-interval helpers.
anchored_calibration_reference <- function(fit, rows) {
  levels <- fit$config$facet_names
  offset <- Reduce(`+`, lapply(levels, function(facet) {
    tab <- fit$facets$others[fit$facets$others$Facet == facet, ]
    fit$config$facet_signs[[facet]] * tab$Estimate[match(rows[[facet]], tab$Level)]
  }))
  K <- fit$config$n_cat - 1L
  steps <- lapply(seq_len(nrow(rows)), function(i) {
    tab <- fit$steps
    if (fit$config$model == 'PCM') tab <- tab[tab$StepFacet == rows[[fit$config$step_facet]][i], ]
    c(0, cumsum(tab$Estimate))
  })
  y <- match(rows$Score, fit$prep$score_map$OriginalScore) - 1L
  logdensity <- function(theta) vapply(theta, function(t) {
    value <- vapply(seq_len(nrow(rows)), function(i) {
      z <- (0:K)*(t + offset[i]) - steps[[i]]
      z[y[i]+1L] - max(z) - log(sum(exp(z-max(z))))
    }, numeric(1))
    sum(value) + dnorm(t, log=TRUE)
  }, numeric(1))
  peak <- optimize(logdensity, c(-12,12), maximum=TRUE)$objective
  density <- function(t) exp(logdensity(t)-peak)
  integrate_density <- function(lower,upper) integrate(density,lower,upper,
    rel.tol=1e-10, abs.tol=0)$value
  mass <- integrate_density(-Inf,Inf)
  mean <- integrate(function(t) t*density(t),-Inf,Inf,rel.tol=1e-10)$value/mass
  sd <- sqrt(integrate(function(t) (t-mean)^2*density(t),-Inf,Inf,rel.tol=1e-10)$value/mass)
  quantile <- function(p) uniroot(function(t) integrate_density(-Inf,t)/mass-p,
    c(-12,12), tol=1e-9)$root
  c(Estimate=mean, SD=sd, Lower=quantile(.025), Upper=quantile(.975),
    Pass=integrate_density(0,Inf)/mass)
}

local({
  fixtures <- readRDS(test_path('fixtures','anchored-facet-information.rds'))
  for (case in names(fixtures)) test_that(paste(case, 'calibration preserves anchors and conditional scoring'), {
    fit <- fixtures[[case]]$fit
    local_mocked_bindings(fit_mfrm=function(...) stop('No refitting'), .package='mfrmr')
    draft <- mfrmr_extract_calibration_draft(fit, scoring_quad_points=161L)
    id <- draft$constraints$identification
    expect_identical(id$ConstraintType[1:2], c('typed_anchor_expansion_facet','uncentered'))
    expect_identical(id$FreeDimension[1:2], as.integer(c(if (case=='RSM_group') 1 else 2,4)))
    expect_equal(nrow(review_mfrm_calibration(draft)),0)
    frozen <- freeze_mfrm_calibration(validate_mfrm_calibration(draft))
    path <- tempfile(fileext='.rds'); on.exit(unlink(path),add=TRUE)
    save_mfrm_calibration(frozen,path)
    reopened <- load_mfrm_calibration(path)
    expect_identical(reopened,frozen)
    expect_identical(reopened$constraints,draft$constraints)
    expect_identical(mfrmr_calibration_find_prohibited(reopened),character())
    coordinates <- reopened$parameters$coordinates
    stored <- coordinates[coordinates$ParameterClass=='facet',]
    source <- fit$facets$others
    expect_equal(stored$Value,unname(source$Estimate[match(paste(stored$OwnerFacet,stored$Level),
      paste(source$Facet,source$Level))]),tolerance=0)
    if (case=='RSM_group') {
      rater <- stored[stored$OwnerFacet=='Rater',]
      expect_equal(rater$Value[4],-.3,tolerance=1e-12)
      expect_equal(mean(rater$Value[1:3]),.1,tolerance=1e-12)
    }
    rows <- anchored_calibration_rows()
    scored <- score_mfrm_calibration(reopened,rows,cut_scores=c(Pass=0))
    fitted <- predict_mfrm_units(fit,rows,scoring_quad_points=161L,cut_scores=c(Pass=0))
    columns <- c('Estimate','SD','Lower','Upper')
    expect_equal(as.matrix(scored$estimates[columns]),
      unname(as.matrix(fitted$estimates[match(scored$estimates$Person,fitted$estimates$Person),columns])),
      ignore_attr=TRUE,tolerance=1e-10)
    expect_identical(scored$estimates$UncertaintyBasis,
      rep('conditional_on_frozen_point_calibration',3L))
    for (person in scored$estimates$Person) {
      ref <- anchored_calibration_reference(fit,rows[rows$Person==person,])
      expect_equal(unname(unlist(scored$estimates[scored$estimates$Person==person,columns])),
        unname(ref[columns]),tolerance=1e-7)
      probability <- scored$decisions$Probability[scored$decisions$Person==person &
        scored$decisions$Type=='above_cut']
      expect_equal(probability,unname(ref['Pass']),tolerance=1e-8)
    }
    bad <- draft
    bad$constraints$identification$FreeDimension[1] <- 20L
    expect_true('IDENTIFICATION_CONTRACT_INVALID' %in% review_mfrm_calibration(bad)$Code)
    bad <- draft
    bad$constraints$identification$ConstraintType[1] <- 'typed_anchor_expansion_remainder'
    expect_true('IDENTIFICATION_CONTRACT_INVALID' %in% review_mfrm_calibration(bad)$Code)
    # Freeze binds the convention as well as the numeric coordinates.
    bad <- frozen; bad$constraints$identification$ConstraintType[1] <- 'typed_anchor_expansion'
    expect_error(score_mfrm_calibration(bad,rows),'IDENTITY_')
    # This separate, narrower scoring entry must retain its anchor exclusion.
    if (case != 'PCM') expect_error(score_mfrm_persons(fit),'unanchored')
  })
})

test_that('explicit extraction distinguishes legacy and uncentered maps without rewriting old artifacts', {
  fit <- readRDS(test_path('fixtures','anchored-facet-information.rds'))$RSM$fit
  for (centered in c(TRUE,FALSE)) {
    spec <- build_facet_constraint(c('A','B','C'),anchors=c(A=.4),
      centered=centered,centering_scope='remainder')
    # Algebra-only configuration, not a claim about a newly estimated fit.
    dummy <- fit; dummy$config$facet_specs$Rater <- spec
    id <- mfrmr_calibration_extract_identification(dummy)
    expect_identical(id$ConstraintType[1],if (centered)
      'typed_anchor_expansion_remainder' else 'typed_anchor_expansion_uncentered')
    dummy$config$facet_specs$Rater$centering_scope <- NULL
    expect_identical(mfrmr_calibration_extract_identification(dummy),id)
  }
  old <- mfrmr_extract_calibration_draft(fit)
  old$constraints$identification$ConstraintType[1:2] <- c('typed_anchor_expansion','sum_to_zero')
  old$integrity$semantic_components <- mfrmr_calibration_semantic_components(old)
  # Old files' coarse labels cannot prove centering; retain the coordinates
  # and labels rather than treating a load as migration or re-estimation.
  old <- freeze_mfrm_calibration(validate_mfrm_calibration(old))
  path <- tempfile(fileext='.rds'); on.exit(unlink(path),add=TRUE)
  save_mfrm_calibration(old,path)
  expect_identical(load_mfrm_calibration(path),old)
})
