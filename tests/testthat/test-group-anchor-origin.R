# Independent category likelihood, using neither package probabilities nor scores.
group_origin_likelihood <- function(data, coords, steps) {
  eta <- coords$Person[data$Person] - coords$Rater[as.character(data$Rater)] -
    coords$Task[as.character(data$Task)]
  ladder <- if (nrow(steps) == 1L) steps[rep(1L, nrow(data)), , drop = FALSE] else
    steps[as.integer(data$Task), , drop = FALSE]
  logits <- cbind(0, eta-ladder[,1], 2*eta-rowSums(ladder[,1:2,drop=FALSE]),
                  3*eta-rowSums(ladder))
  lp <- logits - apply(logits, 1, max)
  lp <- lp - log(rowSums(exp(lp)))
  expected <- drop(exp(lp) %*% (0:3))
  score <- data$Score - expected
  r <- -tapply(score, data$Rater, sum)
  list(loglik = sum(lp[cbind(seq_len(nrow(data)), data$Score+1L)]),
       expected = expected, group_score = unname(r['1']-r['2']))
}

test_that("group origins preserve the declared likelihood and expose native limits", {
  input <- read.csv(test_path('fixtures', 'tam-jml-0.2.4-response.csv'))
  reference <- read.csv(test_path('fixtures', 'facets-group-anchor.csv'))
  for (model in c('RSM','PCM')) for (variant in c('origin','restricted')) {
    d <- input[startsWith(input$Case, model), ]
    group <- data.frame(Facet='Rater', Level=c('1','2'), Group='g', GroupValue=.6)
    if (variant == 'restricted') group <- rbind(group,
      data.frame(Facet='Rater', Level='3', Group='h', GroupValue=-.2))
    fit <- fit_mfrm(d, 'Person', c('Rater','Task'), 'Score', model=model,
      step_facet=if(model=='PCM') 'Task' else NULL, method='JML',
      rating_min=0, rating_max=3, group_anchors=group,
      centering_scope='anchor_origin', optimizer='BFGS', maxit=500, reltol=1e-11)
    pa <- expand_params(fit$opt$par, build_param_sizes(fit$config), fit$config)
    coords <- c(list(Person=pa$theta), pa$facets)
    steps <- if(model=='PCM') pa$steps_mat else matrix(pa$steps,nrow=1)
    own <- group_origin_likelihood(d, coords, steps)
    expect_equal(own$loglik, -fit$opt$value, tolerance=1e-10)
    expect_lt(abs(own$group_score), 1e-4)
    expect_equal(mean(coords$Rater[c('1','2')]), .6, tolerance=1e-12)
    expect_equal(sum(coords$Task), 0, tolerance=1e-12)
    expect_false(fit$config$facet_specs$Rater$centered)
    expect_identical(fit$config$replay_inputs$centering_scope, 'anchor_origin')
    fs <- fit$config$anchor_review$facet_summary
    expect_identical(fs$OriginBasis[fs$Facet=='Rater'], 'group_anchors')
    expect_false(fs$Centered[fs$Facet=='Rater'])
    expect_equal(fit$config$facet_specs$Rater$n_params,
                 if(variant=='origin') 2L else 1L)

    ref <- reference[reference$Case==paste(model,variant,sep='_'), ]
    native <- setNames(lapply(c('Person','Rater','Task'), function(f) {
      rows <- ref[ref$Facet==f, ];setNames(rows$Native,rows$Level)
    }), c('Person','Rater','Task'))
    native_steps <- do.call(rbind,lapply(seq_len(nrow(steps)),function(i)
      ref$Native[ref$Facet==paste0('Step',i)]))
    check <- group_origin_likelihood(d,native,native_steps)
    if (variant == 'origin') {
      for(f in names(coords)) expect_equal(unname(coords[[f]]),
        unname(native[[f]][names(coords[[f]])]), tolerance=1e-4)
      expect_equal(unname(steps),native_steps,tolerance=1e-4)
      # A Person-group origin is another translation of the same predictors.
      person_group <- data.frame(Facet='Person',Level=c('P03','P04'),
        Group='p',GroupValue=-.4)
      pfit <- fit_mfrm(d,'Person',c('Rater','Task'),'Score',model=model,
        step_facet=if(model=='PCM') 'Task' else NULL,method='JML',
        rating_min=0,rating_max=3,group_anchors=person_group,
        noncenter_facet='Rater',centering_scope='anchor_origin',
        optimizer='BFGS',maxit=500,reltol=1e-11)
      pp <- expand_params(pfit$opt$par,build_param_sizes(pfit$config),pfit$config)
      shift <- -.4-mean(coords$Person[c('P03','P04')])
      expect_equal(pp$theta,coords$Person+shift,tolerance=1e-4)
      expect_equal(pp$facets$Rater,coords$Rater+shift,tolerance=1e-4)
      expect_equal(pfit$opt$value,fit$opt$value,tolerance=1e-8)
    } else {
      expect_equal(unname(coords$Rater['3']),-.2,tolerance=1e-12)
      # Both native group targets hold, yet a feasible direction improves LL.
      expect_equal(mean(native$Rater[c('1','2')]),.6,tolerance=1e-12)
      moved <- native
      moved$Rater[c('1','2')] <- moved$Rater[c('1','2')]+c(.001,-.001)
      gain <- group_origin_likelihood(d,moved,native_steps)$loglik-check$loglik
      expect_gt(gain,1e-4)
      expect_gt(abs(check$group_score),.1)
      expect_gt(own$loglik-check$loglik,1e-4)
    }
    saved <- unserialize(serialize(fit,NULL))
    expect_identical(saved$config$facet_specs,fit$config$facet_specs)
    recipe <- build_mfrm_replay_script(saved)
    expect_match(recipe$script,'centering_scope = "anchor_origin"',fixed=TRUE)
    expect_match(recipe$script,'group_anchors',fixed=TRUE)
  }
})

test_that("group-origin review preserves overlap and rejects incompatible declarations", {
  d <- load_mfrmr_data('example_core')
  group <- data.frame(Facet='Rater',Level=c('R01','R02'),Group='g',GroupValue=.6)
  args <- list(data=d,person='Person',facets=c('Rater','Criterion'),score='Score',
    group_anchors=group,centering_scope='anchor_origin')
  review <- do.call(review_mfrm_anchors,args)
  # Review uses the same affine equations without requiring a fit.
  prep <- prepare_mfrm_data(d,'Person',c('Rater','Criterion'),'Score')
  direct <- data.frame(Facet='Rater',Level='R01',Anchor=.2)
  args$anchors <- direct
  both <- do.call(review_mfrm_anchors,args)
  expect_identical(both$facet_summary$OriginBasis[both$facet_summary$Facet=='Rater'],
                   'direct_and_group_anchors')
  specs <- prepare_constraint_specs(prep,anchor_review=both,centering_scope='anchor_origin')
  s <- specs$facet_specs$Rater
  values <- expand_facet_with_constraints(numeric(s$n_params),s)
  expect_equal(unname(values[c('R01','R02')]),c(.2,1),tolerance=1e-12)
  args$anchors <- rbind(direct,data.frame(Facet='Rater',Level='R02',Anchor=.2))
  expect_error(do.call(review_mfrm_anchors,args),'conflict with the mean')
  args$anchors <- NULL; args$noncenter_facet <- 'Rater'
  expect_error(do.call(review_mfrm_anchors,args),'unanchored, active')
  # The existing whole-facet convention continues to balance the group total.
  args$noncenter_facet <- 'Person';args$centering_scope <- 'facet'
  whole <- do.call(review_mfrm_anchors,args)
  s <- prepare_constraint_specs(prep,anchor_review=whole,centering_scope='facet')$facet_specs$Rater
  values <- expand_facet_with_constraints(numeric(s$n_params),s)
  expect_equal(sum(values),0,tolerance=1e-12)
  expect_equal(mean(values[c('R01','R02')]),.6,tolerance=1e-12)
  # Review must expose unequal values as a conflict, not import contributions.
  args$group_anchors$GroupValue <- c(.2,1)
  bad <- do.call(review_mfrm_anchors,args)
  expect_equal(nrow(bad$issues$group_value_conflicts),1L)
})
