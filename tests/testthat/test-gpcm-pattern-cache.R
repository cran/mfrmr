test_that('fixed-grid pattern scoring reuses cumulative probabilities without changing scores', {
  d <- simulate_mfrm_data(n_person=10,n_rater=2,n_criterion=2,score_levels=3,seed=92451)
  prep <- mfrmr:::prepare_mfrm_data(d,person_col='Person',facet_cols=c('Rater','Criterion'),score_col='Score')
  signs <- mfrmr:::build_facet_signs(prep$facet_names)
  cfg <- mfrmr:::build_estimation_config(prep=prep,model='GPCM',method='MML',
    step_facet='Rater',slope_facet='Rater',weight_col=NULL,
    facet_signs=signs$signs,positive_facets=signs$positive_facets,
    noncenter_facet='Person',dummy_facets=character(),anchor_df=NULL,group_anchor_df=NULL,population=NULL)
  par <- mfrmr:::build_initial_param_vector(cfg$config,cfg$sizes)
  idx <- mfrmr:::build_indices(prep,step_facet='Rater',slope_facet='Rater')
  idx <- mfrmr:::mfrmr_subset_observation_indices(idx,which(idx$person==1L))
  patterns <- mfrmr:::mfrmr_enumerate_response_patterns(length(idx$score_k),3)
  quad <- mfrmr:::gauss_hermite_normal(7L)
  evaluate <- function(p=par,patterns_arg=patterns,include_scores=TRUE)
    mfrmr:::mfrmr_mml_evaluate_person_patterns(p,idx,cfg$config,cfg$sizes,quad,patterns_arg,include_scores)
  original <- mfrmr:::compute_P_geq; count <- 0L
  local_mocked_bindings(compute_P_geq=function(x) {count <<- count+1L; original(x)},.package='mfrmr')
  cached <- evaluate()
  expect_equal(count,7L)
  # The previous gradient path is the direct oracle, retaining all posterior
  # calculations but removing the optional cumulative-probability payload.
  gradient <- mfrmr:::mfrm_grad_mml_core
  uncached <- function(...) {
    args <- list(...); args$logprob_bundle$p_geq_list <- NULL
    do.call(gradient,args)
  }
  testthat::with_mocked_bindings({
    before <- count
    reference <- evaluate()
    expect_identical(cached,reference)
    expect_equal(count-before,7L*(nrow(patterns)+1L))
  },mfrm_grad_mml_core=uncached,.package='mfrmr')
  order <- rev(seq_len(nrow(patterns)))
  reversed <- evaluate(patterns_arg=patterns[order,,drop=FALSE])
  expect_identical(reversed$score,cached$score[order,,drop=FALSE])
  expect_identical(reversed$log_marginal,cached$log_marginal[order])
  changed <- par; changed[1] <- changed[1]+.2
  changed_result <- evaluate(p=changed)
  expect_false(isTRUE(all.equal(changed_result$score,cached$score)))
  testthat::with_mocked_bindings({
    expect_identical(evaluate(p=changed),changed_result)
  },mfrm_grad_mml_core=uncached,.package='mfrmr')
  before <- count; probabilities <- evaluate(include_scores=FALSE)
  expect_equal(count,before)
  expect_identical(probabilities$log_marginal,cached$log_marginal)
})
