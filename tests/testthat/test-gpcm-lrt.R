gpcm_lrt_fixture <- function(owner = "Criterion", slope_owner = owner) {
  data <- expand.grid(Person = paste0("P", 1:8), Rater = paste0("R", 1:3),
                      Criterion = paste0("C", 1:4), stringsAsFactors = FALSE)
  data$Score <- (seq_len(nrow(data)) - 1L) %% 3L
  prep <- mfrmr:::prepare_mfrm_data(data, "Person", c("Rater", "Criterion"), "Score")
  lapply(c("PCM", "GPCM"), function(model) {
    config <- mfrmr:::build_estimation_config(
      prep, model, "MML", owner, if (model == "GPCM") slope_owner else NULL,
      weight_col = NULL, positive_facets = character(0), noncenter_facet = "Person",
      dummy_facets = character(0), anchor_df = NULL, group_anchor_df = NULL,
      facet_signs = mfrmr:::build_facet_signs(prep$facet_names)$signs
    )$config
    config$population_spec <- list(active = TRUE, design_matrix = matrix(1, 8, 1),
      design_columns = "(Intercept)", person_lookup = 1:8)
    config$estimation_control <- list(quad_points = 31L, mml_integration = "adaptive")
    sizes <- mfrmr:::build_param_sizes(config)
    par <- rep(0, sum(unlist(sizes)))
    opt <- list(par = par, value = if (model == "PCM") 100 else 97)
    contract <- mfrmr:::build_mfrm_ic_contract(-opt$value, length(par), prep, config, "MML")
    structure(list(config = config, prep = prep, opt = opt, summary = contract), class = "mfrm_fit")
  })
}

test_that("comparison identities change only for the new anchored convention", {
  old <- gpcm_lrt_fixture()[[1]]
  current <- old
  current$config$centering_scope <- "facet"
  expect_identical(normalize_compare_signature(current), normalize_compare_signature(old))
  for (scope in c("remainder", "facet")) {
    fit <- old
    fit$config$centering_scope <- scope
    fit$config$facet_specs$Rater <- build_facet_constraint(fit$prep$levels$Rater,
      anchors = c(R1 = .7), centering_scope = scope)
    signature <- normalize_compare_signature(fit)
    if (scope == "remainder") {
      expect_null(signature$centering_scope)
      legacy <- fit
      legacy$config$centering_scope <- NULL
      expect_identical(normalize_compare_signature(legacy), signature)
    } else {
      expect_identical(signature$centering_scope, "facet")
      expect_match(audit_compare_mfrm_nesting(list(legacy, fit), c("old", "new"))$reason,
        "centering_scope")
    }
  }
})

test_that("PCM/GPCM nesting isolates relative slopes and preserves population design", {
  audit <- function(pair) mfrmr:::audit_compare_mfrm_nesting(pair, c("first", "second"))
  for (owner in c("Criterion", "Rater")) {
    pair <- gpcm_lrt_fixture(owner)
    result <- audit(pair)
    expect_true(result$eligible)
    expect_identical(result$relation, "PCM_in_GPCM")
    expect_identical(result$df, length(pair[[2]]$config$gpcm_spec$levels) - 1L)
    expect_identical(audit(rev(pair))$simpler, "second")
    changed <- pair
    changed[[1]]$config$population_spec$active <- FALSE
    expect_identical(audit(changed)$relation, "population_model_unverified")
    changed <- pair
    changed[[2]]$config$population_spec$design_matrix[1, 1] <- 2
    expect_false(audit(changed)$eligible)
    # Same design can be stored in a different row order when the lookup agrees.
    changed <- pair
    for (i in 1:2) {
      changed[[i]]$config$population_spec$design_matrix <- cbind(1, seq(-1, 1, length.out = 8))
      changed[[i]]$config$population_spec$design_columns <- c("(Intercept)", "Experience")
    }
    changed[[2]]$config$population_spec$design_matrix <-
      changed[[2]]$config$population_spec$design_matrix[8:1, 2:1]
    changed[[2]]$config$population_spec$person_lookup <- 8:1
    changed[[2]]$config$population_spec$design_columns <- c("Experience", "(Intercept)")
    expect_true(audit(changed)$eligible)
    changed <- pair
    changed[[2]]$config$facet_interactions <- "Rater:Criterion"
    expect_match(audit(changed)$reason, "identical facet constraints")
    changed <- pair
    changed[[2]]$config$step_specs[[1]]$target <- .5
    expect_false(audit(changed)$eligible)
    changed <- pair
    changed[[2]]$config$gpcm_spec$n_params <- 99L
    expect_match(audit(changed)$reason, "G-1")
    changed <- pair
    changed[[2]]$config$gpcm_spec$identification <- "unknown"
    expect_false(audit(changed)$eligible)
    changed <- pair
    changed[[1]]$config$method <- "JML"
    expect_match(audit(changed)$reason, "requires MML")
    changed <- pair
    changed[[2]]$config$step_facet <- setdiff(c("Rater", "Criterion"), owner)
    expect_identical(audit(changed)$relation, "PCM_GPCM_owner_mismatch")
  }
})

test_that("unit slopes embed the PCM marginal objective and nuisance gradient in GPCM", {
  for (owner in c("Criterion", "Rater")) {
    pair <- gpcm_lrt_fixture(owner)
    for (i in 1:2) {
      pair[[i]]$config$population_spec$design_matrix <- cbind(1, seq(-1, 1, length.out = 8))
      pair[[i]]$config$population_spec$design_columns <- c("(Intercept)", "Experience")
    }
    sizes <- lapply(pair, function(f) mfrmr:::build_param_sizes(f$config))
    evaluators <- lapply(pair, function(f) {
      idx <- mfrmr:::build_indices(f$prep, f$config$step_facet, f$config$slope_facet,
                                  f$config$interaction_specs)
      size <- mfrmr:::build_param_sizes(f$config)
      cache <- mfrmr:::make_param_cache(size, f$config, idx, is_mml = TRUE)
      mfrmr:::make_mfrm_direct_evaluator("MML", cache, idx, f$config, size,
                                        mfrmr:::gauss_hermite_normal(31L))
    })
    for (offset in c(-.2, 0, .3)) {
      p <- seq(-.4, .4, length.out = sum(unlist(sizes[[1]]))) + offset
      parts <- mfrmr:::split_params(p, sizes[[1]])
      parts$log_slopes <- rep(0, sizes[[2]]$log_slopes)
      g <- unlist(parts[names(sizes[[2]])], use.names = FALSE)
      expect_equal(evaluators[[1]]$value(p), evaluators[[2]]$value(g), tolerance = 1e-9)
      g_gradient <- mfrmr:::split_params(evaluators[[2]]$gradient(g), sizes[[2]])
      expect_equal(evaluators[[1]]$gradient(p),
                   unlist(g_gradient[names(sizes[[1]])], use.names = FALSE), tolerance = 1e-8)
    }
  }
})

test_that("PCM/GPCM LRT uses verified degrees of freedom without granting slope intervals", {
  eligible <- TRUE
  caution <- NULL
  testthat::local_mocked_bindings(
    mfrm_convergence_state = function(...) list(code_converged = TRUE, inference_ready = FALSE, severity = "pass"),
    mfrm_ic_fit_check = function(...) list(eligible = eligible, basis = "test", review = "test solution failure", caution=caution),
    mfrmr_get_readiness_record = function(...) list(fit = data.frame(ReasonCodes = "slope_intervals_unavailable")),
    mfrm_fit_decision_summary = function(...) data.frame(Why = "Slope intervals unavailable"),
    .package = "mfrmr"
  )
  pair <- gpcm_lrt_fixture()
  compare <- function(x = pair, nested = TRUE) compare_mfrm(x[[1]], x[[2]], labels = c("PCM", "GPCM"), nested = nested)
  before <- pair
  result <- compare()
  expect_identical(result$comparison_basis$lrt_status, "computed")
  expect_equal(result$lrt$ChiSq, 6)
  expect_equal(result$lrt$df, 3)
  expect_equal(result$lrt$p_value, pchisq(6, 3, lower.tail = FALSE))
  expect_identical(result$lrt$ReferenceDistribution, "asymptotic_chi_square")
  expect_false(any(result$table$InferenceReady))
  expect_identical(pair, before)
  caution <- "Ill-conditioned information: numerical refinement passed; review interval width."
  expect_warning(cautioned <- compare(),"Ill-conditioned")
  expect_identical(cautioned$table$ICFitCaution,rep(caution,2))
  expect_match(cautioned$comparison_basis$lrt_reason,"Ill-conditioned")
  expect_output(print(cautioned),"Caution: Ill-conditioned")
  expect_equal(cautioned$lrt,result$lrt)
  caution <- NULL
  expect_null(compare(nested = FALSE)$lrt)
  for (invalid in list(NA, "TRUE", c(TRUE, FALSE), 1)) {
    expect_error(compare(nested = invalid), "`nested` must be a single logical value", fixed = TRUE)
  }
  reversed <- compare_mfrm(pair[[2]], pair[[1]], labels = c("GPCM", "PCM"), nested = TRUE)
  expect_equal(result$lrt, reversed$lrt)
  repeated_label <- compare_mfrm(pair[[1]], pair[[2]], labels = c("Model", "Model"), nested = TRUE)
  expect_equal(repeated_label$lrt$p_value, result$lrt$p_value)
  expect_output(print(result), "Reference: asymptotic chi-square", fixed = TRUE)
  tiny_p <- result
  tiny_p$lrt$p_value <- 1e-10
  expect_output(print(tiny_p), "p <")
  expect_false(grepl("p = 0[.]0000|asymptotic_chi_square", paste(capture.output(print(tiny_p)), collapse = " ")))
  eligible <- FALSE
  expect_warning(refused <- compare(), "LRT was not computed")
  expect_null(refused$lrt)
  eligible <- TRUE
  # A negative likelihood gain must remain a numerical failure, not p = 1.
  pair[[2]]$opt$value <- 101
  pair[[2]]$summary <- mfrmr:::build_mfrm_ic_contract(-101, length(pair[[2]]$opt$par), pair[[2]]$prep, pair[[2]]$config, "MML")
  expect_warning(refused <- compare(), "negative likelihood-ratio")
  expect_null(refused$lrt)
  pair <- before
  pair[[2]]$opt$par <- c(pair[[2]]$opt$par, 0)
  pair[[2]]$summary <- mfrmr:::build_mfrm_ic_contract(-97, length(pair[[2]]$opt$par), pair[[2]]$prep, pair[[2]]$config, "MML")
  expect_warning(refused <- compare(), "number of null restrictions")
  expect_null(refused$lrt)
})

test_that("weighting reviews explain a requested but unavailable test", {
  testthat::local_mocked_bindings(
    mfrmr_get_readiness_record = function(...) list(fit = data.frame(NumericalState = "ready")),
    .package = "mfrmr"
  )
  comparison <- list(table = data.frame(LogLik = c(-100, -97)),
    comparison_basis = list(same_data = TRUE, nested_requested = TRUE,
                            lrt_reason = "Population models differ."))
  fit <- function(model) list(config = list(model = model, method = "MML"))
  contract <- mfrmr:::.weighting_review_comparison_contract(fit("PCM"), fit("GPCM"),
    comparison, "PCM", TRUE, "not_computed")
  expect_match(contract$RecommendedUse, "Population models differ", fixed = TRUE)
  expect_false(grepl("set nested = TRUE", contract$RecommendedUse, fixed = TRUE))
})


test_that("separate-owner PCM nesting counts slopes rather than steps", {
  for (owner in c("Rater","Criterion")) {
    slope <- setdiff(c("Rater","Criterion"),owner)
    pair <- gpcm_lrt_fixture(owner,slope)
    audit <- audit_compare_mfrm_nesting(pair,c("PCM","GPCM"))
    expect_true(audit$eligible)
    expect_equal(audit$df,length(pair[[2]]$config$facet_levels[[slope]])-1L)
    expect_match(audit$reason,paste0("slopes on facet '",slope,"'"),fixed=TRUE)
    pair[[2]]$config$gpcm_spec$slope_facet <- owner
    expect_false(audit_compare_mfrm_nesting(pair,c("PCM","GPCM"))$eligible)
  }
})
