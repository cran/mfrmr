test_that("anchor observation counts keep group order and weighted-sum precision", {
  counts <- .mfrmr_anchor_observation_counts(
    c("z", "a", "a", "a", NA_character_, "é"), c(.1, 1e16, 1, 1, .5, .2))
  expect_identical(counts$level, c("a", "z", "é", NA_character_))
  expect_identical(counts$raw, c(3L, 1L, 1L, 1L))
  # R's sum uses extended precision only where the platform provides it.
  expect_identical(counts$weighted, c(sum(c(1e16, 1, 1)), .1, .2, .5))
  expect_identical(.mfrmr_anchor_observation_counts(character(), numeric()),
    list(level = character(), raw = integer(), weighted = numeric()))
})

test_that("empty anchors preserve the full audit and weighted design checks", {
  ratings <- data.frame(Person = c("p1", "p1", "p2", "p2"),
    Rater = c("r1", "r2", "r1", "r2"), Score = c(0, 1, 1, 0),
    Weight = c(1, 2, 3, 4))
  prep <- prepare_mfrm_data(ratings, "Person", "Rater", "Score", weight_col = "Weight")
  # These rows are discarded by normalization. They exercise the general
  # audit path, independently of the empty-input shortcut.
  missing_anchor <- data.frame(Facet = NA_character_, Level = NA_character_,
    Anchor = NA_real_)
  missing_group <- data.frame(Facet = NA_character_, Level = NA_character_,
    Group = NA_character_, GroupValue = NA_real_)
  for (scope in c("facet", "remainder", "anchor_origin")) {
    args <- list(prep = prep, centering_scope = scope, dummy_facets = "Rater",
      min_common_anchors = 2L, min_obs_per_element = 5, min_obs_per_category = 6)
    actual <- do.call(audit_anchor_tables, args)
    reference <- do.call(audit_anchor_tables, c(args,
      list(anchor_df = missing_anchor, group_anchor_df = missing_group)))
    expect_identical(actual, reference, info = scope)
    expect_identical(do.call(audit_anchor_tables, c(args,
      list(anchor_df = data.frame(), group_anchor_df = data.frame()))), reference)
    expect_identical(unname(actual$issue_counts$N), rep(0L, 13L))
    expect_identical(actual$facet_summary$FreeLevels, c(2L, 2L))
    expect_equal(actual$design_checks$category_counts$WeightedN, c(5, 5))
    low <- actual$design_checks$low_observation_levels
    expect_equal(low$WeightedN, c(3, 4))
    expect_identical(low$Level, c("p1", "r1"))
    actual$issues$unknown_anchor_facets$new_column <- logical()
    expect_identical(do.call(audit_anchor_tables, args), reference)
  }
})

local({
  ratings <- load_mfrmr_data("example_core")
  anchor <- data.frame(Facet = "Rater", Level = "R01", Anchor = .7)
  args <- list(data = ratings, person = "Person", facets = c("Rater", "Criterion"),
               score = "Score", anchors = anchor)
  whole <- suppressWarnings(do.call(fit_mfrm, args))
  legacy <- suppressWarnings(do.call(fit_mfrm, c(args, list(centering_scope = "remainder"))))

  test_that("public centering fixes the declared value and preserves legacy reproduction", {
    for (fit in list(whole, legacy)) {
      r <- subset(fit$facets$others, Facet == "Rater")
      expect_equal(unname(r$Estimate[r$Level == "R01"]), .7, tolerance = 1e-12)
      target <- if (fit$config$centering_scope == "facet") -.7 else 0
      expect_equal(sum(r$Estimate[r$Level != "R01"]), target, tolerance = 1e-12)
      expect_identical(fit$config$replay_inputs$centering_scope, fit$config$centering_scope)
      expect_identical(fit$config$anchor_review$centering_scope, fit$config$centering_scope)
      expect_false(isTRUE(fit$config$population_spec$active))
    }
    expect_match(paste(whole$config$anchor_review$recommendations, collapse = " "),
                 "including anchors and group totals")
    expect_match(paste(legacy$config$anchor_review$recommendations, collapse = " "),
                 "Legacy remainder centering")
    nesting <- audit_compare_mfrm_nesting(list(whole, legacy), c("whole", "legacy"))
    expect_false(nesting$eligible)
    expect_match(nesting$reason, "centering_scope")
  })

  test_that("PCM and JML use the same anchored constraint equations", {
    for (model in c("RSM", "PCM")) for (method in c("MML", "JML")) {
      if (model == "RSM" && method == "MML") next
      extra <- list(model = model, method = method)
      if (model == "PCM") extra$step_facet <- "Criterion"
      fit <- suppressWarnings(do.call(fit_mfrm, c(args, extra)))
      r <- subset(fit$facets$others, Facet == "Rater")
      expect_equal(unname(r$Estimate[r$Level == "R01"]), .7, tolerance = 1e-12)
      expect_equal(sum(r$Estimate), 0, tolerance = 1e-12)
      expect_true(is.finite(fit$opt$value))
      expect_identical(fit$config$centering_scope, "facet")
    }
  })

  test_that("refit and export settings preserve saved and metadata-free conventions", {
    old <- legacy
    old$config$centering_scope <- old$config$replay_inputs$centering_scope <- NULL
    for (fit in list(whole, legacy, old)) {
      expected <- fit$config$centering_scope %||% "remainder"
      file <- tempfile(fileext = ".rds")
      saveRDS(fit, file)
      saved <- readRDS(file)
      unlink(file)
      expect_identical(saved$opt$par, fit$opt$par)
      expect_identical(saved$config$facet_specs, fit$config$facet_specs)
      refit <- mfrmr_gqs_refit_arguments(saved, ratings, 31L)
      expect_identical(refit$centering_scope, expected)
      expect_identical(resolve_dff_refit_controls(saved)$centering_scope, expected)
      code <- build_replay_fit_mfrm_lines(saved$config$replay_inputs,
        saved$population, NULL, saved$config$source_columns, saved$config)
      expect_match(paste(code, collapse = "\n"), paste0('centering_scope = "', expected, '"'))
    }
    refit <- suppressWarnings(do.call(fit_mfrm, mfrmr_gqs_refit_arguments(old, ratings, 31L)))
    expect_equal(refit$opt$par, legacy$opt$par, tolerance = 1e-10)
  })

  test_that("review and fitting reject contradictory totals before optimization", {
    fixed <- data.frame(Facet = "Rater", Level = paste0("R0", 1:4), Anchor = .2)
    review_args <- args
    review_args$anchors <- fixed
    expect_error(do.call(review_mfrm_anchors, review_args), "conflict with whole-facet")
    review_args$anchor_policy <- "silent"
    expect_error(do.call(fit_mfrm, review_args), "conflict with whole-facet")
    expect_error(do.call(fit_mfrm, c(args, list(centering_scope = "unknown"))), "arg")
    # A group mean fixes a total, not every member; the remaining level balances it.
    group <- data.frame(Facet = "Rater", Level = paste0("R0", 1:3), Group = "g", GroupValue = .1)
    review <- review_mfrm_anchors(ratings, "Person", c("Rater", "Criterion"), "Score",
      anchors = anchor, group_anchors = group)
    specs <- prepare_constraint_specs(whole$prep, anchor_review = review, centering_scope = "facet")
    b <- expand_facet_with_constraints(.05, specs$facet_specs$Rater)
    expect_equal(mean(b[1:3]), .1, tolerance = 1e-12)
    expect_equal(unname(b[4]), -.3, tolerance = 1e-12)
    expect_equal(unname(b[1]), .7, tolerance = 1e-12)
  })

  test_that("unanchored equations and population assumptions do not change", {
    a <- prepare_constraint_specs(whole$prep, centering_scope = "facet")
    b <- prepare_constraint_specs(whole$prep, centering_scope = "remainder")
    for (name in names(a$facet_specs)) {
      sa <- a$facet_specs[[name]]; sb <- b$facet_specs[[name]]
      free <- seq_len(sa$n_params) / 7
      expect_identical(expand_facet_with_constraints(free, sa),
                       expand_facet_with_constraints(free, sb))
      expect_identical(constraint_jacobian(sa), constraint_jacobian(sb))
    }
  })

  test_that("anchor export preserves values instead of rounding away a constraint", {
    tab <- make_anchor_table(whole, facets = "Rater", readiness_policy = "review")
    r <- subset(whole$facets$others, Facet == "Rater")
    expect_identical(tab$Anchor, unname(r$Estimate[match(tab$Level, r$Level)]))
    expect_equal(sum(tab$Anchor), 0, tolerance = 1e-12)
    expect_s3_class(review_mfrm_anchors(ratings, "Person", c("Rater", "Criterion"),
      "Score", anchors = tab), "mfrm_anchor_review")
    rounded <- make_anchor_table(whole, digits = 6, facets = "Rater", readiness_policy = "review")
    expect_identical(rounded$Anchor, round(tab$Anchor, 6))
    expect_error(make_anchor_table(whole, digits = NA_real_, readiness_policy = "review"), "digits")
  })
})

test_that("compatibility wrappers forward explicit and default centering", {
  captured <- NULL
  local_mocked_bindings(fit_mfrm = function(...) {
    captured <<- list(...)
    stop("captured fitting arguments")
  })
  ratings <- load_mfrmr_data("example_core")
  expect_error(run_mfrm_facets(ratings, "Person", c("Rater", "Criterion"), "Score",
    centering_scope = "remainder"), "captured fitting arguments")
  expect_identical(captured$centering_scope, "remainder")
  expect_error(mfrmRFacets(ratings, "Person", c("Rater", "Criterion"), "Score"),
    "captured fitting arguments")
  expect_identical(captured$centering_scope, "facet")
})
