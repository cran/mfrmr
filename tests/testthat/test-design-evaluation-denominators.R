design_evaluation_denominator_fixture <- local({
  planning_schema <- NULL
  function() {
    # Public evaluation shape with known outcomes; no simulated responses or fits.
    prototype <- mfrmr:::design_eval_result_prototype()
    available <- prototype[rep(NA_integer_, 2), ]
    for (name in names(available)) {
      available[[name]] <- if (is.character(prototype[[name]])) rep("", 2) else if (
        is.logical(prototype[[name]])) rep(FALSE, 2) else rep(0, 2)
    }
    available$design_id <- "D1"
    available$rep <- 1L
    available$Facet <- c("Rater", "Criterion")
    available$n_person <- 30L
    available$n_rater <- 3L
    available$n_criterion <- 3L
    available$raters_per_person <- 2L
    available$Converged <- TRUE
    available$RecoveryComparable <- TRUE
    available$Separation <- 3
    available$Reliability <- .9
    available$SeverityRMSE <- .1
    available$MinCategoryCount <- 20
    complete <- dplyr::bind_rows(lapply(1:10, function(rep) {
      row <- available
      row$design_id <- "D2"
      row$n_person <- 40L
      row$rep <- rep
      row
    }))
    runs <- data.frame(
      design_id = rep(c("D1", "D2", "D3"), each = 10),
      rep = rep(1:10, 3),
      RunOK = c(TRUE, rep(FALSE, 9), rep(TRUE, 10), rep(FALSE, 10)),
      Converged = c(TRUE, rep(FALSE, 9), rep(TRUE, 10), rep(FALSE, 10)),
      RaterComponents = 1L, CriterionComponents = 1L,
      ElapsedSec = 0
    )
    out <- structure(list(results = dplyr::bind_rows(available, complete),
      rep_overview = runs), class = "mfrm_design_evaluation")
    # Real evaluation objects already carry this deterministic design metadata.
    # Reuse only that unchanged input; summaries and decisions still run afresh.
    if (is.null(planning_schema))
      planning_schema <<- mfrmr:::simulation_object_planning_schema(out)
    out$planning_schema <- planning_schema
    out
  }
})

test_that("design recommendations retain failed runs in convergence denominators", {
  x <- design_evaluation_denominator_fixture()
  x$planning_schema <- NULL # Exercise reconstruction for legacy saved results.
  runs <- x$rep_overview
  build_schema <- mfrmr:::simulation_planning_schema
  schema_calls <- 0L
  local_mocked_bindings(simulation_planning_schema = function(...) {
    schema_calls <<- schema_calls + 1L
    build_schema(...)
  }, .package = "mfrmr")

  s <- summary(x, digits = 10)
  expect_identical(schema_calls, 1L)
  partial <- s$design_summary[s$design_summary$design_id == "D1", ]
  full <- s$design_summary[s$design_summary$design_id == "D2", ]
  expect_equal(partial$Reps, c(10, 10))
  expect_equal(partial$AvailableReps, c(1, 1))
  expect_equal(partial$ConvergenceRate, c(.1, .1))
  expect_equal(partial$McseConvergenceRate, rep(sqrt(.1 * .9 / 10), 2), tolerance = 1e-9)
  expect_equal(full$Reps, c(10, 10))
  expect_equal(full$AvailableReps, c(10, 10))
  expect_equal(full$ConvergenceRate, c(1, 1))
  expect_equal(s$design_summary$MeanSeverityRMSE, rep(.1, 4))
  expect_equal(s$overview$Replications, 30)
  expect_equal(s$overview$ConvergedRuns, 11)
  p <- plot(x, facet = "Rater", metric = "convergencerate", draw = FALSE)
  expect_equal(p$data$y, c(.1, 1))
  expect_identical(recommend_mfrm_design(x)$recommended$design_id, "D2")
  expect_identical(recommend_mfrm_design(s)$recommended$design_id, "D2")

  x$rep_overview$Converged[2] <- NA
  expect_equal(summary(x)$design_summary$ConvergenceRate, c(.1, 1, .1, 1))
  x$results <- mfrmr:::design_eval_result_prototype()
  x$rep_overview <- runs[runs$design_id == "D3", ]
  failed <- summary(x)
  expect_equal(failed$overview$Replications, 10)
  expect_equal(failed$overview$ConvergedRuns, 0)
  expect_equal(nrow(failed$design_summary), 0)
  expect_error(recommend_mfrm_design(x), "No design summary rows")
})

test_that("recommendations do not rebuild unrelated planning appendices", {
  x <- design_evaluation_denominator_fixture()
  expected <- recommend_mfrm_design(summary(x), facets = c("Rater", "Criterion"))
  local_mocked_bindings(
    simulation_compact_structural_design_review_summary = function(...)
      stop("Unexpected planning appendix"),
    .package = "mfrmr"
  )
  actual <- recommend_mfrm_design(x, facets = c("Rater", "Criterion"))
  decisions <- c("facet_table", "design_table", "recommended", "thresholds")
  expect_identical(actual[decisions], expected[decisions])
})

test_that("design decisions preserve full precision through summaries and saved results", {
  boundaries <- c(Separation = 1.9999996, Reliability = .7999996,
                  SeverityRMSE = .5000004, MisfitRate = .1000004)
  metrics <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) names(boundaries) else c("Reliability", "SeverityRMSE")
  for (metric in metrics) {
    x <- design_evaluation_denominator_fixture()
    x$results[[metric]][x$results$design_id == "D2"] <- boundaries[[metric]]
    expect_equal(nrow(recommend_mfrm_design(x)$recommended), 0L, info = metric)
    precisions <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) c(0L, 3L, 10L) else 0L
    for (digits in precisions) {
      s <- summary(x, digits = digits)
      column <- paste0("Mean", metric)
      expect_equal(s$design_summary[[column]][s$design_summary$design_id == "D2"],
                   rep(boundaries[[metric]], 2), tolerance = 1e-14)
      expect_equal(nrow(recommend_mfrm_design(s)$recommended), 0L, info = metric)
    }
    saved <- tempfile(fileext = ".rds")
    saveRDS(s, saved)
    expect_identical(recommend_mfrm_design(readRDS(saved)), recommend_mfrm_design(s))
    unlink(saved)
  }
  x <- design_evaluation_denominator_fixture()
  x$results$Reliability[x$results$design_id == "D2"] <- .8
  expect_identical(recommend_mfrm_design(x)$recommended$design_id, "D2")
  s <- summary(x, digits = 0)
  expect_identical(recommend_mfrm_design(s)$recommended$design_id, "D2")
  s$summary_precision <- NULL
  expect_error(recommend_mfrm_design(s), "Rebuild it with summary\\(original_evaluation\\)")
})

test_that("workload limits include failed runs and preserve ranking within limits", {
  x <- design_evaluation_denominator_fixture()
  # D1 has a costly failed run; its mean workload would conceal the overrun.
  x$rep_overview$Observations <- c(rep(180, 9), 480, rep(240, 10), rep(300, 10))
  x$rep_overview$MaxRatingsPerRater <- c(rep(90, 9), 200, rep(80, 10), rep(100, 10))
  s <- summary(x)
  d1 <- s$design_summary[s$design_summary$design_id == "D1", ]
  expect_equal(d1$MaxRatings, c(480, 480))
  expect_equal(d1$MaxRatingsPerRater, c(200, 200))
  expect_identical(recommend_mfrm_design(x, min_convergence_rate = .1)$recommended$design_id, "D1")

  total <- recommend_mfrm_design(x, min_convergence_rate = .1, max_ratings = 300)
  expect_identical(total$recommended$design_id, "D2")
  expect_false(total$design_table$RatingsPass[total$design_table$design_id == "D1"])
  expect_true(all(total$facet_table$Pass))
  per_rater <- recommend_mfrm_design(s, min_convergence_rate = .1, max_ratings_per_rater = 100)
  expect_identical(per_rater$recommended$design_id, "D2")
  expect_false(per_rater$design_table$RaterWorkloadPass[per_rater$design_table$design_id == "D1"])
  boundary <- recommend_mfrm_design(x, min_convergence_rate = .1,
    max_ratings = 480, max_ratings_per_rater = 200)
  expect_identical(boundary$recommended$design_id, "D1")
  expect_equal(boundary$thresholds$max_ratings, 480)
  expect_equal(boundary$thresholds$max_ratings_per_rater, 200)
  expect_equal(nrow(recommend_mfrm_design(x, max_ratings = 0)$recommended), 0)
  expect_equal(nrow(recommend_mfrm_design(x, max_ratings_per_rater = 79)$recommended), 0)

  # Unknown counts cannot be omitted when checking a hard workload limit.
  x$rep_overview$Observations[10] <- NA
  x$rep_overview$MaxRatingsPerRater[10] <- NA
  d1 <- summary(x)$design_summary
  expect_true(all(is.na(d1$MaxRatings[d1$design_id == "D1"])))
  expect_true(all(is.na(d1$MaxRatingsPerRater[d1$design_id == "D1"])))
  expect_identical(recommend_mfrm_design(x, min_convergence_rate = .1,
    max_ratings = 1000, max_ratings_per_rater = 1000)$recommended$design_id, "D2")
})

test_that("a design must provide every requested facet before it can pass", {
  x <- design_evaluation_denominator_fixture()
  x$results <- x$results[!(x$results$design_id == "D1" & x$results$Facet == "Criterion"), ]
  for (object in list(x, summary(x))) {
    rec <- recommend_mfrm_design(object, facets = c("Rater", "Criterion"),
      min_convergence_rate = .1)
    partial <- rec$design_table[rec$design_table$design_id == "D1", ]
    complete <- rec$design_table[rec$design_table$design_id == "D2", ]
    expect_identical(rec$recommended$design_id, "D2")
    expect_identical(partial$FacetsMissing, "Criterion")
    expect_equal(partial$FacetsRequired, 2)
    expect_equal(partial$FacetsPassing, 1)
    expect_false(partial$Pass)
    expect_identical(complete$FacetsMissing, "")
    expect_true(complete$Pass)
    expect_identical(recommend_mfrm_design(object, min_convergence_rate = .1)$recommended$design_id, "D2")
    # An intentional single-facet request remains supported.
    expect_identical(recommend_mfrm_design(object, facets = "Rater",
      min_convergence_rate = .1)$recommended$design_id, "D1")
  }
})

test_that("default facet requests retain stored custom names when results are missing", {
  x <- design_evaluation_denominator_fixture()
  x$planning_schema <- NULL # Custom facet labels need their own schema.
  x$settings <- list(facet_names = c(rater = "Judge", criterion = "Task"))
  x$results$Facet <- ifelse(x$results$Facet == "Rater", "Judge", "Task")
  x$results <- x$results[x$results$Facet == "Judge", ]
  for (object in list(x, summary(x))) {
    expect_error(recommend_mfrm_design(object), "Requested facets not found.*Task")
    expect_identical(recommend_mfrm_design(object, facets = "Judge")$recommended$design_id, "D2")
    for (facets in list(character(0), NA_character_, "", " ")) {
      expect_error(recommend_mfrm_design(object, facets = facets), "`facets` must contain")
    }
  }
})

test_that("legacy workload records are unknown without changing uncapped recommendations", {
  x <- design_evaluation_denominator_fixture()
  x$rep_overview$Observations <- 240
  expect_identical(recommend_mfrm_design(x, max_ratings = 240)$recommended$design_id, "D2")
  expect_equal(nrow(recommend_mfrm_design(x, max_ratings_per_rater = 1000)$recommended), 0)
  s <- summary(x)
  s$design_summary$MaxRatings <- NULL
  s$design_summary$MaxRatingsPerRater <- NULL
  expect_identical(recommend_mfrm_design(s)$recommended$design_id, "D2")
  expect_equal(nrow(recommend_mfrm_design(s, max_ratings = 1000)$recommended), 0)
  expect_equal(nrow(recommend_mfrm_design(s, max_ratings_per_rater = 1000)$recommended), 0)
})

test_that("workload limits require finite non-negative whole numbers", {
  x <- design_evaluation_denominator_fixture()
  for (name in c("max_ratings", "max_ratings_per_rater")) {
    for (value in list(-1, .5, NA_real_, Inf, numeric(0), c(1, 2), "100", TRUE)) {
      args <- c(list(x = x), stats::setNames(list(value), name))
      expect_error(do.call(recommend_mfrm_design, args), paste0("`", name, "` must be"))
    }
  }
})

test_that("actual unweighted assignment counts survive fit and diagnostic failures", {
  calls <- 0L
  testthat::local_mocked_bindings(
    fit_mfrm = function(...) {
      calls <<- calls + 1L
      if (calls == 1L) stop("forced fit failure")
      list()
    },
    diagnose_mfrm = function(...) stop("forced diagnostic failure"),
    .package = "mfrmr"
  )
  skeleton <- rbind(
    expand.grid(TemplatePerson = paste0("P", 1:4), Judge = "J1", Task = c("T1", "T2")),
    expand.grid(TemplatePerson = paste0("P", 1:2), Judge = "J2", Task = c("T1", "T2")),
    data.frame(TemplatePerson = "P3", Judge = "J3", Task = "T1")
  )
  skeleton$TemplatePersonReuse <- TRUE
  skeleton$Weight <- 7
  spec <- build_mfrm_sim_spec(n_person = 4, n_rater = 3,
    n_criterion = 2, raters_per_person = 2, facet_names = c("Judge", "Task"),
    assignment = "skeleton", design_skeleton = skeleton)
  x <- evaluate_mfrm_design(sim_spec = spec, n_person = 4, n_rater = 3,
    n_criterion = 2, raters_per_person = 2, reps = 2, seed = 123, progress = FALSE)
  expect_equal(x$rep_overview$Observations, c(13, 13))
  expect_equal(x$rep_overview$MaxRatingsPerRater, c(8, 8))
  expect_equal(x$rep_overview$RaterComponents, c(1L, 1L))
  expect_equal(x$rep_overview$CriterionComponents, c(1L, 1L))
  expect_identical(x$rep_overview$ErrorComponent, c("fit", "diagnostics"))
  expect_false(any(x$rep_overview$RunOK))
  expect_equal(nrow(x$results), 0)
})

test_that("connectivity and sparse overlap screens retain failed-run evidence", {
  x <- design_evaluation_denominator_fixture()
  x$rep_overview$SparseDesignActive <- TRUE
  x$rep_overview$MinCommonPersonsPerRaterPair <- 3L
  x$rep_overview$ZeroCommonRaterPairs <- 0L
  x$rep_overview$RaterPairsBelowTarget <- 0L
  x$rep_overview$TargetCommonPersonsPerRaterPair <- 2L
  # Only a failed D1 replication falls below the requested overlap target.
  x$rep_overview$MinCommonPersonsPerRaterPair[10] <- 1L
  x$rep_overview$RaterPairsBelowTarget[10] <- 1L
  overlap <- recommend_mfrm_design(x, min_convergence_rate = .1)
  expect_identical(overlap$recommended$design_id, "D1")
  expect_identical(overlap$recommended$ConnectivityStatus, "connected")
  expect_identical(overlap$recommended$LinkReviewStatus, "review")
  expect_match(overlap$recommended$LinkReviewReason, "below the requested")

  x$rep_overview$RaterComponents[10] <- 2L
  for (object in list(x, summary(x))) {
    rec <- recommend_mfrm_design(object, min_convergence_rate = .1)
    partial <- rec$design_table[rec$design_table$design_id == "D1", ]
    expect_equal(partial$MaxRaterComponents, 2)
    expect_equal(partial$MaxCriterionComponents, 1)
    expect_equal(partial$DisconnectedReps, 1)
    expect_identical(partial$ConnectivityStatus, "disconnected")
    expect_false(partial$ConnectivityPass)
    expect_identical(rec$recommended$design_id, "D2")
    explicit <- recommend_mfrm_design(object, min_convergence_rate = .1,
      require_connected = FALSE)
    expect_identical(explicit$recommended$design_id, "D1")
    expect_identical(explicit$recommended$ConnectivityStatus, "disconnected")
    expect_false(explicit$thresholds$require_connected)
    expect_match(explicit$caveats$connectivity, "explicitly disabled")
  }
  # Unknown records must not hide a disconnection already observed elsewhere.
  x$rep_overview$RaterComponents[9] <- NA_integer_
  rec <- recommend_mfrm_design(x, min_convergence_rate = .1)
  expect_identical(rec$design_table$ConnectivityStatus[rec$design_table$design_id == "D1"], "disconnected")
  x$rep_overview$RaterComponents <- 1L
  x$rep_overview$CriterionComponents[10] <- 2L
  rec <- recommend_mfrm_design(x, facets = "Rater", min_convergence_rate = .1)
  expect_identical(rec$recommended$design_id, "D2")
  expect_identical(rec$design_table$ConnectivityStatus[rec$design_table$design_id == "D1"], "disconnected")
})

test_that("missing connectivity records cannot pass the default screen", {
  x <- design_evaluation_denominator_fixture()
  x$rep_overview$RaterComponents <- NULL
  x$rep_overview$CriterionComponents <- NULL
  old_summary <- summary(x)
  old_summary$design_summary$MaxRaterComponents <- NULL
  old_summary$design_summary$MaxCriterionComponents <- NULL
  old_summary$design_summary$DisconnectedReps <- NULL
  for (object in list(x, summary(x), old_summary)) {
    rec <- recommend_mfrm_design(object)
    expect_equal(nrow(rec$recommended), 0)
    expect_true(all(rec$design_table$ConnectivityStatus == "not_assessed"))
    expect_identical(recommend_mfrm_design(object, require_connected = FALSE)$recommended$design_id, "D2")
  }
  for (value in list(NULL, NA, 1, "TRUE", c(TRUE, FALSE))) {
    expect_error(recommend_mfrm_design(x, require_connected = value), "`require_connected` must be")
  }
})

test_that("assignment topology reaches recommendations independently of fit metrics", {
  # Detailed propagation remains in the extended local/CI suite.
  skip_on_cran()
  captured <- NULL
  # Hold fit quality constant to isolate the public design/recommendation path.
  testthat::local_mocked_bindings(
    fit_mfrm = function(data, ...) { captured <<- data; list() },
    mfrm_inference_ready = function(...) TRUE,
    diagnose_mfrm = function(...) list(
      reliability = data.frame(Facet = c("Judge", "Task"), Levels = c(4, 2),
        Separation = 3, Strata = 4, Reliability = .9, MeanInfit = 1, MeanOutfit = 1),
      fit = data.frame(Facet = c("Judge", "Task"), InfitZSTD = 0, OutfitZSTD = 0),
      measures = data.frame(Facet = rep(c("Judge", "Task"), c(4, 2)),
        Level = c(paste0("J", 1:4), "T1", "T2"), Estimate = 0)
    ),
    .package = "mfrmr"
  )
  assignments <- list(
    disconnected = c("J1", "J2", "J1", "J2", "J3", "J4", "J3", "J4"),
    cycle = c("J1", "J2", "J2", "J3", "J3", "J4", "J4", "J1")
  )
  for (name in names(assignments)) {
    roster <- data.frame(TemplatePerson = rep(paste0("P", 1:4), each = 2),
      Judge = assignments[[name]])
    skeleton <- merge(roster, data.frame(Task = c("T1", "T2")), by = NULL)
    skeleton$TemplatePersonReuse <- TRUE
    spec <- build_mfrm_sim_spec(n_person = 4, n_rater = 4, n_criterion = 2,
      raters_per_person = 2, rater_sd = 0, criterion_sd = 0,
      facet_names = c("Judge", "Task"), assignment = "skeleton", design_skeleton = skeleton)
    expect_warning(x <- evaluate_mfrm_design(sim_spec = spec, n_person = 4, n_rater = 4,
      n_criterion = 2, raters_per_person = 2, reps = 1, seed = 123, progress = FALSE), NA)
    review <- describe_mfrm_data(captured, person = "Person", facets = c("Judge", "Task"),
      score = "Score", rating_min = 1, rating_max = 4, keep_original = TRUE,
      include_agreement = FALSE)
    expect_equal(c(x$rep_overview$RaterComponents, x$rep_overview$CriterionComponents),
      review$design_connectivity$Components)
    expect_equal(x$rep_overview$Observations, 16)
    expect_equal(x$rep_overview$MaxRatingsPerRater, 4)
    rec <- recommend_mfrm_design(x)
    expect_true(all(rec$facet_table$Pass))
    expect_identical(rec$design_table$ConnectivityStatus,
      if (name == "cycle") "connected" else "disconnected")
    expect_equal(nrow(rec$recommended), if (name == "cycle") 1L else 0L)
    if (name == "cycle") {
      # J1/J3 and J2/J4 have no direct overlap, but indirect links connect all.
      pairs <- data.frame(Person = captured$Person, Rater = captured$Judge)
      links <- mfrmr:::simulation_sparse_rater_pair_table(pairs, paste0("J", 1:4))
      expect_equal(sum(links$CommonPersons == 0), 2)
    }
  }
})
