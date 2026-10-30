halo_person_fixture <- function(tasks = 1L, persons = 6L) {
  d <- expand.grid(Person = paste0("p", seq_len(persons)),
    Rater = c("r1", "r2"), Criterion = c("c1", "c2"),
    Task = paste0("t", seq_len(tasks)), stringsAsFactors = FALSE)
  profiles <- list(0:5, c(0, 2, 1, 5, 3, 4), c(4, 0, 5, 1, 3, 2), c(2, 5, 0, 4, 1, 3))
  node <- match(paste(d$Rater, d$Criterion), c("r1 c1", "r2 c1", "r1 c2", "r2 c2"))
  d$Observed <- vapply(seq_len(nrow(d)), function(i)
    profiles[[node[i]]][match(d$Person[i], paste0("p", 1:6))], numeric(1))
  list(fit = structure(list(config = list(facet_names = c("Rater", "Criterion", "Task"))),
    class = "mfrm_fit"), diagnostics = list(obs = d))
}

test_that("ordinary pair screens use distinct Persons and the declared reference calculation", {
  d <- halo_person_fixture()$diagnostics$obs
  w <- halo_network_wide_scores(d, c("Person", "Task"), "Rater", "Criterion")
  for (method in c("pearson", "spearman", "kendall")) {
    tab <- halo_pair_correlations(w$wide, w$nodes, method = method, persons = w$persons)
    expect_equal(tab$N, rep(6L, 6))
    expect_equal(tab$Persons, tab$N)
    expect_true(all(tab$PValueStatus == "available"))
    ref <- vapply(seq_len(nrow(tab)), function(i) {
      x <- w$wide[[tab$From[i]]]; y <- w$wide[[tab$To[i]]]
      if (method == "pearson") {
        x <- x - mean(x); y <- y - mean(y)
        r <- sum(x * y) / sqrt(sum(x^2) * sum(y^2))
        2 * stats::pt(-abs(r * sqrt(4 / (1 - r^2))), df = 4)
      } else suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)$p.value)
    }, numeric(1))
    expect_equal(tab$PValue, ref, tolerance = 1e-13)
    expect_equal(tab$PAdjusted, pmin(1, 6 * ref), tolerance = 1e-13)
  }
})

test_that("repeated Person contexts preserve correlations but cannot manufacture evidence", {
  skip_if_not_installed("igraph")
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) invisible(NULL))
  d <- halo_person_fixture()
  one <- rater_halo_network_analysis(d$fit, d$diagnostics, method = "pearson")
  d <- halo_person_fixture(tasks = 5)
  repeated <- rater_halo_network_analysis(d$fit, d$diagnostics, method = "pearson")
  key <- function(x) x[order(x$From, x$To), ]
  expect_equal(key(one$pair_metrics)$Estimate, key(repeated$pair_metrics)$Estimate)
  expect_equal(repeated$pair_metrics$N, rep(30L, 6))
  expect_equal(repeated$pair_metrics$Persons, rep(6L, 6))
  expect_true(all(repeated$pair_metrics$PValueStatus == "repeated_person_contexts"))
  expect_true(all(is.na(repeated$pair_metrics$PValue)))
  expect_true(all(is.na(repeated$pair_metrics$PAdjusted)))
  expect_equal(repeated$summary$AvailablePairs, 6)
  expect_equal(repeated$summary$PValueAvailablePairs, 0)
  expect_equal(repeated$summary$PValueUnavailablePairs, 6)
  expect_equal(repeated$summary$Edges, 0)
  expect_false(any(repeated$halo_summary_by_rater$ReviewStatus == "ok"))
  expect_true(all(repeated$halo_summary_by_rater$UnavailableIncidentPValues == 5))
  expect_match(paste(repeated$halo_summary_by_rater$ReviewReason, collapse = " "), "cannot pass")
  grDevices::pdf(NULL, width = 7, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(repeated, type = "network"))
  expect_no_error(plot(repeated, type = "matrix"))
  # Public no-test mode keeps the descriptive edges, including unavailable tests.
  descriptive <- rater_halo_network_analysis(d$fit, d$diagnostics, method = "pearson",
    alpha = 1, positive_only = FALSE)
  expect_true(all(descriptive$pair_metrics$RetainedByP))
  expect_equal(descriptive$summary$Edges, 6)
  expect_false(descriptive$settings$p_filter_applied)
  expect_true(all(is.na(descriptive$edge_metrics$PValue)))
  expect_identical(descriptive$summary$MeanHaloWeight, repeated$summary$MeanHaloWeight)
  expect_match(paste(plot(descriptive, draw = FALSE)$data$notes, collapse = " "), "filter is disabled")
  # Aggregating to one context per Person deliberately changes the target.
  pooled <- rater_halo_network_analysis(d$fit, d$diagnostics, context_facets = "Person")
  expect_equal(pooled$pair_metrics$N, rep(6L, 6))
  expect_true(all(pooled$pair_metrics$PValueStatus == "available"))
  no_person <- rater_halo_network_analysis(d$fit, d$diagnostics, context_facets = "Task",
    min_pair_n = 2, alpha = 1)
  expect_true(all(is.na(no_person$pair_metrics$Persons)))
  # Nonconstant task profiles without Person are still descriptive.
  d$diagnostics$obs$Observed <- seq_len(nrow(d$diagnostics$obs))
  no_person <- rater_halo_network_analysis(d$fit, d$diagnostics, context_facets = "Task", alpha = 1)
  expect_true(all(no_person$pair_metrics$PValueStatus == "person_context_required"))
  expect_equal(no_person$summary$Edges, 6)
})

test_that("pairwise independence uses only shared finite contexts and preserves row identity", {
  d <- halo_person_fixture(tasks = 2)$diagnostics$obs
  d$Observed[d$Task == "t2" & d$Rater == "r1" & d$Criterion == "c1"] <- NA_real_
  w <- halo_network_wide_scores(d[nrow(d):1, ], c("Person", "Task"), "Rater", "Criterion")
  tab <- halo_pair_correlations(w$wide, w$nodes, persons = w$persons)
  one_context <- tab$From == "r1::c1" | tab$To == "r1::c1"
  expect_equal(sum(one_context), 3)
  expect_true(all(tab$PValueStatus[one_context] == "available"))
  expect_equal(tab$N[one_context], rep(6L, 3))
  expect_true(all(tab$PValueStatus[!one_context] == "repeated_person_contexts"))
  expect_equal(tab$N[!one_context], rep(12L, 3))
  expect_equal(tab$PAdjusted[one_context],
    p.adjust(tab$PValue[one_context], method = "bonferroni"))
})

test_that("two-point and old halo results cannot misrepresent inference", {
  skip_if_not_installed("igraph")
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) invisible(NULL))
  d <- halo_person_fixture(persons = 2)
  out <- rater_halo_network_analysis(d$fit, d$diagnostics, method = "pearson",
    min_pair_n = 2, alpha = 1, positive_only = FALSE)
  expect_true(all(is.finite(out$pair_metrics$Estimate)))
  expect_true(all(is.na(out$pair_metrics$PValue)))
  expect_true(all(out$pair_metrics$PValueStatus == "numerically_unavailable"))
  expect_equal(out$summary$Edges, 6)
  p <- plot(out, type = "matrix", draw = FALSE)
  expect_equal(p$data$summary$PValueUnavailablePairs, 6)
  old <- out; old$summary$PValueAvailablePairs <- NULL
  expect_error(summary(old), "Recreate")
  expect_error(plot(old, draw = FALSE), "Recreate")
  empty <- empty_halo_network_bundle(out$settings)
  expect_no_error(summary(empty))
  expect_equal(empty$summary$PValueAvailablePairs, 0)
})
