test_that("signed pair alpha retains chance disagreement and undefined pairs", {
  skip_if_not_installed("igraph")
  obs <- data.frame(Person = rep(1:3, 4), Rater = rep(c("a", "b", "c", "d"), each = 3),
    Observed = c(0, 1, 2, 2, 1, 0, 1, 1, 1, 1, 1, 1))
  fit <- structure(list(config = list(facet_names = "Rater")), class = "mfrm_fit")
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) NULL,
    compute_prob_matrix = function(...) NULL, .package = "mfrmr")
  a <- rater_network_analysis(fit, list(obs = obs), weight_metric = "KrippendorffAlpha")
  ab <- subset(a$pair_metrics, Rater1 == "a" & Rater2 == "b")
  # Independent interval-alpha reference: Do = 8/3, De = 2*var(pooled) = 8/5.
  expect_equal(ab$KrippendorffAlpha, -2/3)
  expect_true(any(a$edge_metrics$SignedWeight < 0))
  expect_equal(a$node_metrics$SignedStrength,
    a$node_metrics$PositiveStrength - a$node_metrics$NegativeStrength)
  expect_true(all(is.na(a$node_metrics$Betweenness)))
  cd <- subset(a$pair_metrics, Rater1 == "c" & Rater2 == "d")
  expect_true(is.na(cd$KrippendorffAlpha))
  expect_false(cd$AlphaStatus == "available")
  p <- plot(a, draw = FALSE)
  expect_true(any(p$data$edges$SignedWeight < 0))
  expect_equal(plot(a, type = "matrix", draw = FALSE)$data$matrix["a", "b"], -2/3)
  expect_identical(mfrm_report(a)$tables$pair_metrics, a$pair_metrics)
  expect_error(rater_network_analysis(fit, list(obs = obs),
    weight_metric = "KrippendorffAlpha", alpha_level = "ordinal"), "categor")
  a2 <- rater_network_analysis(fit, list(obs = transform(obs, Observed = 10 + 3*Observed)),
    weight_metric = "KrippendorffAlpha")
  expect_equal(a2$pair_metrics$KrippendorffAlpha, a$pair_metrics$KrippendorffAlpha)
})

test_that("magnitude indices obey equations, scaling and undefined-ratio boundaries", {
  w <- data.frame(a = c(0, 3), b = c(2, 1), c = c(1, 2))
  p <- rater_network_direction_pairs(w, names(w)); p$EligiblePair <- TRUE
  z <- .mfrm_network_fit_indices(p, names(w))
  # Symmetric differences give S=1: mean absolute differences minus 1,
  # and mean of (positive difference/2)^2, respectively.
  expect_equal(z$SeverityIndex, rep(0, 3))
  expect_equal(z$SWS, c(.5, .5, 0))
  expect_equal(z$PUS, c(.625, .625, .25))
  p2 <- rater_network_direction_pairs(w * 3 + 7, names(w)); p2$EligiblePair <- TRUE
  z2 <- .mfrm_network_fit_indices(p2, names(w))
  expect_equal(z2$SeverityIndex, z$SeverityIndex)
  expect_equal(z2$SWS + 1, 3 * (z$SWS + 1))
  expect_equal(z2$PUS, 9 * z$PUS)
  reverse <- rater_network_direction_pairs(-w, names(w)); reverse$EligiblePair <- TRUE
  expect_equal(.mfrm_network_fit_indices(reverse, names(w))$PUS, z$PUS)
  q <- rater_network_direction_pairs(data.frame(a = 1:4, b = 1:4, c = 2:5), c("a", "b", "c"))
  q$EligiblePair <- TRUE
  unavailable <- .mfrm_network_fit_indices(q, c("a", "b", "c", "isolate"))
  expect_true(all(is.na(unavailable$SWS)))
  expect_equal(unavailable$PUS[1:3], rep(0, 3))
  expect_true(is.na(unavailable$PUS[4]))
  expect_identical(unavailable$SWSStatus[4], "no_comparisons")
  small <- rater_network_direction_pairs(w/10, names(w)); small$EligiblePair <- TRUE
  expect_true(all(.mfrm_network_fit_indices(small, names(w))$SWS < 0))
  # Unequal overlaps and asymmetric strengths test the actual ratio orientation.
  asymmetric <- data.frame(Rater1 = c("a", "a", "b"), Rater2 = c("b", "c", "c"),
    DirectionN = c(4, 2, 3), Rater1HigherSum = c(6, 2, 0),
    Rater2HigherSum = c(2, 1, 3), EligiblePair = TRUE)
  az <- .mfrm_network_fit_indices(asymmetric, c("a", "b", "c"))
  expect_equal(az$SeverityIndex, -log(c(2.5, .2, 1.5)))
  expect_equal(az$SWS, c((.48 + 25 + 1.2 + 5/3)/6 - 1,
    (25 + .48 + .4)/7 - 1, (5/3 + 1.2 + .4)/5 - 1))
  expect_equal(az$PUS, c(.625, .375, .25))
})

test_that("public magnitude graphs do not filter the network fit calculations", {
  skip_if_not_installed("igraph")
  obs <- data.frame(Person = rep(1:2, 3), Rater = rep(c("a", "b", "c"), each = 2),
    Observed = c(0, 3, 2, 1, 1, 2))
  fit <- structure(list(config = list(facet_names = "Rater")), class = "mfrm_fit")
  local_mocked_bindings(mfrm_results_validate_diagnostics_identity = function(...) NULL,
    compute_prob_matrix = function(...) NULL, .package = "mfrmr")
  a <- rater_network_analysis(fit, list(obs = obs), mode = "severity_direction",
    weight_metric = "DirectionalMagnitude", severity_continuity = 0)
  b <- rater_network_analysis(fit, list(obs = obs), mode = "severity_direction",
    weight_metric = "DirectionalMagnitude", min_weight = 2)
  expect_equal(a$network_fit$SWS, c(.5, .5, 0))
  expect_identical(b$network_fit, a$network_fit)
  expect_equal(nrow(b$edge_metrics), 0)
  expect_identical(mfrm_report(a)$tables$network_fit, a$network_fit)
  expect_equal(sort(plot(a, type = "fit", draw = FALSE)$data$table$SWS), c(0, .5, .5))
  flat <- rater_network_analysis(fit, list(obs = transform(obs, Observed = 1)),
    mode = "severity_direction", weight_metric = "DirectionalMagnitude")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(p <- plot(flat, type = "fit"))
  expect_false(any(p$data$available))
  expect_true(all(is.na(p$data$table$SWS)))
  expect_true(all(p$data$table$SWSStatus == "undefined_strength_ratio"))
})

test_that("network routes and printed metrics retain their specific meanings", {
  guide <- mfrmr_output_guide("network")
  for (fn in c("rater_network_analysis", "rater_halo_network_analysis", "fit_rater_network")) {
    expect_true(any(grepl(fn, guide$MainFunction, fixed = TRUE)))
  }
  shown <- mfrm_descriptive_display(data.frame(
    WeightMetric = c("KrippendorffAlpha", "DirectionalMagnitude")))
  expect_identical(shown$WeightMetric,
    c("Signed pairwise Krippendorff alpha", "Mean positive score difference"))
  plots <- mfrmr_output_guide("plots")
  expect_true(all(c("rater_network_analysis", "fit_rater_network") %in% plots$ResultFunction))
})

test_that("halo bootstrap resamples whole persons and preserves failed targets", {
  n <- 32
  w <- data.frame(a = sin(seq_len(n)), b = cos(seq_len(n)/2),
    c = sin(seq_len(n)/3), d = cos(seq_len(n)))
  nodes <- data.frame(Node = names(w), Rater = rep(c("r1", "r2"), each = 2),
    Criterion = rep(c("c1", "c2"), 2))
  persons <- rep(seq_len(n/2), each = 2)
  pairs <- halo_pair_correlations(w, nodes, method = "pearson", persons = persons)
  set.seed(2); old <- .Random.seed
  z <- .mfrm_halo_contrast(list(wide = w, persons = persons), pairs, "pearson", 5, 199, .95, 31, "different_criteria")
  expect_identical(.Random.seed, old)
  expect_identical(z$table$Status, "available")
  expect_equal(z$table$ReferencePairs, 2)
  # Reconstruct the first resample independently from sampled PERSON IDs.
  set.seed(31); selected <- sample.int(16, 16, replace = TRUE)
  # split() orders character group labels lexicographically.
  blocks <- split(seq_len(n), persons)
  sampled <- w[unlist(blocks[selected]), ]
  corrs <- cor(sampled)
  expected <- (corrs["a","b"] + corrs["c","d"] - corrs["a","d"] - corrs["b","c"]) / 2
  expect_equal(z$draws$Estimate[1], expected)
  expect_equal(c(z$table$Lower, z$table$Upper), unname(quantile(z$draws$Estimate, c(.025, .975))))
  all <- .mfrm_halo_contrast(list(wide = w, persons = persons), pairs, "pearson", 5, 0, .95, 31, "all")
  expect_equal(all$table$ReferencePairs, 4)
  expect_identical(all$table$Status, "not_requested")
  sparse <- w[1:4, ]; sp <- halo_pair_correlations(sparse, nodes, min_pair_n = 2)
  fail <- .mfrm_halo_contrast(list(wide = sparse, persons = 1:4), sp, "pearson", 2, 199, .95, 31, "all")
  expect_identical(fail$table$Status, "incomplete_bootstrap")
  expect_gt(fail$table$Failed, 0)
  expect_equal(fail$table$Valid + fail$table$Failed, 199)
  expect_true(is.na(fail$table$Lower))
  empty <- empty_halo_network_bundle(list(bootstrap = 999, conf_level = .95,
    method = "pearson", min_pair_n = 5, halo_reference = "all"))
  expect_identical(empty$halo_contrast$Status, "no_comparable_edges")
  expect_equal(empty$halo_contrast$Requested, 999)
  expect_equal(empty$halo_contrast$Attempted, 0)
  expect_identical(mfrm_report(empty)$tables$halo_contrast, empty$halo_contrast)
})

network_ergm_fixture <- function() {
  ids <- paste0("r", 1:10)
  p <- as.data.frame(t(combn(ids, 2)), stringsAsFactors = FALSE)
  names(p) <- c("Rater1", "Rater2")
  p$Exact <- rep(c(.1, .6, .9), length.out = nrow(p))
  p$N <- 30L; p$EligiblePair <- TRUE
  structure(list(pair_metrics = p, node_metrics = data.frame(Rater = ids),
    settings = list(weight_metric = "Exact", directed = FALSE, min_pair_n = 1)),
    class = "mfrm_rater_network")
}

test_that("binary ERGM agrees with the independent Bernoulli oracle and exports", {
  skip_if_not_installed("ergm"); skip_if_not_installed("network")
  x <- network_ergm_fixture()
  z <- fit_rater_network(x, ~ edges, threshold = .6, seed = 20,
    gof = ~degree, gof_control = ergm::control.gof.ergm(nsim = 9, seed = 20))
  expect_equal(sum(z$dyads$Tie), 15)
  expect_equal(z$coefficients$Estimate, qlogis(1/3), tolerance = 1e-7)
  # Backend IRLS retains a last-iteration Hessian; absolute SE error < 1e-6.
  expect_lt(abs(z$coefficients$SE - sqrt(1/15 + 1/30)), 1e-6)
  expect_true(z$diagnostics$MPLEisMLE)
  expect_false(z$diagnostics$BackendFailure)
  expect_equal(unname(coef(z)), z$coefficients$Estimate)
  expect_equal(vcov(z), stats::vcov(z$fit))
  ci <- confint(z)
  expect_equal(as.vector(ci), c(z$coefficients$Lower, z$coefficients$Upper))
  expect_equal(as.vector(confint(z, parm = "edges", level = .9)),
    qlogis(1/3) + c(-1,1) * qnorm(.95) * z$coefficients$SE, tolerance = 1e-7)
  expect_identical(confint(z, parm = 1), ci)
  expect_identical(attr(ci, "notes"), z$notes)
  unavailable <- z; unavailable$coefficients$Status <- "unavailable"
  expect_true(all(is.na(confint(unavailable))))
  unavailable <- z; unavailable$diagnostics$BackendFailure <- TRUE
  expect_true(all(is.na(confint(unavailable))))
  expect_identical(attr(confint(unavailable), "diagnostics")$Status, "unavailable")
  for (bad in list("missing", NA_character_, 0, -1, 1.5, 2, TRUE))
    expect_error(confint(z, parm = bad), "parm")
  for (bad in list(0, 1, NA_real_, c(.9,.95), ".95", matrix(.95)))
    expect_error(confint(z, level = bad), "level")
  expect_error(confint(z, unused = TRUE), "must be empty")
  expect_true(nrow(z$gof_draws) > 0)
  expect_identical(mfrm_report(z)$tables$coefficients, z$coefficients)
  expect_identical(build_summary_table_bundle(mfrm_report(z))$tables$dyads, z$dyads)
  expect_identical(plot(z, draw = FALSE)$data$table, z$coefficients)
  expect_identical(plot(z, type = "gof", draw = FALSE)$data$draws, z$gof_draws)
  path <- tempfile(); dir.create(path)
  export_summary_appendix(mfrm_report(z), output_dir = path, prefix = "network")
  f <- list.files(path, "coefficients.*csv$", full.names = TRUE)
  expect_length(f, 1)
  expect_equal(read.csv(f)$Estimate, z$coefficients$Estimate)
  expect_error(fit_rater_network(x, ~edges, threshold = 2), "ties")
  x$pair_metrics$N[1] <- 0; x$pair_metrics$Exact[1] <- NA_real_
  expect_error(fit_rater_network(x, ~edges, threshold = .6), "Unobserved dyads")
  missing_fit <- fit_rater_network(x, ~edges, threshold = .6, missing = "model", seed = 23)
  expect_equal(missing_fit$settings$MissingDyads, 1)
  expect_equal(sum(is.na(missing_fit$dyads$Tie)), 1)
  expect_equal(missing_fit$coefficients$Estimate, qlogis(15/44), tolerance = 1e-6)
  expect_error(fit_rater_network(x, ~edges, threshold = .6, missing = "model",
    node_data = data.frame(Node = "r1", X = 1)), "every Node")
})

test_that("ERGM node covariates match IDs and the independent logistic likelihood", {
  skip_if_not_installed("ergm"); skip_if_not_installed("network")
  x <- network_ergm_fixture()
  nodes <- data.frame(Node = rev(x$node_metrics$Rater), Experience = c(2, 5, 1, 4, 2, 5, 3, 2, 1, 4))
  p <- x$pair_metrics
  covariate <- nodes$Experience[match(p$Rater1, nodes$Node)] + nodes$Experience[match(p$Rater2, nodes$Node)]
  reference <- glm(I(p$Exact > .6) ~ covariate, family = binomial())
  z <- fit_rater_network(x, ~edges + nodecov("Experience"), threshold = .6, node_data = nodes)
  expect_equal(z$coefficients$Estimate, unname(coef(reference)), tolerance = 1e-7)
  expect_identical(as.character(z$nodes$Node), x$node_metrics$Rater)
  expect_equal(network::get.vertex.attribute(z$network, "Experience"),
    nodes$Experience[match(x$node_metrics$Rater, nodes$Node)])
})
