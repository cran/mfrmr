fair_delta_case <- local({
  cached <- list()
  function(model = "RSM") {
    if (is.null(cached[[model]])) {
      if (model == "RSM") {
        d <- load_mfrmr_data("example_core"); facets <- c("Rater", "Criterion")
      } else {
        d <- read.csv(test_path("fixtures", "tam-jml-0.2.4-response.csv"))
        d <- d[d$Case == "PCM_partial_anchor", ]; facets <- c("Rater", "Task")
      }
      f <- fit_mfrm(d, "Person", facets, "Score", method = "JML", model = model,
        step_facet = if (model == "PCM") "Task" else NULL,
        optimizer = "BFGS", maxit = 500, reltol = 1e-12)
      plan <- unique(d[facets]); plan$Weight <- seq_len(nrow(plan))
      a <- data.frame(Facet = "Rater", Level = f$config$facet_levels$Rater[1],
        Anchor = f$facets$others$Estimate[which(f$facets$others$Facet == "Rater")[1]])
      cached[[model]] <<- list(fit = f, data = d, plan = plan, anchors = a)
    }
    cached[[model]]
  }
})

# Independent category polynomial evaluation; no delta design/target helper.
fair_delta_numeric_scores <- function(par, fit, plan, persons, cfg = fit$config) {
  pa <- mfrmr:::expand_params(par, mfrmr:::build_param_sizes(cfg), cfg)
  vapply(persons, function(person) {
    values <- vapply(seq_len(nrow(plan)), function(i) {
      eta <- pa$theta[match(person, cfg$person_levels)]
      for (facet in cfg$facet_names) eta <- eta + cfg$facet_signs[[facet]] *
        pa$facets[[facet]][match(plan[[facet]][i], cfg$facet_levels[[facet]])]
      steps <- if (cfg$model == "RSM") pa$steps else
        pa$steps_mat[match(plan[[cfg$step_facet]][i], cfg$facet_levels[[cfg$step_facet]]), ]
      logits <- (seq_len(cfg$n_cat) - 1L) * eta - c(0, cumsum(steps))
      p <- exp(logits - max(logits)); p <- p / sum(p)
      sum(p * (seq_len(cfg$n_cat) - 1L)) + cfg$rating_min
    }, 0)
    weighted.mean(values, plan$Weight)
  }, 0)
}

fair_delta_fd <- function(fn, x, h = 1e-5) {
  if (!length(x)) return(matrix(0, length(fn(x)), 0L))
  vapply(seq_along(x), function(j) {
    up <- down <- x; up[j] <- up[j] + h; down[j] <- down[j] - h
    (fn(up) - fn(down)) / (2 * h)
  }, numeric(length(fn(x))))
}

fair_delta_constraint_config <- function(fit, anchors) {
  cfg <- fit$config
  spec <- mfrmr:::prepare_constraint_specs(fit$prep, anchor_df = anchors,
    noncenter_facet = cfg$noncenter_facet, dummy_facets = cfg$dummy_facets,
    centering_scope = cfg$centering_scope)
  cfg$theta_spec <- spec$theta_spec; cfg$facet_specs <- spec$facet_specs
  cfg
}

test_that("delta covariance recovers the independent observed-mean identity", {
  case <- fair_delta_case(); f <- case$fit
  equal_plan <- case$plan[f$config$facet_names]
  fair <- fair_average_table(f, reference_plan = equal_plan)
  persons <- f$config$person_levels[1:2]
  set.seed(65); rng <- .Random.seed
  ci <- confint(fair, parm = persons, method = "delta")
  expect_identical(.Random.seed, rng)
  expect_true(all(ci$Available))
  expect_false(any(ci$CIEligible))
  p <- mfrmr:::mfrm_fair_generator_probs(f); k <- 0:(ncol(p) - 1)
  v <- drop(p %*% k^2) - drop(p %*% k)^2
  expected <- vapply(persons, function(person) {
    rows <- as.character(f$prep$data$Person) == person
    sum(v[rows]) / sum(rows)^2
  }, 0)
  expect_equal(diag(attr(ci, "covariance")), expected, tolerance = 1e-9)
  expect_lt(abs(attr(ci, "covariance")[1, 2]), 1e-10)
  expect_equal(ci$SE^2, unname(expected), tolerance = 1e-9)
  expect_equal(ci$CI_Lower, ci$Estimate - qnorm(.975) * ci$SE, tolerance = 1e-12)
  expect_true(all(is.na(ci$Planned)))
  expect_error(confint(fair, method = "delta", nsim = 2), "omit")
  expect_error(confint(fair, method = "delta", seed = 1), "omit")
})

test_that("RSM and PCM information and target derivatives match finite differences", {
  for (model in c("RSM", "PCM")) {
    case <- fair_delta_case(model); f <- case$fit; cfg <- f$config
    plan <- mfrmr:::mfrm_fair_reference_plan(f, case$plan)
    persons <- cfg$person_levels[3:4]
    design <- mfrmr:::mfrm_fair_delta_design(f)
    p <- mfrmr:::mfrm_fair_generator_probs(f)
    D <- mfrmr:::mfrm_fair_delta_center(design$category, p)
    H <- mfrmr:::mfrm_fair_delta_information(D, p)
    sizes <- mfrmr:::build_param_sizes(cfg)
    fdH <- stats::optimHess(f$opt$par,
      fn = function(par) mfrmr:::mfrm_loglik_jml(par, design$idx, cfg, sizes),
      gr = function(par) mfrmr:::mfrm_grad_jml(par, design$idx, cfg, sizes),
      control = list(ndeps = rep(1e-5, length(f$opt$par))))
    expect_equal(unname(as.matrix(H)), unname(fdH), tolerance = 1e-7)
    target <- mfrmr:::mfrm_fair_delta_target(f, plan, persons,
      case$anchors[FALSE, ], Matrix::Matrix(0, 0, 0, sparse = TRUE))
    fdG <- fair_delta_fd(function(par) fair_delta_numeric_scores(par, f, plan, persons), f$opt$par)
    expect_equal(unname(target$free), unname(fdG), tolerance = 1e-8)
    ci <- confint(fair_average_table(f, reference_plan = case$plan), parm = persons, method = "delta")
    direct <- fdG %*% solve(fdH, t(fdG))
    expect_equal(unname(attr(ci, "covariance")), unname(direct), tolerance = 1e-8)
    expect_equal(ci$Estimate, unname(fair_delta_numeric_scores(f$opt$par, f, plan, persons)), tolerance = 1e-12)
  }
})

test_that("anchor derivatives follow each actual centering map", {
  case <- fair_delta_case(); source <- case$fit
  for (scope in c("facet", "remainder", "anchor_origin")) {
    args <- source$config$replay_inputs; args$package_version <- NULL
    args$data <- case$data; args$anchors <- case$anchors; args$centering_scope <- scope
    fit <- do.call(fit_mfrm, args)
    design <- mfrmr:::mfrm_fair_delta_design(fit)
    A <- mfrmr:::mfrm_fair_delta_anchor_map(source, case$anchors)
    O <- mfrmr:::mfrm_fair_delta_anchor_design(fit, design$idx, case$anchors, A)
    p <- mfrmr:::mfrm_fair_generator_probs(fit)
    B <- mfrmr:::mfrm_fair_delta_information(
      mfrmr:::mfrm_fair_delta_center(design$category, p), p,
      mfrmr:::mfrm_fair_delta_center(O, p))
    cfg_at <- function(values) {
      a <- case$anchors; a$Anchor <- values; fair_delta_constraint_config(fit, a)
    }
    fdB <- fair_delta_fd(function(values) mfrmr:::mfrm_grad_jml(fit$opt$par,
      design$idx, cfg_at(values), mfrmr:::build_param_sizes(fit$config)), case$anchors$Anchor)
    expect_equal(unname(as.matrix(B)), unname(fdB), tolerance = 1e-7)
    plan <- mfrmr:::mfrm_fair_reference_plan(fit, case$plan)
    persons <- source$config$person_levels[1:2]
    target <- mfrmr:::mfrm_fair_delta_target(fit, plan, persons, case$anchors, A)
    fdGa <- fair_delta_fd(function(values) fair_delta_numeric_scores(fit$opt$par, fit,
      plan, persons, cfg_at(values)), case$anchors$Anchor)
    expect_equal(unname(target$anchors), unname(fdGa), tolerance = 1e-8)
  }
})

test_that("same-data estimated anchors recover joint covariance including Person anchors", {
  case <- fair_delta_case(); f <- case$fit; persons <- f$config$person_levels[1:2]
  joint <- confint(fair_average_table(f, reference_plan = case$plan), parm = persons, method = "delta")
  all_raters <- subset(f$facets$others, Facet == "Rater", c("Facet", "Level", "Estimate"))
  names(all_raters)[3] <- "Anchor"
  person_anchor <- data.frame(Facet = "Person", Level = persons[1],
    Anchor = f$facets$person$PrimaryEstimate[1])
  for (anchors in list(case$anchors, all_raters, person_anchor)) {
    args <- f$config$replay_inputs; args$package_version <- NULL
    args$data <- case$data; args$anchors <- anchors
    scored <- do.call(fit_mfrm, args)
    fair <- fair_average_table(scored, reference_plan = case$plan)
    n <- nrow(f$prep$data)
    out <- confint(fair, parm = persons, method = "delta", anchor_source = f,
      response_map = list(calibration = seq_len(n), scoring = seq_len(n)))
    expect_true(all(out$Available))
    expect_equal(out$Estimate, joint$Estimate, tolerance = 2e-6)
    # Separate optimizer stops need not give identical coordinates: compare in
    # squared-score units (this bound changes either SE by less than 3e-7).
    expect_lt(max(abs(attr(out, "covariance") - attr(joint, "covariance"))), 1e-7)
    expect_equal(Reduce(`+`, attr(out, "covariance_components")), attr(out, "covariance"), tolerance = 1e-12)
    if (identical(anchors, case$anchors)) {
      # Two distinct response occasions can have identical labels/observed
      # scores. Only shared event indices should induce cross-stage covariance.
      p <- mfrmr:::mfrm_fair_generator_probs(f)
      L <- diag(2); dimnames(L) <- list(persons, persons)
      separate <- mfrmr:::mfrm_fair_delta_covariance(scored, f, fair$reference_plan,
        persons, L, anchors, list(calibration = seq_len(n), scoring = n + seq_len(n)), rbind(p, p))
      expect_equal(unname(separate$components$shared_events), matrix(0, 2, 2))
      expect_equal(separate$covariance, separate$components$direct_scoring +
        separate$components$through_calibration, tolerance = 1e-12)
    }
  }
})

test_that("partial overlap preserves cross-stage covariance and paired contrasts", {
  case <- fair_delta_case(); union <- case$fit
  person <- as.character(union$prep$data$Person)
  calibration_rows <- which(person %in% union$config$person_levels[1:32])
  scoring_rows <- which(person %in% union$config$person_levels[17:48])
  args <- union$config$replay_inputs; args$package_version <- NULL
  args$data <- union$prep$data[calibration_rows, ]
  calibration <- do.call(fit_mfrm, args)
  args$anchors <- case$anchors
  args$anchors$Anchor <- mfrmr:::mfrm_fair_anchor_values(calibration, args$anchors)
  args$data <- union$prep$data[scoring_rows, ]
  scoring <- do.call(fit_mfrm, args)
  fair <- fair_average_table(scoring, reference_plan = case$plan)
  persons <- c("P025", "P040"); L <- matrix(c(1, -1), 1, dimnames = list("difference", persons))
  maps <- list(calibration = calibration_rows, scoring = scoring_rows)
  ci <- confint(fair, parm = persons, method = "delta", anchor_source = calibration,
    generating_fit = union, response_map = maps)
  contrast <- confint(fair, parm = persons, method = "delta", contrasts = L,
    anchor_source = calibration, generating_fit = union, response_map = maps)
  expect_true(all(ci$Available) && all(contrast$Available))
  expect_equal(attr(ci, "settings")$shared_responses, length(intersect(calibration_rows, scoring_rows)))
  expect_gt(max(abs(attr(ci, "covariance_components")$shared_events)), 1e-6)
  expect_equal(attr(contrast, "covariance"), L %*% attr(ci, "covariance") %*% t(L), tolerance = 1e-12)
  expect_equal(contrast$Estimate, drop(L %*% ci$Estimate), ignore_attr = TRUE)
  expect_equal(attr(ci, "response_map"), maps)
  expect_error(confint(fair, method = "delta"), "without the calibration fit")
})

test_that("unavailable delta information preserves points and saved diagnostics", {
  case <- fair_delta_case(); fair <- fair_average_table(case$fit, reference_plan = case$plan)
  testthat::local_mocked_bindings(mfrm_fair_delta_solve = function(...) stop("Controlled information failure"),
    .package = "mfrmr")
  ci <- confint(fair, parm = "P001", method = "delta")
  expect_false(ci$Available)
  expect_identical(ci$Status, "unavailable_information")
  expect_match(ci$Reason, "Controlled information failure")
  expect_true(is.finite(ci$Estimate) && is.na(ci$SE) && is.na(ci$CI_Lower))
  expect_true(all(is.na(attr(ci, "covariance"))))
  expect_identical(summary(ci), ci)
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(ci, path); saved <- readRDS(path)
  p <- plot(saved, draw = FALSE)
  expect_equal(p$data$data$Reason, ci$Reason)
  expect_match(p$data$notes$Text, "unavailable information")
})

test_that("a nonselected extreme generator Person preserves finite target points", {
  case <- fair_delta_case(); d <- case$data
  d$Score[d$Person == tail(case$fit$config$person_levels, 1)] <- max(d$Score)
  fit <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", method = "JML",
    optimizer = "BFGS", maxit = 1000, reltol = 1e-13)
  fair <- fair_average_table(fit, reference_plan = case$plan)
  ci <- confint(fair, parm = "P001", method = "delta")
  expect_true(any(is.infinite(fit$facets$person$PrimaryEstimate)))
  expect_equal(ci$Estimate, fair$reference_scores$Estimate[match("P001", fair$reference_scores$Person)])
  expect_false(ci$Available)
  expect_identical(ci$Status, "unavailable_information")
  expect_match(ci$Reason, "every Person in the generating fit")
  expect_true(is.na(ci$SE) && all(is.na(attr(ci, "covariance"))))
  expect_error(confint(fair, parm = "P001", nsim = 2, seed = 75), "generating fit needs finite")
})
