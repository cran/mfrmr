gmfrm_test_data <- function(max_score = 2L) {
  d <- expand.grid(Person = paste0("p", 1:4), Task = paste0("t", 1:3), Rater = paste0("r", 1:3))
  d$Score <- rep(0:max_score, length.out = nrow(d))
  d
}

gmfrm_test_problem <- function(max_score = 2L) {
  d <- gmfrm_test_data(max_score)
  list(data = d, problem = mfrm_gmfrm_problem(d, max_score, gauss_hermite_normal(9L)))
}

gmfrm_numeric_gradient <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    lo <- hi <- x
    lo[j] <- lo[j] - h
    hi[j] <- hi[j] + h
    (f(hi) - f(lo)) / (2 * h)
  }, numeric(1))
}

test_that("two-family kernel reuse preserves changing E-step counts and invalid trials", {
  evaluations <- 0L
  original <- mfrm_jml_probability_bundle
  local_mocked_bindings(mfrm_jml_probability_bundle = function(...) {
    evaluations <<- evaluations + 1L
    original(...)
  }, .package = "mfrmr")
  p <- gmfrm_test_problem()$problem
  x <- p$start + seq(-.1, .1, length.out = length(p$start))
  initial <- evaluations
  state <- p$marginal(x)
  expect_identical(evaluations - initial, 1L)
  first <- p$complete(x, state$counts)
  doubled <- p$complete(x, 2 * state$counts)
  expect_identical(evaluations - initial, 1L)
  expect_equal(doubled$value, 2 * first$value, tolerance = 0)
  expect_equal(doubled$gradient, 2 * first$gradient, tolerance = 0)
  changed <- x; changed[1] <- changed[1] + .2
  expect_false(identical(p$kernel(changed), p$kernel(x)))
  bad <- x; bad[grepl("^log_slope:Rater", names(bad))] <- 1000
  expect_null(p$kernel(bad))
  expect_identical(p$marginal(bad)$value, Inf)
  expect_identical(p$marginal(x), state)
})

test_that("common product parameters retain owners and reject non-product crossings", {
  p <- gmfrm_test_problem()$problem
  setup <- p$common
  x <- unname(p$start + seq(-.2, .3, length.out = length(p$start)))
  params <- expand_params(x, setup$sizes, setup$config)
  expect_equal(collapse_expanded_params(params, setup$config), x, tolerance = 1e-12)
  expect_named(params$slope_components, c("Task", "Rater"))
  expect_equal(params$slope_components$Task$Level, p$levels$Task)
  expect_equal(params$slope_components$Rater$Level, p$levels$Rater)
  expect_equal(prod(params$slope_components$Task$Slope), 1)
  expect_gt(abs(prod(params$slope_components$Rater$Slope) - 1), .1)
  cells <- setup$config$gpcm_spec$cells
  expect_equal(params$slopes, params$slope_components$Task$Slope[cells$Task] *
    params$slope_components$Rater$Slope[cells$Rater])
  expect_s4_class(setup$config$gpcm_spec$log_slope_design, "sparseMatrix")
  params$log_slopes[1] <- params$log_slopes[1] + .1
  expect_error(collapse_expanded_params(params, setup$config), "product structure")
  expect_error(expand_gpcm_log_slopes(1, setup$config$gpcm_spec), "free coordinate")
})

test_that("common MML likelihood and scores agree with the EM marginal calculation", {
  for (k in c(1L, 2L, 4L)) for (sparse in if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) c(FALSE, TRUE) else k != 2L) {
    d <- gmfrm_test_data(k)
    if (sparse) d <- d[!(d$Task == "t3" & d$Rater == "r2"), ]
    # Existing factor ordering and input order must not relabel the coordinates.
    d$Task <- factor(d$Task, levels = c("t3", "t1", "t2"))
    d <- d[nrow(d):1, ]
    # This is an algebraic evaluator/gradient comparison, not a grid-accuracy
    # claim. All category counts and sparse designs remain in the CRAN run.
    order <- if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) 31L else 9L
    quad <- gauss_hermite_normal(order)
    p <- mfrm_gmfrm_problem(d, k, quad)
    setup <- p$common
    x <- unname(p$start + seq(-.13, .17, length.out = length(p$start)))
    direct <- make_mfrm_direct_evaluator("MML",
      make_param_cache(setup$sizes, setup$config, setup$idx, is_mml = TRUE),
      setup$idx, setup$config, setup$sizes, quad)
    expect_equal(direct$value(x), p$marginal(x)$value * p$n_person, tolerance = 1e-11)
    expect_equal(direct$gradient(x), p$marginal(x)$gradient * p$n_person, tolerance = 1e-10)
    expect_equal(direct$gradient(x), gmfrm_numeric_gradient(direct$value, x), tolerance = 1e-7)
    expect_equal(collapse_expanded_params(expand_params(x, setup$sizes, setup$config),
      setup$config), x, tolerance = 1e-12)
    expect_equal(nrow(setup$config$gpcm_spec$cells), if (sparse) 8L else 9L)
  }
})

test_that("two-family E-step counts retain prior weights at extreme finite likelihoods", {
  d <- expand.grid(Person = c("b", "a", "c"), Task = c("t1", "t2"), Rater = c("r1", "r2"))
  d$Score <- as.integer(d$Person == "a")
  rule <- list(nodes = c(-1, 0, 1), weights = c(.2, .5, .3))
  for (offset in c(1e12, 1e16)) {
    p <- mfrm_gmfrm_problem(d, 1L, rule)
    x <- p$start; x[grepl("location:Rater:", names(x))] <- -offset
    actual <- p$marginal(x)
    # All four ratings of b/c are category zero; a has category one. At
    # 1e16 the node displacement rounds away, so the E step must equal the
    # prior. At 1e12 the likelihood ratios still contribute exp(-4 * theta).
    lower <- if (offset == 1e16) rule$weights else rule$weights * exp(-4 * rule$nodes)
    lower <- lower / sum(lower)
    expected <- rbind(rule$weights, lower, lower)
    expect_equal(unname(actual$posterior), unname(expected), tolerance = 1e-14)
    expect_identical(rownames(actual$posterior), c("a", "b", "c"))
    expect_equal(unname(rowSums(actual$posterior)), rep(1, 3), tolerance = 1e-14)
    expect_equal(sum(actual$counts), nrow(d), tolerance = 1e-14)
    expect_equal(actual$counts,
      cbind(rep(2 * lower, each = 4), rep(rule$weights, each = 4)), tolerance = 1e-14)
    expect_equal(actual$gradient[2:3], rep(-4/3, 2), tolerance = 1e-14)
    expect_equal(p$marginal(x, derivatives = FALSE)$posterior, actual$posterior)
    reordered <- mfrm_gmfrm_problem(d[rev(seq_len(nrow(d))), ], 1L, rule)
    expect_equal(reordered$marginal(x), actual, tolerance = 1e-14)
  }
})

test_that("two-family M-step log-slope scores retain small contrasts at large logits", {
  d <- expand.grid(Person = c("b", "a", "c"), Task = c("t1", "t2"), Rater = c("r1", "r2"))
  d$Score <- rep(0:2, length.out = nrow(d))
  rule <- list(nodes = c(-1, 0, 1), weights = c(.2, .5, .3))
  p <- mfrm_gmfrm_problem(d, 2L, rule)
  x <- p$start
  x[grepl("location:Rater:", names(x))] <- -1e12
  x[grepl("step:Rater:", names(x))] <- -1e12 + 2
  # Category zero has negligible probability; the adjacent log odds between
  # categories 1 and 2 are theta + 2, despite their logits being about 2e12.
  for (category in 1:2) {
    counts <- matrix(0, 12, 3); counts[, category + 1L] <- 1
    odds <- rule$nodes + 2
    row_score <- odds * if (category == 1) plogis(odds) else -plogis(-odds)
    expected <- c(0, rep(2 * sum(row_score) / p$n_person, 2))
    expect_equal(tail(p$complete(x, counts)$gradient, 3), expected, tolerance = 1e-14)
  }
})

test_that("normal population scale changes preserve the identified two-slope likelihood", {
  p <- gmfrm_test_problem()$problem
  x <- p$start + seq(-.2, .3, length.out = length(p$start))
  u <- p$unpack(x)
  d <- gmfrm_test_data()
  for (sd in c(.5, 1, 2)) for (mu in c(-.6, .4)) {
    # Original theta = mu + sd*z. Standardize to z while preserving task GM1.
    # All free rater slopes absorb sd; task/rater locations and steps change units.
    scaled <- x
    locations <- grepl("^(location|step):", names(x))
    raters <- grepl("^location:Rater:", names(x))
    slopes <- grepl("^log_slope:Rater:", names(x))
    scaled[locations] <- x[locations] / sd
    scaled[raters] <- (x[raters] - mu) / sd
    scaled[slopes] <- x[slopes] + log(sd)
    literal <- function(task, rater, z) {
      logits <- c(0, cumsum(u$slopes$Task[task] * u$slopes$Rater[rater] *
        (mu + sd*z - u$locations$Task[task] - u$locations$Rater[rater] - u$steps[rater, ])))
      logp <- logits - max(logits)
      logp - log(sum(exp(logp)))
    }
    grid <- expand.grid(Context = seq_len(nrow(p$cells)), Z = p$nodes)
    expected <- t(vapply(seq_len(nrow(grid)), function(j) {
      cell <- p$cells[grid$Context[j], ]
      literal(p$levels$Task[cell$Task], p$levels$Rater[cell$Rater], grid$Z[j])
    }, numeric(3)))
    expect_equal(p$kernel(scaled)$log_prob, unname(expected), tolerance = 1e-12)
    # Integrate joint Person patterns, not a product of one-rating marginals.
    loglik <- sum(vapply(split(d, d$Person), function(rows) {
      joint <- vapply(p$nodes, function(z) sum(vapply(seq_len(nrow(rows)), function(j)
        literal(as.character(rows$Task[j]), as.character(rows$Rater[j]), z)[rows$Score[j]+1L],
        numeric(1))), numeric(1))
      max(joint) + log(sum(p$weights * exp(joint - max(joint))))
    }, numeric(1)))
    expect_equal(-p$n_person * p$marginal(scaled)$value, loglik, tolerance = 1e-12)
    expect_equal(prod(p$unpack(scaled)$slopes$Task), 1, tolerance = 1e-12)
    expect_equal(unname(p$unpack(scaled)$slopes$Rater), unname(sd*u$slopes$Rater), tolerance = 1e-12)
  }
})

test_that("declared column names preserve the two-slope model and its original identities", {
  fixture <- list(data = gmfrm_test_data())
  fixture$data <- fixture$data[fixture$data$Task != "t3", ]
  reference <- mfrm_gmfrm_problem(fixture$data, 2L, gauss_hermite_normal(9L))
  x <- reference$start + seq(-.2, .3, length.out = length(reference$start))
  grid <- data.frame(Theta = c(-.8, .4), Task = c("t1", "t2"), Rater = c("r2", "r1"))
  for (columns in list(c("Candidate", "Criterion", "Assessor", "Rating"),
                       c("受験者 ID", "観点 名", "Judge-ID", "得点"))) {
    data <- fixture$data
    names(data) <- columns
    data <- data[rev(names(data))] # Model roles cannot depend on column order.
    owners <- columns[2:3]
    p <- mfrm_gmfrm_problem(data, 2L, gauss_hermite_normal(9L),
      slope_facets = owners, person = columns[1], score = columns[4])
    actual <- p$marginal(unname(x)); expected <- reference$marginal(x)
    expect_equal(actual$logLik, expected$logLik, tolerance = 1e-12)
    expect_equal(actual$gradient, expected$gradient, tolerance = 1e-12)
    expect_equal(actual$posterior, expected$posterior, tolerance = 1e-12)
    expect_identical(p$common$config$facet_names, owners)
    expect_identical(p$common$config$step_facet, owners[2])
    expect_identical(names(p$unpack(x)$locations), owners)
    expect_identical(names(p$unpack(x)$slopes), owners)
    expect_identical(p$unpack(x)$step_facet, owners[2])
    query <- grid; names(query)[2:3] <- owners
    response <- p$response(x, query[rev(names(query))])
    base <- reference$response(x, grid)
    expect_equal(response$probabilities, base$probabilities, tolerance = 1e-12)
    expect_equal(response$summary$Information, base$summary$Information, tolerance = 1e-12)
    expect_identical(response$summary$SlopeOwner1, rep(owners[1], 2L))
    expect_identical(response$summary$SlopeOwner2, rep(owners[2], 2L))
    expect_identical(response$summary[[owners[1]]], grid$Task)
    expect_identical(response$summary[[owners[2]]], grid$Rater)
    target <- mfrm_gpcm_product_slope_targets(x, p$common$config)
    expect_identical(names(target$effective$metadata), owners)
    expect_identical(target$components$metadata$Owner, rep(owners, c(2L, 3L)))
    direct <- make_mfrm_direct_evaluator("MML",
      make_param_cache(p$common$sizes, p$common$config, p$common$idx, is_mml = TRUE),
      p$common$idx, p$common$config, p$common$sizes, gauss_hermite_normal(9L))
    expect_equal(direct$value(x), -expected$logLik, tolerance = 1e-12)
    expect_equal(direct$gradient(x), p$n_person*expected$gradient, tolerance = 1e-12)
    fit <- mfrm_gmfrm_em(p, maxit = 2L)
    base_fit <- mfrm_gmfrm_em(reference, maxit = 2L)
    expect_equal(unname(fit$par), unname(base_fit$par), tolerance = 1e-12)
    expect_equal(fit$trace, base_fit$trace, tolerance = 1e-12)
    expect_identical(fit$slopes$SlopeOwner, rep(owners, c(2L, 3L)))
    expect_identical(fit$specification$slope_facets, owners)
    expect_identical(fit$specification$person, columns[1])
    expect_identical(fit$specification$score, columns[4])
  }
  for (owners in list("Task", c("Task", "Task"), c("Task", NA_character_),
                      c("Task", ""), c("Theta", "Rater"), c("Person", "Rater"))) {
    expect_error(mfrm_gmfrm_problem(fixture$data, 2L, gauss_hermite_normal(9L),
      slope_facets = owners), "column|names")
  }
  expect_error(mfrm_gmfrm_problem(fixture$data, 2L, gauss_hermite_normal(9L),
    person = "Score"), "distinct")
})

test_that("facet roles follow the declared order rather than domain labels", {
  d <- gmfrm_test_data()
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(9L), slope_facets = c("Rater", "Task"))
  x <- p$start + seq(-.2, .3, length.out = length(p$start))
  u <- p$unpack(x)
  expect_identical(p$common$config$step_facet, "Task")
  expect_equal(prod(u$slopes$Rater), 1)
  expect_equal(sum(u$locations$Rater), 0)
  expect_gt(abs(prod(u$slopes$Task) - 1), .1)
  expect_identical(rownames(u$steps), p$levels$Task)
  context <- data.frame(Theta = .4, Rater = "r2", Task = "t1")
  logits <- c(0, cumsum(u$slopes$Rater["r2"] * u$slopes$Task["t1"] *
    (.4-u$locations$Rater["r2"]-u$locations$Task["t1"]-u$steps["t1", ])))
  expected <- exp(logits-max(logits)); expected <- expected/sum(expected)
  actual <- p$response(x, context)
  expect_equal(unname(actual$probabilities[1, ]), unname(expected), tolerance = 1e-12)
  expect_identical(actual$summary$SlopeOwner1, "Rater")
  expect_identical(actual$summary$SlopeOwner2, "Task")
})

test_that("adaptive integration differentiates both slope families", {
  d <- gmfrm_test_data()
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(61L))
  setup <- p$common
  setup$config$estimation_control <- list(mml_integration = "adaptive", quad_points = 31L)
  ev <- make_mfrm_direct_evaluator("MML",
    make_param_cache(setup$sizes, setup$config, setup$idx, is_mml = TRUE),
    setup$idx, setup$config, setup$sizes, gauss_hermite_normal(31L))
  x <- unname(p$start + seq(-.13, .17, length.out = length(p$start)))
  u <- p$unpack(x)
  reference <- -sum(vapply(split(d, d$Person), function(rows) {
    integrand <- function(theta) vapply(theta, function(t) {
      lp <- vapply(seq_len(nrow(rows)), function(j) {
        task <- as.character(rows$Task[j]); rater <- as.character(rows$Rater[j])
        adjacent <- u$slopes$Task[task] * u$slopes$Rater[rater] *
          (t - u$locations$Task[task] - u$locations$Rater[rater] - u$steps[rater, ])
        logits <- c(0, cumsum(adjacent))
        logits[rows$Score[j] + 1L] - max(logits) - log(sum(exp(logits - max(logits))))
      }, numeric(1))
      exp(sum(lp)) * dnorm(t)
    }, numeric(1))
    log(integrate(integrand, -Inf, Inf, rel.tol = 1e-11)$value)
  }, numeric(1)))
  expect_equal(ev$value(x), reference, tolerance = 1e-9)
  expect_equal(ev$gradient(x), gmfrm_numeric_gradient(ev$value, x), tolerance = 1e-7)
})

test_that("crossing indices reject relabelled levels and implicit two-family settings", {
  d <- gmfrm_test_data()
  prep <- prepare_mfrm_data(d, "Person", c("Task", "Rater"), "Score")
  owners <- c("Task", "Rater")
  expect_error(build_indices(prep, "Rater", owners), "crossing specification")
  sparse <- prep
  sparse$data <- sparse$data[!(sparse$data$Task == "t3" & sparse$data$Rater == "r2"), ]
  spec <- build_gpcm_product_slope_spec(sparse, owners, "Rater")
  expect_error(build_indices(prep, "Rater", owners, gpcm_spec = spec), "crossings")
  reversed <- sparse
  reversed$levels$Task <- rev(reversed$levels$Task)
  expect_error(build_indices(reversed, "Rater", owners, gpcm_spec = spec), "level ordering")
  expect_error(fit_mfrm(d, person = "Person", facets = owners, score = "Score",
    model = "GPCM", step_facet = "Rater", slope_facet = owners),
    "Two slope families require")
})

test_that("optimizer coordinates distinguish slope owners and their constraints", {
  p <- gmfrm_test_problem()$problem
  s <- p$common
  map <- s$parameter_map
  expect_equal(map$OptimizerIndex, seq_along(p$start))
  slopes <- map[map$Block == "log_slopes", ]
  expect_equal(slopes$Facet, c(rep("Task", 2), rep("Rater", 3)))
  expect_equal(slopes$Level, c("t1", "t2", "r1", "r2", "r3"))
  expect_equal(slopes$ReferenceLevel, c("t3", "t3", "", "", ""))
  expect_equal(slopes$Constraint, c(rep("sum_zero_log_slopes", 2), rep("free_log_slope", 3)))
  expect_equal(anyDuplicated(map$Coordinate), 0L)
  wrong <- s$config
  wrong$slope_facet <- rev(wrong$slope_facet)
  expect_error(mfrmr_gpcm_log_slope_map(wrong, s$sizes), "fitted owners")
  x <- p$start + seq(-.2, .3, length.out = length(p$start))
  audit <- mfrmr_nonlinear_transformation_audit(x, s$sizes, s$config, "log_slopes")
  expect_identical(audit$status, "evaluated_diagnostic_only")
  expect_true(audit$block_summary$LogFullColumnRank)
  expect_true(audit$block_summary$NaturalFullColumnRank)
  expect_lt(audit$block_summary$MaxScaledNaturalJacobianDifference, 1e-8)
  expect_lt(audit$block_summary$InvariantResidual, 1e-14)
  expect_false(audit$structural_identification_classified)
  x[build_param_slices(s$sizes)$log_slopes] <- 1000
  extreme <- mfrmr_nonlinear_transformation_audit(x, s$sizes, s$config, "log_slopes")
  expect_identical(extreme$status, "partially_unavailable")
  expect_false(extreme$structural_identification_classified)
})

test_that("component and product slope targets retain cross-owner covariance", {
  p <- gmfrm_test_problem()$problem
  s <- p$common
  x <- unname(p$start + seq(-.2, .3, length.out = length(p$start)))
  targets <- mfrm_gpcm_product_slope_targets(x, s$config)
  u <- p$unpack(x)
  expect_equal(targets$components$estimate, unname(c(u$slopes$Task, u$slopes$Rater)))
  contexts <- targets$effective$metadata
  expect_equal(targets$effective$estimate,
    unname(u$slopes$Task[contexts$Task] * u$slopes$Rater[contexts$Rater]))
  for (target in c("components", "effective")) {
    numeric <- mfrmr_numeric_transformation_jacobian(
      function(z) mfrm_gpcm_product_slope_targets(z, s$config)[[target]]$value, x)
    expect_true(numeric$valid)
    expect_equal(targets[[target]]$jacobian, numeric$jacobian, tolerance = 1e-8)
  }
  # Positive-definite joint covariance with nonzero between-owner covariance.
  covariance <- diag(length(x)) * .04
  coords <- build_param_slices(s$sizes)$log_slopes
  covariance[coords[1], coords[3]] <- covariance[coords[3], coords[1]] <- -.01
  a <- targets$components$jacobian
  v <- a %*% covariance %*% t(a)
  z <- targets$effective$jacobian
  product_covariance <- z %*% covariance %*% t(z)
  expect_equal(product_covariance[1, 1], v[1, 1] + v[4, 4] + 2 * v[1, 4])
  expect_equal(product_covariance[1, 1], .06)
  expect_gt(abs(product_covariance[1, 1] - (v[1, 1] + v[4, 4])), .01)
  changed <- s$config; changed$population_spec$active <- TRUE
  expect_error(mfrm_gpcm_product_slope_targets(x, changed), "fixed-standard-normal")
  expect_error(mfrm_gpcm_product_slope_targets(x[-1], s$config), "complete finite")
})

test_that("component tables and readiness distinguish shared labels across owners", {
  d <- gmfrm_test_data()
  d$Task <- paste0("level:", as.integer(d$Task))
  d$Rater <- paste0("level:", as.integer(d$Rater))
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(9L))
  common <- p$common
  par <- p$start + seq(-.2, .3, length.out = length(p$start))
  expanded <- expand_params(par, common$sizes, common$config)
  table <- build_slope_table(common$config, common$prep, expanded)
  expect_identical(table$SlopeOwner, rep(c("Task", "Rater"), each = 3L))
  expect_identical(table$SlopeFacet, rep(p$levels$Task, 2L))
  expect_equal(table$Estimate, unname(c(p$unpack(par)$slopes$Task, p$unpack(par)$slopes$Rater)))
  expect_equal(prod(table$Estimate[table$SlopeOwner == "Task"]), 1)
  expect_gt(abs(prod(table$Estimate[table$SlopeOwner == "Rater"]) - 1), .1)
  target <- mfrm_gpcm_product_slope_targets(par, common$config)$components
  expect_identical(table$ScaleReference, target$metadata$ScaleReference)
  expect_identical(table$Identification, target$metadata$Identification)
  expect_identical(table$ScaleReference, rep(c("geometric_mean_one", "fixed_standard_normal"), each = 3L))
  changed <- expanded
  changed$slope_components$Task <- changed$slope_components$Task[3:1, ]
  expect_error(build_slope_table(common$config, common$prep, changed), "level ordering")
  shuffled <- table[c(6, 1, 5, 2, 4, 3), ]
  record <- mfrmr_readiness_gpcm_slope_parameters(common$config, shuffled)
  expect_identical(record$Facet, table$SlopeOwner)
  expect_identical(record$Level, table$SlopeFacet)
  expect_equal(record$OptimizerEstimate, table$Estimate)
  expect_true(all(record$ParameterStatus == "not_evaluated"))
  expect_true(all(!record$SEEligible & !record$CIEligible))
  expect_true(all(is.na(record$PrimaryEstimate)))
  # A single-family boundary certificate must not promote the product model.
  altered <- common$config
  altered$boundary_audit$gpcm_slope_boundary <- list(state = "no_free_log_slope_coordinates")
  expect_equal(mfrmr_readiness_gpcm_slope_parameters(altered, shuffled), record)
  ready <- apply_mfrm_slope_readiness(shuffled, list(parameters = record[6:1, ]))
  expect_equal(ready$OptimizerEstimate, shuffled$Estimate)
  expect_identical(ready$ScaleReference, shuffled$ScaleReference)
  invalid_record <- record
  invalid_record$ParameterId[1] <- "Slope:Rater:wrong"
  expect_error(apply_mfrm_slope_readiness(shuffled, list(parameters = invalid_record)), "identities")
  expect_error(apply_mfrm_slope_readiness(shuffled[, -1], list(parameters = record)), "SlopeOwner")
  duplicate <- rbind(shuffled, shuffled[1, ])
  expect_error(mfrmr_readiness_gpcm_slope_parameters(common$config, duplicate), "unique within")
  # Updating interval fields must use both owner and level; locations stay untouched.
  ready$PrimaryEstimate <- seq_len(nrow(ready))
  ready$PrimaryLogEstimate <- log(ready$PrimaryEstimate)
  ready$CIEligible <- c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
  extra <- record[1, ]; extra$ParameterClass <- "facet_location"
  updated <- mfrm_update_slope_readiness_parameters(rbind(record, extra), ready)
  index <- match(seq_len(nrow(table)), c(6, 1, 5, 2, 4, 3))
  expect_equal(updated$PrimaryEstimate[1:6], ready$PrimaryEstimate[index])
  expect_identical(updated$CIEligible[1:6], ready$CIEligible[index])
  expect_identical(updated[7, ], rbind(record, extra)[7, ])
  expect_identical(mfrm_update_slope_readiness_parameters(extra, ready), extra)
})

test_that("joint observed information and Person scores use the product specification", {
  p <- gmfrm_test_problem()$problem
  s <- p$common
  s$config$estimation_control <- list(quad_points = 9L, mml_integration = "fixed")
  x <- unname(p$start + seq(-.2, .3, length.out = length(p$start)))
  retained <- list(prep = s$prep, config = s$config, opt = list(par = x))
  information <- compute_mml_parameter_covariance(retained)
  reference <- p$n_person * optimHess(x, function(z) p$marginal(z)$value,
                                     function(z) p$marginal(z)$gradient)
  expect_equal(information$hessian, reference, tolerance = 1e-7)
  scores <- mfrm_mml_person_scores_numeric(retained)
  expect_equal(unname(colSums(scores)), -p$n_person * p$marginal(x)$gradient,
               tolerance = 1e-7)
})

test_that("connected crossings alone do not identify marginal two-slope models", {
  # Hold all 16 task/rater rows fixed. Change only whether each Person has one
  # rating or four joint ratings. This is a local design check, not a recovery study.
  ranks <- integer(2)
  for (shared in c(FALSE, TRUE)) {
    d <- expand.grid(Task = c("t1", "t2"), Rater = c("r1", "r2"), Repeat = 1:4)
    d$Person <- if (shared) paste0("p", d$Repeat) else paste0("p", 1:16)
    d$Score <- rep(c(0L, 1L, 0L, 1L), 4)
    p <- mfrm_gmfrm_problem(d, 1L, gauss_hermite_normal(31L))
    s <- p$common
    x <- unname(p$start + seq(-.3, .4, length.out = length(p$start)))
    preflight <- audit_mfrm_estimability(s$prep, s$idx, s$config, s$sizes)
    expect_false(preflight$complete)
    expect_equal(preflight$design$nullity, 0L)
    expect_equal(preflight$readiness$OptimizerFreeDimension, 6L)
    audit <- audit_mfrm_mml_all_pattern_information(list(par = x), s$prep,
      s$idx, s$config, s$sizes, 31L, "log_slopes")
    expect_identical(audit$status, "evaluated_all_patterns_local_diagnostic_only")
    expect_lt(audit$probability_normalization$MaximumAbsoluteMassError, 1e-12)
    expect_lt(audit$expected_score_identity$MaximumAggregateScoreAbs, 1e-11)
    expect_lt(audit$numerical_differentiation$max_scaled_difference, 1e-8)
    expect_equal(audit$parameter_map, s$parameter_map)
    expect_false(audit$structural_identification_classified)
    expect_identical(audit$readiness_effect, "none_all_patterns_local_diagnostic_only")
    ranks[1L + shared] <- audit$local_rank
    observed <- audit_mfrm_mml_observed_pattern_score(list(par = x), s$prep,
      s$idx, s$config, s$sizes, 31L, "log_slopes")
    expect_identical(observed$status, "evaluated_observed_pattern_diagnostic_only")
    expect_lt(observed$reconstruction$ScoreSumVsNegativeGradientMaxAbs, 1e-10)
    expect_lt(observed$numerical_differentiation$max_scaled_difference, 1e-8)
    expect_false(observed$all_possible_response_patterns_evaluated)
    preflight$mml_all_pattern_information <- audit
    preflight$mml_observed_pattern_score <- observed
    classified <- mfrmr_classify_nonlinear_local_estimability(preflight, s$config)
    expect_identical(classified$state, if (shared) "locally_full_rank_sufficient" else
      "locally_first_order_rank_deficient")
    expect_true(classified$parameter_map_complete)
    expect_false(classified$global_identification_classified)
    expect_false(classified$continuous_integral_identification_classified)
    expect_false(classified$boundary_classified)
    expect_identical(classified$readiness_effect, "none_local_property_only")
  }
  expect_equal(ranks, c(4L, 6L))
})

test_that("two slope probabilities match the literal adjacent-category model", {
  p <- gmfrm_test_problem()$problem
  par <- p$start + seq(-.25, .3, length.out = length(p$start))
  u <- p$unpack(par)
  actual <- p$kernel(par)$prob
  expected <- do.call(rbind, lapply(p$nodes, function(theta) {
    t(vapply(seq_len(nrow(p$cells)), function(c) {
      i <- p$cells$Task[c]; r <- p$cells$Rater[c]
      adjacent <- u$slopes$Task[i] * u$slopes$Rater[r] *
        (theta - u$locations$Task[i] - u$locations$Rater[r] - u$steps[r, ])
      z <- c(0, cumsum(adjacent))
      exp(z - max(z)) / sum(exp(z - max(z)))
    }, numeric(3)))
  }))
  expect_equal(actual, unname(expected), tolerance = 1e-13)
  expect_equal(prod(u$slopes$Task), 1)
  expect_equal(sum(u$locations$Task), 0)
  expect_equal(unname(rowSums(u$steps)), rep(0, 3))
  expect_gt(abs(prod(u$slopes$Rater) - 1), .1)
})

test_that("M-step and observed marginal scores differentiate their own objectives", {
  p <- gmfrm_test_problem()$problem
  x <- p$start + seq(-.2, .2, length.out = length(p$start))
  counts <- p$marginal(p$start)$counts
  expect_equal(p$complete(x, counts)$gradient,
    gmfrm_numeric_gradient(function(v) p$complete(v, counts, FALSE)$value, x), tolerance = 1e-7)
  expect_equal(p$marginal(x)$gradient,
    gmfrm_numeric_gradient(function(v) p$marginal(v, FALSE)$value, x), tolerance = 1e-7)
  # Curvature of frozen Q cannot be used as observed-data information.
  observed <- stats::optimHess(x, function(v) p$marginal(v)$value,
    function(v) p$marginal(v)$gradient)
  counts <- p$marginal(x)$counts
  complete <- stats::optimHess(x, function(v) p$complete(v, counts)$value,
    function(v) p$complete(v, counts)$gradient)
  expect_gt(max(abs(observed - complete)), .01)
})

test_that("E step uses all ratings of a Person and is invariant to row order", {
  fixture <- gmfrm_test_problem()
  p <- fixture$problem
  x <- p$start + seq(-.2, .2, length.out = length(p$start))
  pr <- p$kernel(x)$prob
  cells <- p$cells
  d <- fixture$data
  cell <- match(paste(d$Task, d$Rater), paste(p$levels$Task[cells$Task], p$levels$Rater[cells$Rater]))
  likelihood <- vapply(p$levels$Person, function(id) {
    rows <- which(d$Person == id)
    sum(vapply(seq_along(p$nodes), function(q)
      prod(pr[cbind(cell[rows] + nrow(cells) * (q - 1L), d$Score[rows] + 1L)]) * p$weights[q], numeric(1)))
  }, numeric(1))
  expect_equal(p$marginal(x)$logLik, sum(log(likelihood)), tolerance = 1e-12)
  other <- mfrm_gmfrm_problem(d[nrow(d):1, ], 2L, gauss_hermite_normal(9L))
  expect_equal(other$marginal(x)$value, p$marginal(x)$value)
  expect_equal(other$marginal(x)$gradient, p$marginal(x)$gradient)
  expect_equal(other$marginal(x)$posterior, p$marginal(x)$posterior)
  expect_equal(unname(rowSums(p$marginal(x)$posterior)), rep(1, 4))
})

test_that("binary reduction and larger category sets retain valid gradients", {
  for (k in c(1L, 4L)) {
    p <- gmfrm_test_problem(k)$problem
    x <- p$start + seq(-.1, .1, length.out = length(p$start))
    expect_equal(p$marginal(x)$gradient,
      gmfrm_numeric_gradient(function(v) p$marginal(v)$value, x), tolerance = 1e-7)
    expect_equal(rowSums(p$kernel(x)$prob), rep(1, 81))
  }
})

test_that("unit slope reductions agree with the existing PCM and GPCM kernels", {
  p <- gmfrm_test_problem()$problem
  base <- p$start + seq(-.2, .2, length.out = length(p$start))
  for (fixed in list("log_slope:Task", "log_slope:Rater", "log_slope")) {
    par <- base
    par[grepl(fixed, names(par))] <- 0
    u <- p$unpack(par)
    cell <- rep(seq_len(nrow(p$cells)), length(p$nodes))
    i <- p$cells$Task[cell]; r <- p$cells$Rater[cell]
    eta <- rep(p$nodes, each = nrow(p$cells)) - u$locations$Task[i] - u$locations$Rater[r]
    cum <- t(apply(cbind(0, u$steps), 1L, cumsum))
    prob <- p$kernel(par)$log_prob
    for (k in 0:2) {
      if (fixed == "log_slope") {
        old <- loglik_pcm(eta, rep(k, length(eta)), cum, r)
      } else {
        slopes <- if (fixed == "log_slope:Task") u$slopes$Rater else u$slopes$Task
        owner <- if (fixed == "log_slope:Task") r else i
        old <- loglik_gpcm(eta, rep(k, length(eta)), cum, r, slopes, owner)
      }
      expect_equal(sum(prob[, k + 1L]), old, tolerance = 1e-12)
    }
  }
})

test_that("zero or overflowing effective slopes are not admitted as finite fits", {
  p <- gmfrm_test_problem()$problem
  for (extreme in c(-1000, 1000)) {
    par <- p$start
    par[grepl("log_slope:Rater", names(par))] <- extreme
    expect_null(p$kernel(par))
    expect_identical(p$marginal(par)$value, Inf)
  }
})

test_that("two-slope score sensitivity and information match probability derivatives", {
  fixture <- gmfrm_test_problem()
  p <- fixture$problem
  par <- p$start + seq(-.25, .3, length.out = length(p$start))
  prob <- p$kernel(par)$prob
  u <- p$unpack(par)
  slope <- rep(u$slopes$Task[p$cells$Task] * u$slopes$Rater[p$cells$Rater], length(p$nodes))
  score <- 0:2
  mean <- drop(prob %*% score)
  centered <- matrix(score, nrow(prob), length(score), byrow = TRUE) - mean
  variance <- rowSums(prob * centered^2)
  h <- 1e-5
  # Shift evaluation points solely to differentiate conditional probabilities;
  # do not use either shifted rule to evaluate the fixed-normal MML objective.
  shifted <- lapply(c(-h, h), function(shift) {
    q <- list(nodes = p$nodes + shift, weights = p$weights)
    mfrm_gmfrm_problem(fixture$data, 2L, q)$kernel(par)$prob
  })
  derivative <- (shifted[[2]] - shifted[[1]]) / (2 * h)
  expect_equal(derivative, unname(slope) * centered * prob, tolerance = 1e-8)
  expect_equal(drop(derivative %*% score), unname(slope) * variance, tolerance = 1e-8)
  expect_equal(rowSums(derivative^2 / prob), unname(slope)^2 * variance, tolerance = 1e-8)
})

test_that("unsupported inputs and disconnected task-rater blocks are explicit", {
  d <- gmfrm_test_data()
  quad <- gauss_hermite_normal(9L)
  expect_error(mfrm_gmfrm_problem(rbind(d, d[1, ]), 2L, quad), "Repeated")
  d$Score[1] <- NA
  expect_error(mfrm_gmfrm_problem(d, 2L, quad), "missing")
  d <- gmfrm_test_data()
  d <- d[as.integer(d$Task) == as.integer(d$Rater), ]
  expect_error(mfrm_gmfrm_problem(d, 2L, quad), "crossings")
})

test_that("EM traces ascent and does not claim success at the iteration limit", {
  p <- gmfrm_test_problem()$problem
  result <- mfrm_gmfrm_em(p, maxit = 2L, score_tol = 1e-12)
  expect_gte(min(diff(result$trace$logLik)), -1e-10)
  expect_gte(min(result$trace$q_gain[-1]), -1e-12)
  expect_false(result$converged)
  expect_identical(result$reason, "iteration_limit")
  expect_equal(result$parameter_map, p$common$parameter_map)
  expect_identical(result$nonlinear_transformation$status, "evaluated_diagnostic_only")
  expect_false(result$nonlinear_transformation$structural_identification_classified)
  expect_error(mfrm_gmfrm_em(p, maxit = 0L), "controls")
})

test_that("saved two-owner results replay in a fresh process without fitting", {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  skip_if_not_installed("callr")
  fixture <- gmfrm_test_problem()
  # Scrambled rows/levels make positional relabelling detectable on replay.
  data <- fixture$data[nrow(fixture$data):1, ]
  data$Task <- factor(data$Task, levels = c("t3", "t1", "t2"))
  data <- data[!(data$Task == "t3" & data$Rater == "r3"), ]
  names(data) <- c("Candidate", "観点 名", "Judge-ID", "Rating")
  owners <- c("観点 名", "Judge-ID")
  rule <- gauss_hermite_normal(7L)
  problem <- mfrm_gmfrm_problem(data, 2L, rule, slope_facets = owners,
    person = "Candidate", score = "Rating")
  fit <- mfrm_gmfrm_em(problem, maxit = 2L, score_tol = 1e-12)
  expect_false(fit$converged)
  expect_identical(fit$reason, "iteration_limit")
  # Keep all modelling columns and labels; expand.grid's display metadata is irrelevant.
  for (column in names(data)) expect_identical(fit$specification$data[[column]], data[[column]])
  expect_identical(fit$specification$quadrature, rule[c("nodes", "weights")])
  expect_identical(rownames(fit$posterior), problem$levels$Candidate)
  expect_identical(colnames(fit$posterior), paste0("node", seq_along(rule$nodes)))
  expect_identical(fit$controls$score_tol, 1e-12)
  # The missing crossing can be predicted using known levels, but stays labelled.
  grid <- data.frame(Theta = c(-.5, .7), Task = c("t1", "t3"), Rater = c("r1", "r3"))
  names(grid)[2:3] <- owners
  native <- mfrm_gmfrm_fit_result(problem, fit)
  path <- tempfile(fileext = ".rds"); withr::defer(unlink(path))
  saveRDS(list(numerical = fit, native = native), path)
  root <- normalizePath(find.package("mfrmr"))
  worker <- function(root, path, grid) {
    if (file.exists(file.path(root, "R", "core-gmfrm-em.R"))) {
      pkgload::load_all(root, quiet = TRUE, compile = FALSE)
    } else library("mfrmr", lib.loc = dirname(root))
    ns <- asNamespace("mfrmr")
    testthat::local_mocked_bindings(
      mfrm_gmfrm_em = function(...) stop("unexpected estimation"),
      fit_mfrm = function(...) stop("unexpected estimation"), .package = "mfrmr")
    saved <- readRDS(path)
    fit <- saved$numerical
    problem <- do.call(get("mfrm_gmfrm_problem", ns), fit$specification)
    state <- problem$marginal(fit$par)
    evaluator <- get("mfrm_gpcm_response_evaluator", ns)(saved$native$prep, saved$native$config, grid)
    list(fit = fit, native = saved$native,
      native_summary = summary(saved$native, compute = "never"),
      native_response = evaluator$evaluate(saved$native$opt$par), contexts = evaluator$contexts,
      logLik = state$logLik, posterior = state$posterior,
      max_score = max(abs(state$gradient)), parameters = problem$unpack(fit$par),
      parameter_map = problem$common$parameter_map,
      response = problem$response(fit$par, grid),
      targets = get("mfrm_gpcm_product_slope_targets", ns)(fit$par, problem$common$config))
  }
  environment(worker) <- baseenv()
  replay <- callr::r(worker, args = list(root, path, grid), libpath = .libPaths())
  expect_identical(replay$fit, fit)
  expect_identical(replay$native, native)
  expect_identical(replay$native_summary, summary(native, compute = "never"))
  expect_equal(replay$native_response$probabilities, problem$response(fit$par, grid)$probabilities)
  expect_identical(replay$contexts$ObservedContext, c(TRUE, FALSE))
  restored <- do.call(mfrm_gmfrm_problem, replay$fit$specification)
  restored_slopes <- build_slope_table(restored$common$config, restored$common$prep,
    expand_params(fit$par, restored$common$sizes, restored$common$config))
  restored_record <- mfrmr_readiness_gpcm_slope_parameters(restored$common$config, restored_slopes)
  expect_identical(restored_record, fit$slope_readiness)
  expect_identical(apply_mfrm_slope_readiness(restored_slopes, list(parameters = restored_record)), fit$slopes)
  expect_true(all(!fit$slopes$CIEligible))
  expect_equal(replay$logLik, fit$logLik, tolerance = 1e-12)
  expect_equal(replay$posterior, fit$posterior, tolerance = 1e-12)
  expect_equal(replay$max_score, fit$max_score, tolerance = 1e-12)
  expect_equal(replay$parameters, fit$parameters, tolerance = 1e-12)
  expect_identical(replay$parameter_map, fit$parameter_map)
  expect_equal(replay$response, problem$response(fit$par, grid), tolerance = 1e-12)
  expect_identical(replay$response$summary$ObservedContext, c(TRUE, FALSE))
  expect_equal(replay$targets, mfrm_gpcm_product_slope_targets(fit$par, problem$common$config))
  # Reconstructing predictions does not turn an unconverged fit into a qualified one.
  expect_false(replay$fit$converged)
  expect_identical(replay$fit$reason, "iteration_limit")
})

test_that("prediction and estimation share probabilities without merging slope owners", {
  p <- gmfrm_test_problem()$problem
  par <- p$start + seq(-.25, .3, length.out = length(p$start))
  cell <- rep(seq_len(nrow(p$cells)), length(p$nodes))
  grid <- data.frame(Theta = rep(p$nodes, each = nrow(p$cells)),
    Task = p$levels$Task[p$cells$Task[cell]], Rater = p$levels$Rater[p$cells$Rater[cell]])
  actual <- p$response(par, grid)
  expect_equal(unname(actual$probabilities), p$kernel(par)$prob, tolerance = 1e-13)
  expect_equal(unname(actual$log_probabilities), p$kernel(par)$log_prob, tolerance = 1e-13)
  expect_equal(actual$summary$EffectiveSlope, actual$summary$Slope1 * actual$summary$Slope2)
  expect_equal(actual$summary$Information,
    actual$summary$EffectiveSlope * actual$summary$ResponseSensitivity)
  expect_true(all(actual$summary$ObservedContext))
  reorder <- nrow(grid):1
  reverse <- p$response(par, grid[reorder, ])
  expect_equal(reverse$probabilities, actual$probabilities[reorder, ])
  expect_equal(reverse$summary$ExpectedScore, actual$summary$ExpectedScore[reorder])
  expect_error(p$response(par, transform(grid, Rater = "new")), "present in the fit")
  expect_error(p$response(par, transform(grid, Theta = NA_real_)), "finite numeric")
  expect_error(p$response(par, grid[, c("Theta", "Task")]), "every non-Person facet column")
})

test_that("known but unobserved crossings retain their extrapolation identity", {
  d <- gmfrm_test_data()
  d <- d[!(d$Task == "t3" & d$Rater == "r3"), ]
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(9L))
  par <- p$start + seq(-.25, .3, length.out = length(p$start))
  grid <- data.frame(Theta = c(.4, -.2), Task = c("t3", "t1"), Rater = c("r3", "r1"))
  actual <- p$response(par, grid)
  expect_identical(actual$summary$ObservedContext, c(FALSE, TRUE))
  u <- p$unpack(par)
  expected <- t(vapply(seq_len(nrow(grid)), function(row) {
    i <- grid$Task[row]; r <- grid$Rater[row]
    z <- c(0, cumsum(u$slopes$Task[i] * u$slopes$Rater[r] *
      (grid$Theta[row] - u$locations$Task[i] - u$locations$Rater[r] - u$steps[r, ])))
    exp(z - max(z)) / sum(exp(z - max(z)))
  }, numeric(3)))
  expect_equal(unname(actual$probabilities), unname(expected), tolerance = 1e-13)
  expect_equal(nrow(actual$summary), nrow(grid))
})

test_that("shared log probabilities and information remain stable in response tails", {
  p <- gmfrm_test_problem(1L)$problem
  grid <- data.frame(Theta = c(-1000, -40, 40, 1000), Task = "t1", Rater = "r1")
  actual <- p$response(p$start, grid)
  expect_true(all(is.finite(actual$log_probabilities)))
  expect_equal(actual$log_probabilities[1, 2], -1000)
  expect_equal(actual$log_probabilities[4, 1], -1000)
  reference <- exp(-40) / (1 + exp(-40))^2
  expect_equal(actual$summary$Information[2:3] / reference, c(1, 1), tolerance = 1e-12)
  expect_equal(actual$summary$ResponseSensitivity[2:3] / reference, c(1, 1), tolerance = 1e-12)
  for (model in c("RSM", "PCM", "GPCM")) {
    step <- if (model == "RSM") c(0, 0) else matrix(c(0, 0), 1)
    b <- mfrm_jml_probability_bundle(c(-1000, 1000), c(1L, 0L), model, step,
      criterion_idx = c(1L, 1L), slopes = 1, include_log_probs = TRUE)
    expect_equal(b$log_prob_obs, c(-1000, -1000))
    expect_equal(b$log_probs[cbind(1:2, c(2, 1))], b$log_prob_obs)
    empty <- mfrm_jml_probability_bundle(numeric(), integer(), model, step, include_log_probs = TRUE)
    expect_identical(dim(empty$log_probs), c(0L, 2L))
    default <- mfrm_jml_probability_bundle(c(-1000, 1000), c(1L, 0L), model, step,
      criterion_idx = c(1L, 1L), slopes = 1)
    expect_null(default$log_probs)
    expect_identical(default$probs, b$probs)
  }
})

test_that("shared response derivatives retain both slope families and joint covariance", {
  d <- gmfrm_test_data()
  d <- d[!(d$Task == "t3" & d$Rater == "r3"), ]
  names(d)[match(c("Task", "Rater"), names(d))] <- c("観点 名", "Judge-ID")
  p <- mfrm_gmfrm_problem(d, 2L, gauss_hermite_normal(9L),
    slope_facets = c("観点 名", "Judge-ID"))
  par <- p$start + seq(-.25, .3, length.out = length(p$start))
  grid <- data.frame(Theta = c(.4, -.2), check.names = FALSE)
  grid[["観点 名"]] <- c("t3", "t1"); grid[["Judge-ID"]] <- c("r3", "r1")
  evaluator <- mfrm_gpcm_response_evaluator(p$common$prep, p$common$config, grid)
  expect_identical(evaluator$newdata, grid)
  expect_identical(evaluator$contexts$ObservedContext, c(FALSE, TRUE))
  literal <- function(x, type) {
    u <- p$unpack(x)
    values <- lapply(seq_len(nrow(grid)), function(row) {
      i <- grid[["観点 名"]][row]; r <- grid[["Judge-ID"]][row]
      a <- u$slopes[[1]][i] * u$slopes[[2]][r]
      eta <- grid$Theta[row] - u$locations[[1]][i] - u$locations[[2]][r]
      z <- c(0, cumsum(a * (eta - u$steps[r, ])))
      pr <- exp(z - max(z)); pr <- pr / sum(pr)
      if (type == "probability") pr else
        unname(a^2 * sum(pr * ((0:2) - sum(pr * (0:2)))^2))
    })
    unlist(values, use.names = FALSE)
  }
  # A positive-definite covariance with explicit cross-family covariance.
  # This checks propagation, not a sampling or coverage claim.
  covariance <- diag(.01, length(par))
  first <- grep("^log_slope:観点 名:", names(par))[1]
  second <- grep("^log_slope:Judge-ID:", names(par))[1]
  covariance[first, second] <- covariance[second, first] <- .006
  for (type in c("probability", "information")) {
    evaluate <- function(x) {
      value <- evaluator$evaluate(x)
      if (type == "probability") as.vector(t(value$probabilities)) else value$information
    }
    jac <- mfrmr_numeric_transformation_jacobian(evaluate, par, relative_step = 1e-5)
    reference <- t(vapply(seq_along(evaluate(par)), function(row)
      gmfrm_numeric_gradient(function(x) literal(x, type)[row], par), numeric(length(par))))
    expect_true(jac$valid)
    expect_equal(evaluate(par), literal(par, type), tolerance = 1e-12)
    expect_equal(jac$jacobian, reference, tolerance = 1e-8, ignore_attr = TRUE)
    variance <- rowSums((jac$jacobian %*% covariance) * jac$jacobian)
    expect_equal(variance, diag(reference %*% covariance %*% t(reference)),
      tolerance = 1e-9, ignore_attr = TRUE)
    diagonal_only <- rowSums(sweep(jac$jacobian^2, 2L, diag(covariance), `*`))
    expect_gt(max(abs(variance - diagonal_only)), 1e-7)
    if (type == "probability") {
      expect_equal(colSums(jac$jacobian[1:3, , drop = FALSE]), rep(0, length(par)),
        tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
})
