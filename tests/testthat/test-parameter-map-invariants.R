.map_setup <- function(model = "RSM", method = "MML", facet = "Rater",
                       constrained = FALSE, population = FALSE) {
  d <- expand.grid(Person = paste0("p", 1:6), Rater = paste0("r", 1:4),
    Item = c("i1", "i2"), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d$Score <- rep(0:3, length.out = nrow(d))
  names(d)[names(d) == "Rater"] <- facet
  facets <- c(facet, "Item")
  prep <- prepare_mfrm_data(d, "Person", facets, "Score", keep_original = TRUE)
  interactions <- if (constrained) build_facet_interaction_specs(prep,
    list(facets), min_obs_per_interaction = 0) else list()
  pop <- if (population) list(active = TRUE, person_id = "Person",
    person_table = data.frame(Person = prep$levels$Person),
    design_matrix = cbind("(Intercept)" = 1, x = seq(-1, 1, length.out = prep$n_person))) else NULL
  built <- build_estimation_config(prep, model, method,
    if (model != "RSM") "Item", if (model == "GPCM") "Item",
    interaction_specs = interactions, weight_col = NULL,
    facet_signs = setNames(c(1, -1), facets), positive_facets = facet,
    noncenter_facet = "Person", dummy_facets = NULL,
    anchor_df = if (constrained) data.frame(Facet = facet, Level = "r1", Anchor = .7),
    group_anchor_df = if (constrained) data.frame(Facet = facet,
      Level = paste0("r", 1:3), Group = "g", GroupValue = .1),
    population = pop, centering_scope = "facet")
  cfg <- built$config
  if (constrained) cfg$step_specs <- lapply(cfg$step_specs, function(spec)
    build_step_constraint(3L, anchors = c("2" = .15), scope = spec$scope))
  sizes <- build_param_sizes(cfg)
  list(prep = prep, config = cfg, sizes = sizes,
    idx = build_indices(prep, cfg$step_facet, cfg$slope_facet, cfg$interaction_specs))
}

test_that("parameter-block names are refused before fitting and on saved configurations", {
  d <- load_mfrmr_data("example_core")
  cfg <- .map_setup()$config
  for (name in c("theta", "steps", "interactions", "log_slopes", "beta", "log_sigma2")) {
    changed <- d
    names(changed)[names(changed) == "Criterion"] <- name
    for (method in c("JML", "MML", "CML")) {
      expect_error(fit_mfrm(changed, "Person", c("Rater", name), "Score",
        method = method, category_policy = "preserve"), "reserved internal names")
    }
    expect_error(describe_mfrm_data(changed, "Person", c("Rater", name), "Score"),
      "reserved internal names")
    saved <- cfg; saved$facet_names[1] <- name
    expect_error(build_param_sizes(saved), "saved fits require refitting")
  }
})

test_that("similar facet names cannot create absent parameter blocks", {
  withr::local_options(warnPartialMatchDollar = TRUE)
  for (population in c(FALSE, TRUE)) for (model in if (isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) c("RSM", "PCM", "GPCM") else if (population) c("PCM", "GPCM") else "RSM") {
    b <- .map_setup(model, population = population)
    par <- .2 * sin(seq_len(sum(unlist(b$sizes))))
    reference <- mfrmr_make_adaptive_mml_evaluator(b$idx, b$config, b$sizes, 7L)(par)
    for (name in c("log_sigma2_group", "log_slopes_group", "beta_group", "interactions_group")) {
      a <- .map_setup(model, facet = name, population = population)
      expect_identical(unname(a$sizes), unname(b$sizes))
      map <- mfrmr_mml_optimizer_parameter_map(a$prep, a$idx, a$config, a$sizes)
      expect_identical(map$OptimizerIndex, seq_along(par))
      if (model != "GPCM") {
        adjacent <- mfrmr_estimability_adjacent_design(a$prep, a$idx, a$config, a$sizes)
        rank <- qr(as.matrix(adjacent$design))$rank
        quotient <- mfrmr_comparison_quotient(adjacent, a$sizes, rank)
        expect_identical(quotient$variance_index,
          build_param_slices(a$sizes)[["log_sigma2"]] %||% integer())
        expect_equal(ncol(quotient$transform), rank + as.integer(population))
      }
      value <- mfrmr_make_adaptive_mml_evaluator(a$idx, a$config, a$sizes, 7L)(par)
      expect_equal(value, reference, tolerance = 1e-12)
    }
  }
})

test_that("free-coordinate maps preserve constrained logits and their derivatives", {
  for (method in c("JML", "MML")) for (model in c("RSM", "PCM", "GPCM")) {
    b <- .map_setup(model, method, constrained = TRUE, population = method == "MML")
    cfg <- b$config; sizes <- b$sizes; idx <- b$idx
    par <- .2 * sin(seq_len(sum(unlist(sizes))))
    params <- expand_params(par, sizes, cfg)
    expect_length(build_initial_param_vector(cfg, sizes), length(par))
    expect_identical(unname(unlist(build_param_slices(sizes))), seq_along(par))
    expect_equal(collapse_expanded_params(params, cfg), par, tolerance = 1e-12)
    expect_equal(unname(params$facets$Rater[1]), .7)
    expect_equal(mean(params$facets$Rater[1:3]), .1, tolerance = 1e-12)
    expect_equal(sum(params$facets$Rater), 0, tolerance = 1e-12)
    step <- if (model == "RSM") matrix(params$steps, nrow = 1L) else params$steps_mat
    expect_equal(unname(step[, 2]), rep(.15, nrow(step)))
    expect_equal(unname(rowSums(step)), rep(0, nrow(step)), tolerance = 1e-12)
    expect_equal(rowSums(params$interactions[[1]]), rep(0, 4), ignore_attr = TRUE, tolerance = 1e-12)
    expect_equal(colSums(params$interactions[[1]]), rep(0, 2), ignore_attr = TRUE, tolerance = 1e-12)
    if (model == "GPCM") expect_equal(sum(params$log_slopes), 0, tolerance = 1e-12)

    # Independently assemble physical adjacent logits from labelled cells.
    adjacent <- function(x) {
      p <- expand_params(x, sizes, cfg)
      out <- matrix(0, nrow(b$prep$data), 3L)
      for (i in seq_len(nrow(out))) {
        row <- b$prep$data[i, ]
        eta <- p$facets$Rater[as.character(row$Rater)] - p$facets$Item[as.character(row$Item)] +
          p$interactions[[1]][as.character(row$Rater), as.character(row$Item)]
        if (method == "JML") eta <- eta + p$theta[as.character(row$Person)] else {
          person <- match(as.character(row$Person), b$prep$levels$Person)
          eta <- eta + sum(cfg$population_spec$design_matrix[person, ] * p$population$coefficients)
        }
        s <- if (model == "RSM") p$steps else p$steps_mat[as.character(row$Item), ]
        out[i, ] <- eta - s
      }
      as.vector(out)
    }
    design <- mfrmr_estimability_adjacent_design(b$prep, idx, cfg, sizes)
    at <- design$map$OptimizerIndex
    expect_equal(as.numeric(design$design %*% par[at]), adjacent(par) - adjacent(par * 0), tolerance = 1e-12)

    # Differentiate the direct likelihood in optimizer coordinates, including
    # fixed anchors, signed facets, interactions and typed step constraints.
    quad <- gauss_hermite_normal(21L)
    objective <- if (method == "JML") function(x) mfrm_loglik_jml(x, idx, cfg, sizes) else
      function(x) mfrm_loglik_mml(x, idx, cfg, sizes, quad)
    analytic <- if (method == "JML") mfrm_grad_jml(par, idx, cfg, sizes) else
      mfrm_grad_mml(par, idx, cfg, sizes, quad)
    numerical <- vapply(seq_along(par), function(j) {
      delta <- par * 0; delta[j] <- 1e-5
      (objective(par + delta) - objective(par - delta)) / 2e-5
    }, numeric(1))
    expect_equal(analytic, numerical, tolerance = 1e-7)
    expect_error(expand_params(c(par, 0), sizes, cfg), "free-coordinate blocks")
    expect_error(expand_params(par[-1L], sizes, cfg), "free-coordinate blocks")
  }
})

test_that("marginal score aggregation preserves weighted constrained derivatives", {
  for (model in c("RSM", "PCM", "GPCM")) {
    b <- .map_setup(model, "MML", constrained = TRUE, population = TRUE)
    idx <- b$idx
    idx$weight <- rep(c(0, .25, 1, 2.5), length.out = length(idx$score_k))
    par <- .2 * sin(seq_len(sum(unlist(b$sizes))))
    quad <- gauss_hermite_normal(21L)
    objective <- function(x) mfrm_loglik_mml(x, idx, b$config, b$sizes, quad)
    numerical <- vapply(seq_along(par), function(j) {
      delta <- par * 0
      delta[j] <- 1e-5
      (objective(par + delta) - objective(par - delta)) / 2e-5
    }, numeric(1))
    analytic <- mfrm_grad_mml(par, idx, b$config, b$sizes, quad)
    expect_equal(analytic, numerical, tolerance = 2e-7, info = model)
  }
})
