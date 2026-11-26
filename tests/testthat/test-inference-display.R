test_that("tidy intervals follow core readiness and log-slope intervals", {
  d <- mfrmr::load_mfrmr_data("example_core")
  f <- mfrmr::fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = "GPCM", step_facet = "Criterion")
  for (component in c("facets", "persons", "steps", "slopes")) {
    core <- confint(f, component = if (component == "persons") "person" else component, se_type = "real")
    tb <- tidy(f, component)
    expect_equal(unname(as.matrix(tb[c("conf.low", "conf.high")])), unname(core))
  }
  f$converged <- FALSE
  for (component in c("facets", "persons", "steps", "slopes")) {
    tb <- tidy(f, component)
    expect_true(all(is.na(tb$conf.low) & is.na(tb$conf.high)))
    expect_true(all(tb$inference_status == "review_required"))
  }
})

test_that("FairM cutoffs invert weighted expectations for each person's task mix", {
  d <- mfrmr::load_mfrmr_data("example_core")
  d <- d[!(as.integer(sub("P", "", d$Person)) %% 2 == 0 & d$Criterion == "Content"), ]
  d$Weight <- rep(c(.5, 2), length.out = nrow(d))
  for (setting in c("RSM", "PCM", "full_predictor", "ability_only")) {
    model <- if (setting %in% c("RSM", "PCM")) setting else "GPCM"
    f <- mfrmr::fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score", model = model,
      step_facet = if (model == "RSM") NULL else "Criterion", weight = "Weight",
      spec = if (model == "GPCM") list(slope_type = setting) else list())
    # Independent matrix softmax, with raw data weights rather than encoded mixtures.
    offset <- sum(vapply(split(f$facets, f$facets$Facet), function(tb) tb$Sign[1] * mean(tb$Estimate), numeric(1)))
    a <- exp(weighted.mean(log(f$alpha_obs), d$Weight))
    means <- function(th, id) {
      ix <- which(as.character(d$Person) == id)
      groups <- if (model == "RSM") rep(1L, length(ix)) else match(d$Criterion[ix], f$encoding$step_owners)
      cum <- f$step_cumulative[groups, , drop = FALSE]
      eta <- if (setting == "ability_only") a * th + offset else a * (th + offset)
      logits <- matrix(rep(0:3 * eta, each = length(ix)), nrow = length(ix)) -
        if (setting == "ability_only") cum else a * cum
      probs <- exp(logits - apply(logits, 1, max))
      probs <- probs / rowSums(probs)
      weighted.mean(as.vector(probs %*% (1:4)), d$Weight[ix])
    }
    reconstructed <- vapply(seq_len(nrow(f$person)), function(i) means(f$person$Estimate[i], f$person$Person[i]), numeric(1))
    expect_equal(reconstructed, f$person$FairM, tolerance = 1e-10)
    cuts <- vapply(f$person$Person, function(id) uniroot(function(th) means(th, id) - 2.5, c(-20, 20), tol = 1e-10)$root, numeric(1))
    out <- mfrmr:::.tidymfrm_cut_decisions(f, 2.5)
    expect_equal(out$CutLogitMin, min(cuts), tolerance = 1e-8)
    expect_equal(out$CutLogitMax, max(cuts), tolerance = 1e-8)
    p <- pnorm((f$person$Estimate - cuts) / f$person$ModelSE)
    expect_true(is.na(out$Borderline90_ModelSE))
    expect_equal(out$Available_ModelSE, 0L)
    expect_match(out$InferenceReason, "nonunit_observation_weights")
    # Verify the normal-probability arithmetic separately from its eligibility.
    arithmetic <- mfrmr:::.tidymfrm_cut_decisions(f, 2.5,
      readiness = list(ready = TRUE, reason_codes = character()))
    expect_equal(arithmetic$Borderline90_ModelSE, sum(p > .05 & p < .95))
    expect_equal(arithmetic$Available_ModelSE, nrow(f$person))
    if (model != "RSM") expect_true(is.na(out$CutLogit))
    outside <- mfrmr:::.tidymfrm_cut_decisions(f, 4)
    expect_true(is.na(outside$Borderline90_ModelSE))
    expect_equal(outside$Available_ModelSE, 0)
    expect_match(outside$InferenceReason, "outside")
  }
})

test_that("saved workflows preserve unavailable screens in tables, plots and HTML", {
  d <- expand.grid(Person = paste0("P", 1:20), Rater = paste0("R", 1:3), stringsAsFactors = FALSE)
  d$Score <- ifelse(d$Rater == "R1", 3, rep(0:2, length.out = nrow(d)))
  wf <- mfrm_workflow(d, "Person", "Rater", "Score")
  wf <- suppressWarnings(step_mfrm_measure(wf, bootstrap = 2L))
  expect_null(summary(wf)$bootstrap)
  expect_match(summary(wf)$bootstrap_reason, "source calibration requires review")
  wf <- step_mfrm_diagnose(wf, cut_scores = 1.5)
  expect_true(all(tidy(wf)$flag == "Review required"))
  expect_true(all(is.na(tidy(wf)$conf.low)))
  expect_true(is.na(glance(wf)$flagged_raters))
  expect_equal(glance(wf)$unavailable_raters, 3)
  # A saved pre-fix workflow is rechecked, including stale Acceptable labels.
  wf$diagnostics$readiness <- NULL
  wf$diagnostics$raters$DiagnosticFlag <- "Acceptable"
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path))
  saveRDS(wf, path); wf <- readRDS(path)
  expect_true(all(summary(wf)$rater_diagnostics$DiagnosticFlag == "Review required"))
  expect_true(is.na(summary(wf)$decisions$Borderline90_ModelSE))
  p <- plot(wf, "raters", draw = FALSE)
  expect_true(all(is.na(p$data$RealLow)))
  expect_true(all(p$data$DiagnosticFlag == "Review required"))
  html <- render_mfrm_report(wf, format = "html")
  expect_match(html, "structural_likelihood_boundary")
  expect_match(html, "Unavailable rater screens")
})

test_that("undefined G coefficients never become recommended designs", {
  d <- mfrmr::load_mfrmr_data("example_core"); d$Score <- 2
  wf <- step_mfrm_gtheory(mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score", rating_min = 1, rating_max = 4))
  expect_equal(nrow(wf$gtheory$recommended_for_g), 0L)
  expect_true(all(is.na(tidy(wf, "d_study")$G)))
  expect_true(is.na(glance(wf)$g_coef))
  expect_s3_class(plot(wf, "dstudy", draw = FALSE), "ggplot")
  html <- render_mfrm_report(wf, "html")
  expect_match(html, "unavailable Phi (Zero total variance)", fixed = TRUE)
  expect_match(html, "Attainment of the 0.80 target cannot be assessed.", fixed = TRUE)
  expect_false(grepl("evaluable plans reaches Phi", html, fixed = TRUE))
})

test_that("reports preserve equal-count plans and original rating row numbers", {
  d <- mfrmr::load_mfrmr_data("example_core")
  wf <- tidymfrm(d, "Person", c("Rater", "Criterion"), "Score")
  md <- render_mfrm_report(wf, "markdown")
  expect_match(md, "Rater 4 x Criterion 3, projected Phi = 0.831", fixed = TRUE)
  expect_match(md, "| 3 | 4 | 12.000 | 0.860 | 0.828 | Yes | Yes |", fixed = TRUE)

  d$Score[2] <- NA_real_
  wf_csv <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score") |>
    step_mfrm_measure() |>
    step_mfrm_diagnose()
  u <- wf_csv$diagnostics$unexpected[1, ]
  expect_gt(u$.row_id, 2)
  expect_equal(u$Observed, d$Score[u$.row_id])
  expect_equal(u$Expected, fitted(wf_csv$fit, type = "response")[u$.row_id - 1L])
  html <- render_mfrm_report(wf_csv, "html", max_unexpected = 1L)
  expect_match(html, "<th>InputRow</th>", fixed = TRUE)
  expect_match(html, sprintf("<td>%d</td><td>%s</td><td>%s</td>", u$.row_id, u$Person, u$Rater), fixed = TRUE)
  expect_match(html, "No G-Study result is attached.", fixed = TRUE)
  planned <- step_mfrm_gtheory(wf_csv)
  reference <- mfrmr::mfrm_g_study(d[-2, ], scores = "Score", person = "Person", facets = c("Rater", "Criterion"), method = "urgenova")
  expect_equal(planned$gtheory$g_study_summary, reference$g_study_summary)
  expect_match(render_mfrm_report(planned), "The G-Study used 767 of 768 input rows.", fixed = TRUE)
})

test_that("unidentifiable planning retains calibration and rating review", {
  d <- mfrmr::load_mfrmr_data("example_core")
  lev <- unique(d$Criterion)
  d <- d[d$Criterion == lev[(match(d$Person, unique(d$Person)) - 1L) %% length(lev) + 1L], ]
  calibration <- mfrm_workflow(d, "Person", c("Rater", "Criterion"), "Score") |>
    step_mfrm_measure()
  expect_false(glance(calibration)$inference_ready)
  expect_true("population_assumption_linked" %in% summary(calibration$fit)$readiness$reason_codes)
  expect_error(step_mfrm_gtheory(calibration), "not identifiable")
  expect_warning(wf <- tidymfrm(d, "Person", c("Rater", "Criterion"), "Score"), "G-Study unavailable")
  wf <- unserialize(serialize(wf, NULL))
  expect_equal(coef(wf), coef(calibration))
  expect_null(summary(wf)$d_study)
  expect_match(summary(wf)$gtheory_reason, "not identifiable")
  expect_equal(nrow(summary(wf)$rater_diagnostics), 4L)
  html <- render_mfrm_report(wf, "html")
  expect_match(html, "Inference readiness is review_required (population_assumption_linked)", fixed = TRUE)
  expect_match(html, "G-Study unavailable:.*not identifiable")
  expect_error(tidymfrm(d, "Person", c("Rater", "Criterion"), "Score", g_facets = "Absent"),
               "Missing column")
})

test_that("conditional bootstrap refits use the original measurement specification", {
  d <- mfrmr::load_mfrmr_data("example_core")
  facets <- c("Rater", "Criterion")
  pd <- unique(d["Person"])
  pd$Cohort <- rep(c("A", "B"), length.out = nrow(pd))
  specs <- list(
    list(model = "GPCM", step_facet = "Criterion", slope_facet = "Rater",
         slope_type = "ability_only", gpcm_mml_identification = "fixed_standard_normal"),
    list(population_formula = ~ Cohort, person_data = pd, person_score = "WLE",
         mml_integration = "adaptive")
  )
  for (spec in specs) {
    set.seed(101)
    wf <- do.call(step_mfrm_measure, c(list(wf = mfrm_workflow(d, "Person", facets, "Score"),
                                            bootstrap = 3L), spec))
    probs <- fitted(wf$fit, type = "probs")
    set.seed(101)
    draws <- t(vapply(seq_len(nrow(probs)), function(i) {
      sample(as.integer(colnames(probs)), 3L, replace = TRUE, prob = probs[i, ])
    }, integer(3L)))
    fits <- lapply(seq_len(ncol(draws)), function(i) {
      simulated <- d
      simulated$Score <- draws[, i]
      do.call(mfrmr::fit_mfrm, c(list(data = simulated, person = "Person",
                                     facets = facets, score = "Score",
                                     rating_min = 1L, rating_max = 4L), spec))
    })
    usable <- vapply(fits, function(fit) isTRUE(summary(fit)$readiness$ready), logical(1))
    b <- summary(wf)$bootstrap
    expect_equal(b$BootstrapRequested, rep(3L, nrow(b)))
    expect_equal(b$BootstrapAvailable, rep(sum(usable), nrow(b)))
    expect_true(sum(usable) >= 2L)
    expected <- vapply(fits[usable], function(fit) fit$facets$Estimate, numeric(nrow(b)))
    expect_equal(b$BootstrapSE, apply(expected, 1L, sd), tolerance = 1e-7)
    expect_equal(b$BootLower95, apply(expected, 1L, quantile, probs = .025, names = FALSE),
                 tolerance = 1e-7)
  }

  wf$bootstrap$BootstrapAvailable <- NULL
  saved <- unserialize(serialize(wf, NULL))
  expect_null(summary(saved)$bootstrap)
  expect_match(render_mfrm_report(saved, "html"), "predates model-specification and refit checks")
})

test_that("G-Study uses the full declaration or an explicit subset", {
  d <- expand.grid(Person = 1:12, Rater = 1:2, Criterion = 1:2, Task = 1:2, Occasion = 1:2)
  d$Score <- pmin(4L, 1L + (d$Person + d$Rater + d$Criterion * d$Task) %% 3L +
    (d$Occasion - 1L))
  facets <- c("Rater", "Criterion", "Task", "Occasion")
  reordered <- facets[c(1L, 2L, 4L, 3L)]
  a <- step_mfrm_gtheory(mfrm_workflow(d, "Person", facets, "Score"))
  b <- step_mfrm_gtheory(mfrm_workflow(d, "Person", reordered, "Score"))
  expect_equal(a$gtheory$facets, facets)
  expect_equal(b$gtheory$facets, reordered)
  column_names <- c(facets, "G", "Phi")
  a_plan <- as.data.frame(summary(a)$d_study[column_names])
  b_plan <- as.data.frame(summary(b)$d_study[column_names])
  expect_equal(a_plan[do.call(order, a_plan[facets]), , drop = FALSE],
               b_plan[do.call(order, b_plan[facets]), , drop = FALSE],
               ignore_attr = TRUE, tolerance = 1e-10)
  selected <- tidymfrm(d, "Person", facets, "Score", g_facets = facets[1L:2L])
  expect_equal(selected$gtheory$facets, facets[1L:2L])
})
