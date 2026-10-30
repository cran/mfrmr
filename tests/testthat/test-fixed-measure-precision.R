test_that("fixed coordinates include direct, group and singleton constraints", {
  specs <- list(
    direct = build_facet_constraint(c("A", "B", "C"), anchors = c(A = .2)),
    group = build_facet_constraint(c("A", "B", "C", "D"),
      anchors = c(A = .2), groups = c(A = "g", B = "g"), group_values = c(g = .5)),
    singleton = build_facet_constraint("A"),
    free = build_facet_constraint(c("A", "B"), centered = FALSE),
    all = build_facet_constraint(c("A", "B"), anchors = c(A = .2, B = -.3))
  )
  expected <- list(c(TRUE, FALSE, FALSE), c(TRUE, TRUE, FALSE, FALSE),
    TRUE, c(FALSE, FALSE), c(TRUE, TRUE))
  for (i in seq_along(specs)) {
    tbl <- data.frame(Facet = "Rater", Level = specs[[i]]$levels, SE = .1)
    out <- apply_fixed_measure_precision(tbl,
      list(method = "JML", facet_specs = list(Rater = specs[[i]])))
    expect_identical(out$Fixed, expected[[i]])
    expect_true(all(is.na(out$SE[out$Fixed])))
    expect_equal(out$SE[!out$Fixed], rep(.1, sum(!out$Fixed)))
  }
})

test_that("older measure tables retain their unconstrained precision calculations", {
  tbl <- tibble::tibble(Facet = "Rater", Level = c("A", "B"),
    Estimate = c(-1, 1), SE = c(.5, .5))
  expect_equal(calc_facets_chisq(tbl)$FixedChiSq, 8)
  expect_equal(summarize_precision_basis(tbl, "SE")$Reliability, .875)
})

test_that("PCM and GPCM MML separate fixed locations from derived-score uncertainty", {
  d <- load_mfrmr_data("example_core")
  for (model in c("PCM", "GPCM")) {
    fit <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
      method = "MML", model = model, step_facet = "Criterion",
      anchors = data.frame(Facet = "Rater", Level = "R01", Anchor = .1))
    dx <- diagnose_mfrm(fit)
    a <- dx$measures$Facet == "Rater" & dx$measures$Level == "R01"
    expect_true(dx$measures$Fixed[a])
    expect_true(is.na(dx$measures$SE[a]) && is.na(dx$measures$CI_Lower[a]))
    for (type in c("wright", "fit_pathway")) {
      if (model == "GPCM") {
        # This anchored GPCM has a source review restriction; drawing must retain it.
        expect_warning(p <- plot(fit, type = type, diagnostics = dx, show_ci = TRUE,
          draw = FALSE), "^Review-only display:")
      } else p <- plot(fit, type = type, diagnostics = dx, show_ci = TRUE, draw = FALSE)
      tbl <- if (type == "wright") p$data$locations else p$data$table
      a <- if (type == "wright") tbl$Group == "Rater" & tbl$Label == "R01" else
        tbl$Facet == "Rater" & tbl$Level == "R01"
      expect_equal(sum(a), 1L)
      expect_true(is.na(tbl$CI_Lower[a]) && is.na(tbl$CI_Upper[a]))
    }
    if (model == "GPCM") {
      fair <- fair_average_table(fit, dx, facets = "Rater", fair_se = TRUE)$raw_by_facet$Rater
      a <- fair$Level == "R01"
      expect_true(is.na(fair$ModelSE[a]))
      expect_gt(fair$FairMSE[a], 0)
      expect_true(is.finite(fair$FairM_CI_Lower[a]))
    }
  }
})

test_that("anchored JML and MML estimates never acquire diagnostic intervals", {
  d <- load_mfrmr_data("example_core")
  for (method in c("JML", "MML")) {
    anchors <- data.frame(Facet = "Rater", Level = "R01", Anchor = .1)
    if (method == "JML") anchors <- rbind(anchors,
      data.frame(Facet = "Person", Level = as.character(d$Person[1]), Anchor = .2))
    fit <- fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
      method = method, anchors = anchors)
    dx <- diagnose_mfrm(fit)
    fixed <- dx$measures[dx$measures$Fixed, ]
    expect_equal(nrow(fixed), nrow(anchors))
    expect_equal(unname(fixed$Estimate),
      anchors$Anchor[match(paste(fixed$Facet, fixed$Level), paste(anchors$Facet, anchors$Level))])
    for (nm in c("SE", "ModelSE", "RealSE", "CI_Lower", "CI_Upper"))
      expect_true(all(is.na(fixed[[nm]])))
    expect_false(any(fixed$CIEligible | fixed$SupportsFormalInference))
    expect_true(all(fixed$CIUse == "not_applicable"))
    expect_true(all(grepl("Fixed", fixed$CILabel)))
    free <- dx$measures$Facet == "Rater" & !dx$measures$Fixed
    expect_true(all(is.finite(dx$measures$SE[free])))
    expect_equal(dx$measures$CI_Lower[free],
      dx$measures$Estimate[free] - qnorm(.975) * dx$measures$SE[free])
    expect_true(all(is.na(dx$facets_chisq$FixedProb[dx$facets_chisq$Facet == "Rater"])))
    expect_true(all(is.na(dx$reliability$Reliability[dx$reliability$Facet == "Rater"])))
    expect_equal(dx$precision_review$Status[dx$precision_review$Check == "SE source labels"], "pass")

    for (saved in list(NULL, dx)) {
      se <- compute_se_for_plot(fit, ci_level = .9, diagnostics = saved)
      anchor <- se$Facet == "Rater" & se$Level == "R01"
      expect_true(se$Fixed[anchor])
      expect_true(is.na(se$SE[anchor]) && is.na(se$CI_Lower[anchor]))
    }
    p <- plot(fit, type = "wright", diagnostics = dx, show_ci = TRUE, draw = FALSE)
    expect_s3_class(p, "mfrm_plot_data")
    loc <- p$data$locations
    anchor <- loc$Group == "Rater" & loc$Label == "R01"
    expect_true(loc$Fixed[anchor])
    expect_true(is.na(loc$CI_Lower[anchor]) && is.na(loc$CI_Upper[anchor]))
    report <- fit_measures_table(fit, diagnostics = dx)
    anchor <- report$table$Facet == "Rater" & report$table$Level == "R01"
    expect_true(report$table$Fixed[anchor])
    expect_true(is.na(report$table$CI_Lower[anchor]))
    expect_true(report$facets_table$Fixed[report$facets_table$Facet == "Rater" &
      report$facets_table$Level == "R01"])
    # Saved diagnostic objects preserve both the values and their explanation.
    path <- tempfile(fileext = ".rds")
    saveRDS(dx, path)
    expect_identical(readRDS(path)$measures, dx$measures)
    unlink(path)
    attached <- attach_diagnostics_to_fit(fit)
    expect_true(is.na(attached$facets$others$SE[
      attached$facets$others$Facet == "Rater" & attached$facets$others$Level == "R01"]))
  }
})
