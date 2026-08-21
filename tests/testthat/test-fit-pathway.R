test_that("fit pathway keeps fit on x and logits on y", {
  fit <- make_toy_fit(maxit = 20)
  diag <- make_toy_diagnostics(fit)
  out <- suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    diagnostics = diag,
    fit_stat = "Outfit",
    draw = FALSE
  ))

  expect_s3_class(out, "mfrm_plot_data")
  expect_identical(out$name, "fit_pathway")
  expect_identical(out$data$view, "fit_measure")
  expect_identical(out$data$fit_scale, "mnsq")
  expect_identical(out$data$fit_column, "Outfit")
  expect_equal(out$data$table$FitValue, out$data$table$Outfit)
  expect_equal(out$data$reference_lines$value, c(0.5, 1, 1.5))
  expect_false("Person" %in% out$data$table$Facet)
  expect_true(out$data$show_ci)
  expect_true(all(c(
    "CI_Lower", "CI_Upper", "SE_Method", "PrecisionTier", "CIBasis"
  ) %in% names(out$data$table)))
  expect_true(is.data.frame(out$data$uncertainty_basis))
  expect_true(all(c("WhiskerDefinition", "CI_Level") %in%
                    names(out$data$uncertainty_basis)))
})

test_that("fit pathway defaults to Infit on the horizontal axis", {
  fit <- make_toy_fit(maxit = 20)
  diag <- make_toy_diagnostics(fit)
  out <- suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    diagnostics = diag,
    draw = FALSE
  ))

  expect_identical(out$data$fit_stat, "Infit")
  expect_identical(out$data$fit_column, "Infit")
  expect_equal(out$data$table$FitValue, out$data$table$Infit)
})

test_that("fit pathway supports person selection and labels", {
  fit <- make_toy_fit(maxit = 20)
  diag <- make_toy_diagnostics(fit)
  person_ids <- as.character(diag$measures$Level[diag$measures$Facet == "Person"])[1:4]
  out <- suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    diagnostics = diag,
    include_person = TRUE,
    person_subset = person_ids,
    top_n_person = 2,
    person_labels = "all",
    facet_labels = "none",
    panel = "facet",
    draw = FALSE
  ))

  persons <- out$data$table[out$data$table$Facet == "Person", , drop = FALSE]
  expect_lte(nrow(persons), 2L)
  expect_true(all(persons$Level %in% person_ids))
  expect_true(all(nzchar(persons$LabelText)))
  expect_true(all(persons$Shape == 15L))
  expect_true(all(out$data$table$LabelText[out$data$table$Facet != "Person"] == ""))
  expect_identical(out$data$facet_labels, "none")
  expect_true(all(out$data$table$Panel == out$data$table$Facet))
})

test_that("fit pathway offers engine and FACETS companion ZSTD", {
  fit <- make_toy_fit(maxit = 20)
  facets <- suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    fit_scale = "zstd",
    zstd_method = "facets",
    draw = FALSE
  ))
  engine <- suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    fit_scale = "zstd",
    zstd_method = "engine",
    draw = FALSE
  ))

  expect_identical(facets$data$fit_column, "InfitZSTD_FACETS")
  expect_identical(engine$data$fit_column, "InfitZSTD_ENGINE")
  expect_equal(facets$data$table$FitValue, facets$data$table$InfitZSTD_FACETS)
  expect_equal(engine$data$table$FitValue, engine$data$table$InfitZSTD_ENGINE)
  expect_equal(facets$data$reference_lines$value, c(-2, 0, 2))
})

test_that("expected-score pathway remains distinct from fit pathway", {
  fit <- make_toy_fit(maxit = 20)
  expected <- suppressWarnings(plot(fit, type = "pathway", draw = FALSE))
  expect_identical(expected$name, "pathway_map")
  expect_match(expected$data$title, "Expected-score")
  expect_true(all(c("Theta", "ExpectedScore") %in% names(expected$data$expected)))
  expect_false("FitValue" %in% names(expected$data$expected))
})

test_that("fit pathway base renderer draws CI whiskers", {
  fit <- make_toy_fit(maxit = 20)
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  expect_silent(suppressWarnings(plot(
    fit,
    type = "fit_pathway",
    include_person = TRUE,
    top_n_person = 3
  )))
})

test_that("fit pathway label policies are validated", {
  fit <- make_toy_fit(maxit = 20)
  expect_error(
    plot(fit, type = "fit_pathway", facet_labels = "everything", draw = FALSE),
    "arg"
  )
})
