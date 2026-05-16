test_that("read_facets_fit_table standardizes delimited FACETS-style columns", {
  path <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      Facet = c("Rater", "Rater"),
      Level = c("R1", "R2"),
      Measure = c(-0.25, 0.25),
      `S.E.` = c(0.11, 0.12),
      `Infit MnSq` = c(0.98, 1.05),
      `Outfit MnSq` = c(1.01, 1.07),
      `Infit ZStd` = c(-0.2, 0.6),
      `Outfit ZStd` = c(0.1, 0.7),
      `T.Count` = c(24, 25),
      `Infit df` = c(21.5, 22.5),
      `Outfit df` = c(23.5, 24.5),
      check.names = FALSE
    ),
    path,
    row.names = FALSE
  )

  out <- mfrmr::read_facets_fit_table(path)

  expect_s3_class(out, "tbl_df")
  expect_equal(out$Facet, c("Rater", "Rater"))
  expect_equal(out$Level, c("R1", "R2"))
  expect_equal(out$Estimate, c(-0.25, 0.25))
  expect_equal(out$SE, c(0.11, 0.12))
  expect_equal(out$N, c(24, 25))
  expect_equal(out$Infit, c(0.98, 1.05))
  expect_equal(out$Outfit, c(1.01, 1.07))
  expect_equal(out$InfitZSTD, c(-0.2, 0.6))
  expect_equal(out$OutfitZSTD, c(0.1, 0.7))
  expect_equal(out$DF_Infit, c(21.5, 22.5))
  expect_equal(out$DF_Outfit, c(23.5, 24.5))
})

test_that("read_facets_fit_table parses FACETS score.N.txt files", {
  dir <- tempfile()
  dir.create(dir)
  path <- file.path(dir, "score.2.txt")
  writeLines(c(
    "Some FACETS heading",
    "\"Measure\",\"S.E.\",\"Infit MS\",\"Infit Z\",\"Outfit MS\",\"Outfit Z\",\"T.Count\",\"Infit df\",\"Outfit df\",\"Rater\",\"F-Number\"",
    "0.10,0.05,0.99,-0.10,1.01,0.20,18,16.5,17.5,\"Rater_A\",2",
    "-0.20,0.06,1.04,0.50,1.06,0.70,19,17.2,18.1,\"Rater_B\",2"
  ), path)

  out <- mfrmr::read_facets_fit_table(
    path,
    facet_map = c("1" = "Person", "2" = "Rater")
  )

  expect_equal(out$Facet, c("Rater", "Rater"))
  expect_equal(out$Level, c("Rater_A", "Rater_B"))
  expect_equal(out$Estimate, c(0.10, -0.20))
  expect_equal(out$SE, c(0.05, 0.06))
  expect_equal(out$Infit, c(0.99, 1.04))
  expect_equal(out$Outfit, c(1.01, 1.06))
  expect_equal(out$InfitZSTD, c(-0.10, 0.50))
  expect_equal(out$OutfitZSTD, c(0.20, 0.70))
  expect_equal(out$DF_Infit, c(16.5, 17.2))
  expect_equal(out$DF_Outfit, c(17.5, 18.1))
  expect_equal(out$RawFacetNumber, c("2", "2"))
})

test_that("read_facets_fit_table parses fixed-field FACETS score files", {
  dir <- tempfile()
  dir.create(dir)
  path <- file.path(dir, "score.2.txt")
  fixed_row <- function(raw_score, count, measure, se, infit, infit_z,
                        outfit, outfit_z, infit_df, outfit_df, number, label) {
    chars <- rep(" ", 250 + nchar(label))
    put <- function(start, end, value) {
      width <- end - start + 1L
      txt <- sprintf(paste0("%", width, "s"), as.character(value))
      substring <- strsplit(txt, "", fixed = TRUE)[[1]]
      chars[start:end] <<- utils::tail(substring, width)
    }
    put(1, 10, raw_score)
    put(11, 20, count)
    put(41, 50, measure)
    put(51, 60, se)
    put(61, 70, infit)
    put(71, 80, infit_z)
    put(81, 90, outfit)
    put(91, 100, outfit_z)
    put(191, 200, infit_df)
    put(221, 230, outfit_df)
    put(241, 250, number)
    chars[251:(250 + nchar(label))] <- strsplit(label, "", fixed = TRUE)[[1]]
    paste(chars, collapse = "")
  }
  writeLines(c(
    "FACETS fixed-field score file",
    fixed_row("9.7", "18.0", "-0.38", "1.04", "1.00", "-0.2",
              "0.70", "1.0", "16.5", "17.5", "1", "Rater_A"),
    fixed_row("10.0", "19.0", "0.25", "0.06", "1.04", "0.5",
              "1.06", "0.7", "17.2", "18.1", "2", "Rater_B")
  ), path)

  out <- mfrmr::read_facets_fit_table(
    path,
    facet_map = c("1" = "Person", "2" = "Rater")
  )

  expect_equal(out$Facet, c("Rater", "Rater"))
  expect_equal(out$Level, c("Rater_A", "Rater_B"))
  expect_equal(out$Estimate, c(-0.38, 0.25))
  expect_equal(out$SE, c(1.04, 0.06))
  expect_equal(out$N, c(18, 19))
  expect_equal(out$Infit, c(1.00, 1.04))
  expect_equal(out$Outfit, c(0.70, 1.06))
  expect_equal(out$InfitZSTD, c(-0.2, 0.5))
  expect_equal(out$OutfitZSTD, c(1.0, 0.7))
  expect_equal(out$DF_Infit, c(16.5, 17.2))
  expect_equal(out$DF_Outfit, c(17.5, 18.1))
  expect_equal(out$RawFacetNumber, c("2", "2"))
})

test_that("imported FACETS fit table can feed facets_fit_review", {
  d <- mfrmr:::sample_mfrm_data(seed = 323)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none", fit_df_method = "both")
  facets_like <- diag$fit |>
    dplyr::transmute(
      Facet = .data$Facet,
      Level = .data$Level,
      Infit = .data$Infit,
      Outfit = .data$Outfit,
      InfitZSTD = .data$InfitZSTD_FACETS,
      OutfitZSTD = .data$OutfitZSTD_FACETS,
      DF_Infit = .data$DF_Infit_FACETS,
      DF_Outfit = .data$DF_Outfit_FACETS,
      `T.Count` = .data$N
    )
  path <- tempfile(fileext = ".csv")
  utils::write.csv(facets_like, path, row.names = FALSE)

  imported <- mfrmr::read_facets_fit_table(path)
  review <- mfrmr::facets_fit_review(fit, diagnostics = diag, facets_fit = imported)

  expect_gt(nrow(imported), 0)
  expect_true(all(review$external_comparison$ExternalMatched))
  expect_true(all(review$external_comparison$ExternalStatus == "same"))
  expect_equal(review$summary$ExternalNeedsReview, 0)
})

test_that("simulation-style partial FACETS tables import and report quality", {
  d <- mfrmr:::sample_mfrm_data(seed = 324)
  fit <- suppressWarnings(
    mfrmr::fit_mfrm(
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score",
      method = "JML",
      model = "RSM",
      maxit = 20
    )
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none", fit_df_method = "both")
  fit_tbl <- as.data.frame(diag$fit, stringsAsFactors = FALSE)
  partial <- fit_tbl |>
    dplyr::filter(.data$Facet == "Rater") |>
    dplyr::slice_head(n = 3) |>
    dplyr::transmute(
      Facet = .data$Facet,
      Level = .data$Level,
      InfitZSTD = .data$InfitZSTD_FACETS,
      OutfitZSTD = .data$OutfitZSTD_FACETS,
      TCount = .data$N
    )
  duplicate_row <- partial[1, , drop = FALSE]
  partial <- dplyr::bind_rows(partial, duplicate_row)
  path <- tempfile(fileext = ".csv")
  utils::write.csv(partial, path, row.names = FALSE)

  imported <- mfrmr::read_facets_fit_table(path)
  review <- mfrmr::facets_fit_review(fit, diagnostics = diag, facets_fit = imported)

  expect_equal(nrow(imported), 4L)
  expect_true(all(is.na(imported$Infit)))
  expect_true(all(is.na(imported$DF_Infit)))
  expect_equal(review$external_table_quality$Rows, 4L)
  expect_equal(review$external_table_quality$CompleteMnSqRows, 0L)
  expect_equal(review$external_table_quality$CompleteDFRows, 0L)
  expect_equal(review$external_table_quality$CompleteZSTDRows, 4L)
  expect_equal(review$external_table_quality$DuplicateFacetLevelRows, 2L)
  expect_equal(review$summary$ExternalDuplicateKeyRows, 2L)
  expect_gt(nrow(review$external_comparison), 0)
})
