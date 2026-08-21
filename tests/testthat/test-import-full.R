# Round-trip tests for the full mirt / TAM import path.
#
# - Each test fits the upstream package on synthetic data, calls the
#   importer with `compute_fit = TRUE`, and asserts that the
#   measurement-side bundle is populated and that the synthetic
#   `mfrm_diagnostics` slot has the canonical column shape.
# - Tests skip gracefully when the upstream package is not
#   installed; this keeps R CMD check passing on any reviewer's
#   machine.

skip_if_no <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    skip(paste0("`", pkg, "` (Suggests) not installed."))
  }
}

# --- mirt full import ---------------------------------------------------

test_that("import_mirt_fit returns a populated bundle", {
  skip_if_no("mirt")
  set.seed(42L)
  n_p <- 80L; n_i <- 6L
  theta <- stats::rnorm(n_p)
  beta <- seq(-1, 1, length.out = n_i)
  Y <- matrix(0L, n_p, n_i)
  for (i in seq_len(n_i)) {
    Y[, i] <- stats::rbinom(n_p, 1, stats::plogis(theta - beta[i]))
  }
  colnames(Y) <- paste0("I", seq_len(n_i))
  mirt_fit <- suppressMessages(suppressWarnings(
    mirt::mirt(as.data.frame(Y), 1, itemtype = "Rasch", verbose = FALSE)
  ))
  imp <- import_mirt_fit(mirt_fit, model = "RSM", compute_fit = TRUE)
  expect_s3_class(imp, "mfrm_imported_fit")
  expect_s3_class(imp, "mfrm_fit")
  expect_equal(nrow(imp$facets$person), n_p)
  expect_equal(nrow(imp$facets$others), n_i)
  expect_true("Slope" %in% names(imp$facets$others))
  expect_true(all(c("Infit", "Outfit") %in% names(imp$facets$others)))
  expect_s3_class(imp$diagnostics, "mfrm_diagnostics")
  expect_true(all(c("Facet", "Level", "Estimate") %in%
                    names(imp$diagnostics$measures)))
})

test_that("import_mirt_fit honours compute_fit = FALSE (skeleton)", {
  skip_if_no("mirt")
  set.seed(43L)
  Y <- matrix(stats::rbinom(80L * 4L, 1, 0.5), nrow = 80L)
  colnames(Y) <- paste0("I", seq_len(ncol(Y)))
  mirt_fit <- suppressMessages(suppressWarnings(
    mirt::mirt(as.data.frame(Y), 1, itemtype = "Rasch", verbose = FALSE)
  ))
  imp <- import_mirt_fit(mirt_fit, model = "RSM")
  expect_null(imp$diagnostics)
  expect_false("Infit" %in% names(imp$facets$others))
})

# --- TAM single-facet ---------------------------------------------------

test_that("import_tam_fit (single-facet) returns a populated bundle", {
  skip_if_no("TAM")
  set.seed(44L)
  n_p <- 80L; n_i <- 6L
  theta <- stats::rnorm(n_p)
  beta <- seq(-1, 1, length.out = n_i)
  Y <- matrix(0L, n_p, n_i)
  for (i in seq_len(n_i)) {
    Y[, i] <- stats::rbinom(n_p, 1, stats::plogis(theta - beta[i]))
  }
  colnames(Y) <- paste0("I", seq_len(n_i))
  tam_fit <- suppressMessages(suppressWarnings(
    TAM::tam.mml(as.data.frame(Y), verbose = FALSE)
  ))
  imp <- import_tam_fit(tam_fit, model = "RSM", compute_fit = TRUE)
  expect_s3_class(imp, "mfrm_imported_fit")
  expect_equal(nrow(imp$facets$person), n_p)
  expect_equal(nrow(imp$facets$others), n_i)
  expect_s3_class(imp$diagnostics, "mfrm_diagnostics")
  expect_identical(imp$summary$ICStatus[1], "imported_native_descriptive")
  expect_false(imp$summary$ICEligible[1])
  expect_false(imp$summary$ICSelectable[1])
  expect_true(imp$summary$NativeAICFormulaVerified[1])
  expect_true(imp$summary$NativeBICFormulaVerified[1])
  expect_true(imp$summary$NativeABICFormulaVerified[1])
  expect_true(imp$summary$NativeLogLikDevianceConsistent[1])
  expect_true(imp$summary$NativePersonCountConsistent[1])
  expect_true(is.finite(imp$summary$Iterations[1]))
  expect_true(is.finite(imp$summary$IterationCeiling[1]))
  expect_identical(
    imp$summary$Iterations[1],
    imp$source$convergence$Iterations
  )
  expect_identical(imp$source$dimensions, 1L)
  expect_identical(
    imp$summary$NativeABICFormula[1],
    "tam_deviance_plus_log_n_minus_2_over_24_k"
  )
})

# --- TAM multi-facet (tam.mml.mfr) --------------------------------------

test_that("import_tam_fit detects multi-facet fits", {
  skip_if_no("TAM")
  # tam.mml.mfr requires a wide layout with rater facets. Build a
  # tiny rater x item dataset.
  set.seed(45L)
  n_p <- 30L; n_r <- 3L; n_i <- 4L
  d <- expand.grid(pid = paste0("P", seq_len(n_p)),
                    rater = paste0("R", seq_len(n_r)),
                    stringsAsFactors = FALSE)
  d$pid <- as.character(d$pid)
  resp_mat <- matrix(0L, nrow(d), n_i)
  theta <- stats::rnorm(n_p)
  delta <- stats::rnorm(n_r, sd = 0.5)
  beta <- seq(-1, 1, length.out = n_i)
  for (i in seq_len(n_i)) {
    eta <- theta[match(d$pid, paste0("P", seq_len(n_p)))] -
            delta[match(d$rater, paste0("R", seq_len(n_r)))] -
            beta[i]
    resp_mat[, i] <- stats::rbinom(nrow(d), 1, stats::plogis(eta))
  }
  colnames(resp_mat) <- paste0("I", seq_len(n_i))
  resp_df <- cbind(d, as.data.frame(resp_mat))

  facets_df <- data.frame(rater = d$rater, stringsAsFactors = FALSE)
  formulaA <- ~ item + rater
  tam_mfr <- tryCatch(
    suppressMessages(suppressWarnings(
      TAM::tam.mml.mfr(resp = resp_mat, facets = facets_df,
                        pid = d$pid, formulaA = formulaA,
                        verbose = FALSE)
    )),
    error = function(e) NULL
  )
  if (is.null(tam_mfr)) skip("TAM::tam.mml.mfr failed on the toy data.")

  imp <- import_tam_fit(tam_mfr, model = "RSM", compute_fit = FALSE)
  expect_s3_class(imp, "mfrm_imported_fit")
  expect_true(isTRUE(imp$source$multi_facet))
  expect_gte(length(unique(imp$facets$others$Facet)), 1L)
})

# --- eRm person-parameter schemas ---------------------------------------

test_that("import_erm_fit reads the current eRm person-parameter schema", {
  skip_if_no("eRm")
  set.seed(47L)
  response_matrix <- matrix(sample(0:2, 120L, replace = TRUE), nrow = 30L)
  colnames(response_matrix) <- paste0("I", seq_len(ncol(response_matrix)))
  rownames(response_matrix) <- paste0("Candidate", seq_len(nrow(response_matrix)))

  erm_fit <- suppressMessages(suppressWarnings(eRm::PCM(response_matrix)))
  imported <- import_erm_fit(erm_fit, model = "PCM")

  expect_s3_class(imported, "mfrm_imported_fit")
  expect_identical(nrow(imported$facets$person), nrow(response_matrix))
  expect_identical(
    imported$facets$person$Person,
    rownames(eRm::person.parameter(erm_fit)$theta.table)
  )
  expect_true(all(is.finite(imported$facets$person$Estimate)))
  expect_true(all(is.finite(imported$facets$person$SE)))
})

test_that("eRm person extraction rejects ambiguous or misaligned schemas", {
  ambiguous <- list(theta.table = data.frame(theta = 1:2, thetapar = 1:2))
  expect_error(
    mfrmr:::.erm_extract_person_table(ambiguous),
    "Expected exactly one person-estimate column"
  )

  misaligned <- list(
    theta.table = data.frame(theta = c(0.1, 0.2, 0.3)),
    se.theta = list(c(0.4, 0.5))
  )
  expect_error(
    mfrmr:::.erm_extract_person_table(misaligned),
    "3 person estimates but 2 standard errors"
  )
})

# --- diagnostics integration -------------------------------------------

test_that("synthetic diagnostics has the slots downstream helpers expect", {
  skip_if_no("mirt")
  set.seed(46L)
  Y <- matrix(stats::rbinom(80L * 4L, 1, 0.5), nrow = 80L)
  colnames(Y) <- paste0("I", seq_len(ncol(Y)))
  mirt_fit <- suppressMessages(suppressWarnings(
    mirt::mirt(as.data.frame(Y), 1, itemtype = "Rasch", verbose = FALSE)
  ))
  imp <- import_mirt_fit(mirt_fit, model = "RSM", compute_fit = TRUE)
  diag <- imp$diagnostics
  expect_s3_class(diag, "mfrm_diagnostics")
  expect_true(all(c("measures", "fit", "reliability", "facets_chisq",
                     "overall_fit") %in% names(diag)))
  expect_true("Person" %in% diag$measures$Facet)
})
