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
