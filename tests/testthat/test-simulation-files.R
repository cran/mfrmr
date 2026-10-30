simulation_file_test_args <- function() {
  list(n_person = 12, n_rater = 2, n_criterion = 2, reps = 2,
       maxit = 5, seed = 713)
}

simulation_file_without_time <- function(x) {
  if (is.data.frame(x)) return(x[, setdiff(names(x), "ElapsedSec"), drop = FALSE])
  if (is.list(x)) {
    x$reproducibility <- NULL
    return(lapply(x, simulation_file_without_time))
  }
  x
}

test_that("complete simulation files load the ordinary result without refitting", {
  path <- tempfile()
  on.exit(unlink(paste0(path, ".rds")), add = TRUE)
  args <- simulation_file_test_args()
  fit <- mfrmr::fit_mfrm
  calls <- 0L
  local_mocked_bindings(fit_mfrm = function(...) {
    calls <<- calls + 1L
    fit(...)
  }, .package = "mfrmr")
  set.seed(809)
  rng <- .Random.seed
  baseline <- suppressWarnings(do.call(simulate_mfrm_recovery, args))
  saved <- suppressWarnings(do.call(simulate_mfrm_recovery, c(args, list(file = path))))
  expect_identical(.Random.seed, rng)
  expect_equal(simulation_file_without_time(saved), simulation_file_without_time(baseline))
  before <- calls
  expect_message(again <- do.call(simulate_mfrm_recovery, c(args, list(file = path))),
    "Loaded completed simulation")
  expect_identical(calls, before)
  expect_identical(again, saved)
  expect_identical(readRDS(paste0(path, ".rds")), saved)
  expect_identical(saved$reproducibility$seeds, saved$rep_overview$Seed)
  expect_identical(saved$reproducibility$inputs$arguments$theta_sd, 1)
  expect_false(dir.exists(paste0(path, ".rds.lock")))

  for (change in list(list(seed = 714), list(reps = 3), list(maxit = 6),
                     list(theta_sd = 1.1))) {
    different <- utils::modifyList(args, change)
    expect_error(do.call(simulate_mfrm_recovery, c(different, list(file = path))),
      "does not match")
    expect_identical(readRDS(paste0(path, ".rds")), saved)
  }
  different <- utils::modifyList(args, list(reps = 1))
  expect_message(changed <- suppressWarnings(do.call(simulate_mfrm_recovery,
    c(different, list(file = path, file_refit = "on_change")))), "starting a new run")
  expect_equal(nrow(changed$rep_overview), 1)
  before <- calls
  suppressWarnings(do.call(simulate_mfrm_recovery,
    c(different, list(file = path, file_refit = "always"))))
  expect_equal(calls, before + 1L)
})

test_that("recovery resumes completed failures and retains the full trial RNG", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  args <- utils::modifyList(simulation_file_test_args(), list(reps = 3))
  write <- mfrmr:::simulation_file_write
  fit <- mfrmr::fit_mfrm
  writes <- 0L
  calls <- 0L
  draws <- numeric()
  first_seed <- mfrmr:::with_preserved_rng_seed(args$seed,
    sample.int(.Machine$integer.max, 3))[1]
  fail_draw <- mfrmr:::with_preserved_rng_seed(first_seed, runif(1))
  local_mocked_bindings(
    simulation_file_write = function(state) {
      writes <<- writes + 1L
      if (writes == 3L) stop("simulated interruption before second trial write")
      write(state)
    },
    fit_mfrm = function(...) {
      calls <<- calls + 1L
      draw <- runif(1)
      draws <<- c(draws, draw)
      if (identical(draw, fail_draw)) stop("recorded trial failure")
      fit(...)
    }, .package = "mfrmr")
  set.seed(1001)
  rng <- .Random.seed
  expect_error(suppressWarnings(do.call(simulate_mfrm_recovery,
    c(args, list(file = path)))), "simulated interruption")
  partial <- readRDS(path)
  expect_s3_class(partial, "mfrm_simulation_checkpoint")
  expect_identical(vapply(partial$cells, is.null, logical(1)), c(FALSE, TRUE, TRUE))
  expect_false(partial$cells[[1]]$rep_row$RunOK)
  expect_false(dir.exists(paste0(path, ".lock")))
  expect_identical(.Random.seed, rng)
  expect_message(resumed <- suppressWarnings(do.call(simulate_mfrm_recovery,
    c(args, list(file = path)))), "Resuming saved simulation")
  expect_equal(calls, 4)
  expect_identical(draws[2], draws[3])
  expect_identical(.Random.seed, rng)
  expect_identical(resumed$rep_overview$Error[1], "recorded trial failure")
  expect_equal(nrow(resumed$rep_overview), 3)
  baseline <- suppressWarnings(do.call(simulate_mfrm_recovery, args))
  expect_equal(simulation_file_without_time(resumed), simulation_file_without_time(baseline))
})

test_that("serial design files resume and future results reuse the same identity", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  plan <- future::plan(future::sequential)
  on.exit(future::plan(plan), add = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  args <- c(simulation_file_test_args(), list(progress = FALSE))
  write <- mfrmr:::simulation_file_write
  writes <- 0L
  local_mocked_bindings(simulation_file_write = function(state) {
    writes <<- writes + 1L
    if (writes == 3L) stop("design interruption")
    write(state)
  }, .package = "mfrmr")
  expect_error(suppressWarnings(do.call(evaluate_mfrm_design,
    c(args, list(file = path)))), "design interruption")
  partial <- readRDS(path)
  expect_identical(vapply(partial$cells, is.null, logical(1)), c(FALSE, TRUE))
  expect_message(resumed <- suppressWarnings(do.call(evaluate_mfrm_design,
    c(args, list(file = path)))), "Resuming saved simulation")
  baseline <- suppressWarnings(do.call(evaluate_mfrm_design, args))
  expect_equal(simulation_file_without_time(resumed), simulation_file_without_time(baseline))
  expect_message(again <- suppressMessages(do.call(evaluate_mfrm_design,
    c(args, list(file = path, parallel = "future")))), NA)
  expect_identical(again, resumed)
})

test_that("file errors cannot overwrite previous evidence or leave live locks", {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(path, paste0(path, ".lock")), recursive = TRUE), add = TRUE)
  open <- mfrmr:::simulation_file_open
  expect_error(open(character(), "never", "test", list(), 2, 1), "non-empty path")
  expect_error(open(file.path(tempfile(), "study"), "never", "test", list(), 2, 1),
    "directory does not exist")
  dir.create(paste0(path, ".lock"))
  expect_error(open(path, "always", "test", list(), 2, 1), "locked")
  unlink(paste0(path, ".lock"), recursive = TRUE)
  writeLines("incomplete RDS", path)
  before <- readLines(path)
  expect_error(open(path, "on_change", "test", list(), 2, 1), "Cannot read")
  expect_identical(readLines(path), before)
  expect_false(dir.exists(paste0(path, ".lock")))
  state <- open(path, "always", "test", list(a = 1), 2, NULL)
  unlink(state$lock, recursive = TRUE)
  expect_identical(readRDS(path)$reproducibility$seeds, state$seeds)
  expect_message(reopened <- open(path, "never", "test", list(a = 1), 2, NULL),
    "Resuming")
  expect_identical(reopened$seeds, state$seeds)
  unlink(reopened$lock, recursive = TRUE)
  withr::local_options(mfrmr.max_information_bytes = 1024)
  expect_error(open(path, "never", "test", list(a = 1), 2, NULL), "does not match")
  saved <- readRDS(path)
  saved$reproducibility$identity$engine$code <- "changed source"
  saveRDS(saved, path)
  expect_error(open(path, "never", "test", list(a = 1), 2, NULL), "does not match")
  expect_identical(readRDS(path), saved)
})

test_that("future checkpoints retain whole completed design rows", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  plan <- future::plan(future::sequential)
  on.exit(future::plan(plan), add = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  args <- c(utils::modifyList(simulation_file_test_args(), list(n_person = c(12, 14))),
    list(parallel = "future", progress = FALSE))
  write <- mfrmr:::simulation_file_write
  writes <- 0L
  local_mocked_bindings(simulation_file_write = function(state) {
    writes <<- writes + 1L
    if (writes == 3L) stop("second row interrupted before checkpoint")
    write(state)
  }, .package = "mfrmr")
  expect_error(suppressMessages(suppressWarnings(do.call(evaluate_mfrm_design,
    c(args, list(file = path))))), "second row interrupted")
  partial <- readRDS(path)
  expect_identical(vapply(partial$cells, is.null, logical(1)), c(FALSE, FALSE, TRUE, TRUE))
  resumed <- suppressMessages(suppressWarnings(do.call(evaluate_mfrm_design,
    c(args, list(file = path)))))
  expect_identical(resumed$rep_overview$ElapsedSec[1:2],
    vapply(partial$cells[1:2], function(x) x$rep_row$ElapsedSec, numeric(1)))
  baseline <- suppressMessages(suppressWarnings(do.call(evaluate_mfrm_design, args)))
  expect_equal(simulation_file_without_time(resumed), simulation_file_without_time(baseline))
})

test_that("design cells receive numerical options independently of worker defaults", {
  withr::local_options(mfrmr.use_cpp11_backend = FALSE)
  cell <- mfrmr:::simulation_evaluate_design_cell
  local_mocked_bindings(
    simulation_evaluate_design_cell = function(...) {
      old <- options(mfrmr.use_cpp11_backend = TRUE)
      on.exit(options(old), add = TRUE)
      result <- cell(...)
      expect_true(getOption("mfrmr.use_cpp11_backend"))
      result
    },
    fit_mfrm = function(...) {
      expect_false(getOption("mfrmr.use_cpp11_backend"))
      stop("intentional fit failure after checking worker options")
    }, .package = "mfrmr")
  value <- do.call(evaluate_mfrm_design, c(simulation_file_test_args(), list(progress = FALSE)))
  expect_true(all(grepl("checking worker options", value$rep_overview$Error)))
  expect_false(getOption("mfrmr.use_cpp11_backend"))
})
