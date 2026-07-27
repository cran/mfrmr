test_that("packaged example datasets are available and well-formed", {
  keys <- list_mfrmr_data()
  expect_true(all(c("example_core", "example_bias", "example_operational") %in% keys))

  core <- load_mfrmr_data("example_core")
  bias <- load_mfrmr_data("example_bias")
  operational <- load_mfrmr_data("example_operational")

  for (dat in list(core, bias, operational)) {
    expect_true(is.data.frame(dat))
    expect_true(all(c("Study", "Person", "Rater", "Criterion", "Score", "Group") %in%
                      names(dat)))
    expect_gte(length(unique(dat$Person)), 30)
    expect_true(all(table(dat$Score) >= 10))
  }
})

test_that("example datasets are available through utils::data", {
  env <- new.env(parent = emptyenv())
  utils::data(list = "mfrmr_example_core", package = "mfrmr", envir = env)
  utils::data(list = "mfrmr_example_bias", package = "mfrmr", envir = env)
  utils::data(list = "mfrmr_example_operational", package = "mfrmr", envir = env)
  utils::data(list = "mfrmr_example_operational_design", package = "mfrmr", envir = env)

  expect_true(exists("mfrmr_example_core", envir = env, inherits = FALSE))
  expect_true(exists("mfrmr_example_bias", envir = env, inherits = FALSE))
  expect_true(exists("mfrmr_example_operational", envir = env, inherits = FALSE))
  expect_true(exists("mfrmr_example_operational_design", envir = env, inherits = FALSE))
  expect_s3_class(env$mfrmr_example_core, "data.frame")
  expect_s3_class(env$mfrmr_example_bias, "data.frame")
  expect_s3_class(env$mfrmr_example_operational, "data.frame")
  expect_s3_class(env$mfrmr_example_operational_design, "data.frame")
  expect_equal(nrow(env$mfrmr_example_operational_design), 288L)
  expect_false("Score" %in% names(env$mfrmr_example_operational_design))
})

test_that("operational example has a connected incomplete assignment", {
  dat <- load_mfrmr_data("example_operational")
  env <- new.env(parent = emptyenv())
  utils::data(
    list = "mfrmr_example_operational_design",
    package = "mfrmr",
    envir = env
  )
  planned <- env$mfrmr_example_operational_design

  expect_equal(nrow(dat), 282L)
  expect_equal(as.integer(table(dat$Rater)), c(47L, 56L, 50L, 47L, 44L, 38L))
  expect_true(all(table(dat$Criterion) == 94L))
  expect_true(all(table(dat$Person) %in% c(5L, 6L)))
  expect_equal(
    as.integer(table(unique(dat[c("Person", "Group")])$Group)),
    c(24L, 24L)
  )
  expect_equal(as.integer(table(dat$Group)), c(141L, 141L))
  expect_equal(sum(table(dat$Person) == 5L), 6L)
  key <- c("Person", "Rater", "Criterion")
  expect_equal(nrow(planned), 288L)
  expect_equal(nrow(dplyr::anti_join(planned[key], dat[key], by = key)), 6L)
  expect_equal(nrow(dplyr::anti_join(dat[key], planned[key], by = key)), 0L)

  person_rater <- unique(dat[c("Person", "Rater")])
  edges <- lapply(split(person_rater$Rater, person_rater$Person), unique)
  pair_key <- function(x) paste(sort(x), collapse = "|")
  observed_pair_keys <- vapply(edges, pair_key, character(1))
  planned_pair_keys <- vapply(seq_len(6L), function(i) {
    pair_key(sprintf("R%02d", c(i, (i %% 6L) + 1L)))
  }, character(1))
  expect_true(all(lengths(edges) == 2L))
  expect_true(all(observed_pair_keys %in% planned_pair_keys))
  expect_setequal(observed_pair_keys, planned_pair_keys)
  adjacency <- matrix(FALSE, nrow = 6L, ncol = 6L,
                      dimnames = list(sprintf("R%02d", 1:6), sprintf("R%02d", 1:6)))
  for (edge in edges) {
    adjacency[edge[1], edge[2]] <- TRUE
    adjacency[edge[2], edge[1]] <- TRUE
  }
  reached <- "R01"
  repeat {
    expanded <- union(reached, colnames(adjacency)[colSums(adjacency[reached, , drop = FALSE]) > 0])
    if (length(expanded) == length(reached)) break
    reached <- expanded
  }
  expect_setequal(reached, rownames(adjacency))
})

test_that("dataset catalog explains teaching roles", {
  catalog <- list_mfrmr_data(details = TRUE)

  expect_s3_class(catalog, "data.frame")
  expect_setequal(catalog$Key, list_mfrmr_data())
  expect_true("CountBasis" %in% names(catalog))
  expect_true(all(!catalog$Empirical))
  operational <- catalog[catalog$Key == "example_operational", , drop = FALSE]
  expect_equal(operational$PrimaryUse, "Beginner applied workflow")
  expect_match(operational$Design, "Connected two-rater assignment")

  combined <- catalog[catalog$Key == "combined", , drop = FALSE]
  expect_match(combined$PrimaryUse, "not direct fit", ignore.case = TRUE)
  expect_match(combined$Design, "explicit anchors/linking", ignore.case = TRUE)
  expect_match(combined$CountBasis, "513 persons and 30 raters", fixed = TRUE)
  expect_error(list_mfrmr_data(details = NA), "TRUE or FALSE", fixed = TRUE)
})

test_that("combined study objects expose their identity and linking hazard", {
  dat <- load_mfrmr_data("combined")
  ids_by_study <- split(dat, dat$Study)

  expect_gt(length(intersect(
    unique(ids_by_study[[1]]$Person),
    unique(ids_by_study[[2]]$Person)
  )), 0L)
  expect_gt(length(intersect(
    unique(ids_by_study[[1]]$Rater),
    unique(ids_by_study[[2]]$Rater)
  )), 0L)

  prefixed <- transform(
    dat,
    Person = paste(Study, Person, sep = "::"),
    Rater = paste(Study, Rater, sep = "::")
  )
  expect_length(unique(prefixed$Study), 2L)
  prefixed_raters <- split(prefixed$Rater, prefixed$Study)
  expect_length(intersect(
    unique(prefixed_raters[[1]]),
    unique(prefixed_raters[[2]])
  ), 0L)
})

test_that("example_bias supports non-null DIF and bias help examples", {
  dat <- load_mfrmr_data("example_bias")
  fit <- fit_mfrm(dat, "Person", c("Rater", "Criterion"), "Score",
                  method = "JML", maxit = 50)
  diag <- diagnose_mfrm(fit, residual_pca = "none")

  dif <- analyze_dif(fit, diag, facet = "Criterion",
                     group = "Group", data = dat, method = "residual")
  expect_true(max(abs(dif$dif_table$Contrast), na.rm = TRUE) > 0.3)

  bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion",
                        max_iter = 2)
  expect_true(max(abs(suppressWarnings(as.numeric(bias$table$`Bias Size`))), na.rm = TRUE) > 0.5)
  expect_true(all(as.character(bias$table$InferenceTier) == "screening"))
})
