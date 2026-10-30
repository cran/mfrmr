# Independently assemble the declared equalities C b = d. Their null space
# must be exactly the free space used by expansion, gradients and covariance.
facet_equalities <- function(spec) {
  n <- length(spec$levels)
  rows <- list()
  target <- numeric(0)
  for (i in which(is.finite(spec$anchors))) {
    rows[[length(rows) + 1L]] <- as.numeric(seq_len(n) == i)
    target <- c(target, spec$anchors[i])
  }
  for (group in unique(stats::na.omit(spec$groups))) {
    members <- which(spec$groups == group)
    rows[[length(rows) + 1L]] <- as.numeric(seq_len(n) %in% members)
    target <- c(target, spec$group_values[[group]] * length(members))
  }
  if (isTRUE(spec$centered)) {
    rows[[length(rows) + 1L]] <- rep(1, n)
    target <- c(target, 0)
  }
  list(C = if (length(rows)) do.call(rbind, rows) else matrix(0, 0, n),
       d = unname(target))
}

test_that("whole-facet expansion spans exactly the declared affine space", {
  level <- LETTERS[1:5]
  declarations <- list(
    list(),
    list(anchors = c(A = .7)),
    list(anchors = c(E = .7)),
    list(anchors = c(A = .2, B = -.1, C = .3, D = .4)),
    list(anchors = c(A = .2, B = -.1, C = .3, D = .4, E = -.8)),
    list(anchors = c(A = .2), groups = c(A = "g", B = "g", C = "g"),
         group_values = c(g = .1)),
    list(groups = c(A = "g", B = "g", C = "h", D = "h", E = "h"),
         group_values = c(g = .3, h = -.2)),
    list(anchors = c(A = .2, B = -.2), groups = c(A = "g", B = "g"),
         group_values = c(g = 0)),
    list(centered = FALSE, anchors = c(A = .7)),
    list(centered = FALSE, groups = c(A = "g", B = "g"), group_values = c(g = .4))
  )
  for (args in declarations) {
    spec <- do.call(build_facet_constraint,
      c(list(levels = level, centering_scope = "facet"), args))
    eq <- facet_equalities(spec)
    J <- constraint_jacobian(spec)
    b0 <- expand_facet_with_constraints(numeric(spec$n_params), spec)
    expected_rank <- if (nrow(eq$C)) qr(t(eq$C))$rank else 0L
    expect_equal(spec$n_params, length(level) - expected_rank)
    expect_equal(drop(eq$C %*% b0), eq$d, tolerance = 1e-12)
    expect_equal(unname(eq$C %*% J), matrix(0, nrow(eq$C), spec$n_params),
                 tolerance = 1e-12)
    expect_equal(unname(as.matrix(mfrmr_constraint_jacobian_sparse(spec, "Rater")$jacobian)),
                 unname(J), tolerance = 1e-12)
    free <- seq_len(spec$n_params) / 7
    expanded <- expand_facet_with_constraints(free, spec)
    expect_equal(unname(expanded), unname(b0 + drop(J %*% free)), tolerance = 1e-12)
    expect_equal(collapse_facet_with_constraints(expanded, spec), free, tolerance = 1e-12)
    if (length(free)) {
      expect_equal(qr(J)$rank, length(free))
      # A quadratic in expanded coordinates has an independently known score.
      objective <- function(x) sum(expand_facet_with_constraints(x, spec)^2) / 2
      numerical <- vapply(seq_along(free), function(i) {
        delta <- rep(0, length(free)); delta[i] <- 1e-5
        (objective(free + delta) - objective(free - delta)) / 2e-5
      }, 0)
      expect_equal(constraint_grad_project(expanded, spec), numerical, tolerance = 1e-8)
      V <- J %*% diag(seq_along(free), nrow = length(free)) %*% t(J)
      expect_equal(unname(eq$C %*% V %*% t(eq$C)), matrix(0, nrow(eq$C), nrow(eq$C)),
                   tolerance = 1e-12)
    }
  }
})

test_that("the singleton remainder absorbs anchors rather than becoming zero", {
  spec <- build_facet_constraint(c("A", "B"), anchors = c(A = .7),
                                centering_scope = "facet")
  expect_equal(spec$n_params, 0)
  expect_equal(unname(expand_facet_with_constraints(numeric(0), spec)), c(.7, -.7))
  expect_error(collapse_facet_with_constraints(c(.7, 0), spec), "whole-facet")
})

test_that("affine offsets cannot create uncertainty in fixed constraint targets", {
  for (scope in c("facet", "remainder")) for (anchor in c(.4, 1.2, 1e16)) {
    spec <- build_facet_constraint(c("A", "B", "C"), anchors = c(A = anchor),
      centering_scope = scope)
    expect_identical(unname(constraint_jacobian(spec)), matrix(c(0, 1, -1), 3, 1))
    fit <- list(config = list(facet_names = "Rater", facet_specs = list(Rater = spec)))
    C <- rbind(Sum = c(1, 1, 1), Mean = c(1/3, 1/3, 1/3),
      Difference = c(0, 1, -1))
    colnames(C) <- spec$levels
    target <- mfrm_facet_contrasts(fit, "Rater", C)
    expect_identical(unname(target$fixed), c(TRUE, TRUE, FALSE))
    expect_identical(unname(target$jacobian), matrix(c(0, 0, 2), 3, 1))
  }
  # A fixed group mean leaves its two unanchored members free to differ.
  spec <- build_facet_constraint(c("A", "B", "C", "D"), anchors = c(A = .7),
    groups = c(A = "g", B = "g", C = "g"), group_values = c(g = .1),
    centering_scope = "facet")
  expect_identical(unname(constraint_jacobian(spec)), matrix(c(0, 1, -1, 0), 4, 1))
})

test_that("overlapping anchor messages do not imply their removal", {
  ratings <- load_mfrmr_data("example_core")
  review <- review_mfrm_anchors(ratings, "Person", c("Rater", "Criterion"), "Score",
    anchors = data.frame(Facet = "Rater", Level = "R01", Anchor = .7),
    group_anchors = data.frame(Facet = "Rater", Level = c("R01", "R02", "R03"),
      Group = "g", GroupValue = .1))
  expect_equal(nrow(review$anchors), 1L)
  expect_equal(nrow(review$group_anchors), 3L)
  expect_match(format_anchor_review_message(review), "retain both constraints")
  expect_match(format_anchor_review_message(review), "Invalid rows, if any")
})

test_that("whole-facet and group contradictions fail before optimization", {
  expect_error(build_facet_constraint(c("A", "B"), anchors = c(A = .7, B = .1),
    centering_scope = "facet"), "conflict with whole-facet")
  expect_error(build_facet_constraint(c("A", "B"), groups = c(A = "g", B = "g"),
    group_values = c(g = .2), centering_scope = "facet"), "conflict with whole-facet")
  expect_error(build_facet_constraint(c("A", "B", "C"), anchors = c(A = .7, B = .1),
    groups = c(A = "g", B = "g"), group_values = c(g = .2),
    centering_scope = "facet"), "conflict with the mean")
  expect_error(build_facet_constraint(c("A", "B", "C"),
    groups = c(A = "g", B = "g"), group_values = c(g = Inf),
    centering_scope = "facet"), "finite target")
  spec <- build_facet_constraint(c("A", "B"), anchors = c(A = .7, B = .1),
    centered = FALSE, centering_scope = "facet")
  expect_equal(unname(expand_facet_with_constraints(numeric(0), spec)), c(.7, .1))
})

test_that("legacy specifications retain their original expansion after serialization", {
  legacy <- list(levels = c("A", "B", "C"),
    anchors = c(A = .7, B = NA_real_, C = NA_real_),
    groups = c(A = NA_character_, B = NA_character_, C = NA_character_),
    group_values = numeric(0), centered = TRUE, n_params = 1)
  current <- build_facet_constraint(legacy$levels, anchors = c(A = .7))
  expect_identical(current, legacy)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(legacy, path)
  saved <- readRDS(path)
  expect_identical(expand_facet_with_constraints(.2, saved), c(A = .7, B = .2, C = -.2))
  whole <- build_facet_constraint(legacy$levels, anchors = c(A = .7), centering_scope = "facet")
  expect_equal(unname(expand_facet_with_constraints(.2, whole)), c(.7, .2, -.9))
  saveRDS(whole, path)
  expect_identical(expand_facet_with_constraints(.2, readRDS(path)),
                   expand_facet_with_constraints(.2, whole))
})

