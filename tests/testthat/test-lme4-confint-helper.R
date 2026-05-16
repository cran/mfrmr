# Tests for the internal .lme4_confint_components() helper.
# This helper maps lme4::confint() row names to VarCorr component
# positions. lme4 uses two conventions for the SD components:
#   (a) terse: ".sig01"..".sigNN" plus ".sigma"
#   (b) verbose: "sd_(Intercept)|<group>" plus "sigma"
# The tests below enumerate both formats and a few edge cases, so
# that a future lme4 release which narrows or re-orders the row
# labels fails loudly here rather than silently returning NA_real_
# ICC CI bounds in compute_facet_icc(ci_method = "profile").

.lme4_confint_components <- getFromNamespace(".lme4_confint_components", "mfrmr")

test_that("terse '.sig01'/'.sigma' row names map to VarCorr positions", {
  # Random-intercept model with two RE terms and a residual.
  rn <- c(".sig01", ".sig02", ".sigma", "(Intercept)")
  vc_grp <- c("Person", "Rater", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(1L, 2L, 3L))
})

test_that("terse form handles three RE components plus residual", {
  rn <- c(".sig01", ".sig02", ".sig03", ".sigma", "(Intercept)")
  vc_grp <- c("Person", "Rater", "Criterion", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(1L, 2L, 3L, 4L))
})

test_that("verbose 'sd_...|<group>' row names map by group name", {
  rn <- c("sd_(Intercept)|Person", "sd_(Intercept)|Rater",
          "sd_(Intercept)|Criterion", "sigma", "(Intercept)")
  vc_grp <- c("Person", "Rater", "Criterion", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(1L, 2L, 3L, 4L))
})

test_that("verbose form tolerates a different component order in row names", {
  # VarCorr order is not always the same as confint() row order.
  rn <- c("sd_(Intercept)|Rater", "sd_(Intercept)|Person", "sigma")
  vc_grp <- c("Person", "Rater", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(2L, 1L, 3L))
})

test_that("missing residual row returns NA for the Residual position", {
  rn <- c(".sig01", ".sig02")
  vc_grp <- c("Person", "Rater", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(1L, 2L, NA_integer_))
})

test_that("empty row names or vc_grp return an all-NA integer", {
  expect_equal(.lme4_confint_components(character(0), character(0)),
               integer(0))
  expect_equal(.lme4_confint_components(character(0),
                                         c("Person", "Residual")),
               c(NA_integer_, NA_integer_))
  expect_equal(.lme4_confint_components(c(".sig01", ".sigma"),
                                         character(0)),
               integer(0))
})

test_that("group names containing regex metacharacters are matched literally", {
  rn <- c("sd_(Intercept)|Person.A", "sd_(Intercept)|Rater[1]", "sigma")
  vc_grp <- c("Person.A", "Rater[1]", "Residual")
  out <- .lme4_confint_components(rn, vc_grp)
  expect_equal(out, c(1L, 2L, 3L))
})

test_that("helper aligns with a real lme4 profile CI call", {
  skip_if_not_installed("lme4")
  toy <- load_mfrmr_data("example_core")
  for (col in c("Rater", "Criterion", "Person")) {
    toy[[col]] <- as.factor(as.character(toy[[col]]))
  }
  fit <- suppressWarnings(suppressMessages(
    lme4::lmer(Score ~ 1 + (1 | Person) + (1 | Rater) + (1 | Criterion),
               data = toy, REML = TRUE)
  ))
  ci <- suppressWarnings(suppressMessages(
    stats::confint(fit, level = 0.95, method = "profile")
  ))
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc <- vc[is.na(vc$var2), , drop = FALSE]
  ordered <- .lme4_confint_components(rownames(ci), as.character(vc$grp))
  # Every VarCorr component should resolve to a real row index.
  expect_false(any(is.na(ordered)))
  expect_true(all(ordered >= 1L & ordered <= nrow(ci)))
  # Row widths should be finite (not NA), confirming the mapping
  # picks actual SD rows rather than the fixed-effect intercept row.
  widths <- ci[ordered, 2L] - ci[ordered, 1L]
  expect_true(all(is.finite(widths) & widths > 0))
})

test_that("compute_facet_icc(ci_method = 'profile') uses the helper path", {
  skip_if_not_installed("lme4")
  toy <- load_mfrmr_data("example_core")
  icc <- compute_facet_icc(toy, facets = c("Rater", "Criterion"),
                            score = "Score", person = "Person",
                            ci_method = "profile", ci_level = 0.95)
  # CI columns populated for every VarCorr component.
  expect_true(all(c("ICC_CI_Lower", "ICC_CI_Upper", "ICC_CI_Method") %in%
                    names(icc)))
  expect_true(all(icc$ICC_CI_Method == "profile"))
  expect_true(all(is.finite(icc$ICC_CI_Lower)))
  expect_true(all(is.finite(icc$ICC_CI_Upper)))
  # Point estimates should lie inside their CIs (tolerance 1e-4).
  expect_true(all(icc$ICC_CI_Lower - 1e-4 <= icc$ICC &
                    icc$ICC <= icc$ICC_CI_Upper + 1e-4))
})
