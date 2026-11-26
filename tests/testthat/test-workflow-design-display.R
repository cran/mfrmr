test_that("anchored calibrations pass through a workflow with their retained rows", {
  d <- mfrmr::load_mfrmr_data("example_core")
  level <- as.character(d$Rater[1L])
  fit <- mfrmr::fit_mfrm(d, "Person", c("Rater", "Criterion"), "Score",
    anchors = data.frame(Facet = "Rater", Level = level, Anchor = 0.25))
  work <- function(data) mfrm_workflow(data, "Person", c("Rater", "Criterion"), "Score", fit = fit)
  wf <- step_mfrm_diagnose(work(d))
  expect_identical(wf$fit, fit)
  expect_equal(tidy(wf, "raters")$estimate[tidy(wf, "raters")$term == level], 0.25)
  expect_match(render_mfrm_report(wf), "Calibration includes 1 fixed elements", fixed = TRUE)
  expect_error(work(d[nrow(d):1L, ]), "original order")
  display <- plot(wf, type = "pathway", facet = "Rater", fit_statistic = "both", draw = FALSE)
  expect_setequal(display$data$Statistic, c("Infit", "Outfit"))
  expect_equal(display$data$FitMnSq[display$data$Statistic == "Outfit"],
               tidy(wf, "raters")$outfit_mnsq)
})

test_that("reports count absent and missing ratings against the planned roster", {
  ratings <- mfrmr::load_mfrmr_data("example_core")
  roster <- ratings[c("Person", "Rater", "Criterion")]
  ratings$Score[2L] <- NA
  wf <- mfrm_workflow(ratings[-1L, ], "Person", c("Rater", "Criterion"), "Score",
                      expected_design = roster)
  expect_equal(summary(wf)$planned_coverage$MissingExpectedCells, 2L)
  expect_match(render_mfrm_report(wf), "2 not observed", fixed = TRUE)
})

test_that("workflow connectivity uses the fitted category support", {
  set.seed(261126L)
  facets <- c("Rater", "Criterion")
  d <- expand.grid(Person = sprintf("P%02d", 1:30), Rater = paste0("R", 1:4),
                   Criterion = LETTERS[1:3], stringsAsFactors = FALSE)
  support <- data.frame(Level = LETTERS[1:3], Min = c(1, 3, 1), Max = c(4, 5, 4))
  owner <- match(d$Criterion, support$Level)
  d$Score <- mapply(function(lo, hi) sample(seq.int(lo, hi), 1L),
                    support$Min[owner], support$Max[owner])
  d$Score[d$Person == "P01"] <- support$Max[owner[d$Person == "P01"]]
  original <- mfrm_workflow(d, "Person", facets, "Score", rating_min = 1, rating_max = 5)
  wf <- step_mfrm_measure(original, model = "PCM", step_facet = "Criterion",
                          category_support = support)
  direct <- check_mfrm_connectivity(d, "Person", facets, "Score", model = "PCM",
    step_facet = "Criterion", category_support = support, rating_min = 1, rating_max = 5)
  extreme <- function(w) {
    persons <- summary(w)$connectivity$person_summary
    unlist(persons[persons$Person == "P01", c("RawScore", "MaxScore", "ExtremeFlag")])
  }
  expect_equal(extreme(wf), c(RawScore = 52, MaxScore = 52, ExtremeFlag = 1))
  expect_equal(summary(wf)$connectivity$person_summary, direct$person_summary)
  expect_identical(wf$connectivity$rater_linkage, original$connectivity$rater_linkage)
  attached <- mfrm_workflow(d, "Person", facets, "Score", fit = wf$fit)
  expect_equal(extreme(attached), extreme(wf))
  legacy <- wf
  legacy$connectivity <- original$connectivity
  legacy$connectivity_calibrated <- NULL
  expect_equal(extreme(unserialize(serialize(legacy, NULL))), extreme(wf))
  expect_equal(extreme(step_mfrm_measure(wf)), c(RawScore = 52, MaxScore = 60, ExtremeFlag = 0))
})

test_that("explicit over-limit G-Study retains the measurement workflow", {
  set.seed(6517L)
  d <- expand.grid(Person = 1:16, Rater = 1:2, Criterion = 1:2, Task = 1:2,
                   Site = 1:2, Mode = 1:2)
  d$Score <- sample(1:4, nrow(d), replace = TRUE)
  facets <- setdiff(names(d), c("Person", "Score"))
  expect_warning(
    wf <- mfrm_analysis(d, "Person", facets, "Score", g_facets = facets),
    "supports up to four non-person facets")
  expect_s3_class(wf$fit, "mfrm_fit")
  expect_equal(nrow(augment(wf)), nrow(d))
  expect_match(summary(wf)$gtheory_reason, "four non-person facets")
  expect_match(render_mfrm_report(wf, format = "html"), "G-Study unavailable", fixed = TRUE)
  expect_null(step_mfrm_gtheory(wf, g_facets = "Rater")$gtheory_reason)
  expect_error(step_mfrm_gtheory(wf, g_facets = c(facets, "Typo")), "Missing column.*Typo")
})

test_that("calibration and prediction reject ambiguous covariate IDs", {
  d <- mfrmr::load_mfrmr_data("example_core")
  facets <- c("Rater", "Criterion")
  persons <- unique(as.character(d$Person))
  raters <- unique(as.character(d$Rater))
  pd <- data.frame(Person = persons, Prior = seq_along(persons) / length(persons))
  rd <- data.frame(Rater = raters, Tenure = seq_along(raters))
  estimate <- function(pd, rd) fit_mfrm(d, "Person", facets, "Score",
    population_formula = ~ Prior, person_data = pd,
    facet_formulas = list(Rater = ~ Tenure), facet_data = list(Rater = rd))
  bad_p <- rbind(pd, transform(pd[1L, ], Prior = 9))
  bad_r <- rbind(rd, transform(rd[1L, ], Tenure = 9))
  expect_error(estimate(bad_p, rd), "person_data.*Duplicate ID")
  expect_error(estimate(pd, bad_r[nrow(bad_r):1L, ]), "facet_data.*Duplicate ID")
  fit <- estimate(pd[nrow(pd):1L, ], rd[nrow(rd):1L, ])
  expect_equal(fitted(fit), fitted(estimate(pd, rd)))
  new <- d[1L, c("Person", facets)]
  new$Person <- "unseen"
  new$Rater <- "new"
  np <- data.frame(Person = "unseen", Prior = 0.5)
  nr <- data.frame(Rater = "new", Tenure = 2)
  expect_true(is.finite(predict(fit, new, person_data = np, facet_data = nr)))
  expect_error(predict(fit, new, person_data = rbind(np, transform(np, Prior = 1)),
                       facet_data = nr), "person_data.*Duplicate ID")
  expect_error(predict(fit, new, person_data = np,
                       facet_data = rbind(nr, transform(nr, Tenure = 3))),
               "facet_data.*Duplicate ID")
})
