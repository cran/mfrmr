test_that("linked MML rejects separation with one or both directions of cross-rating", {
  persons <- data.frame(Person = paste0("P", 1:24),
    Pool = factor(rep(c("A", "B"), each = 12)))
  cells <- expand.grid(Id = 1:24, Criterion = c("C1", "C2"),
    stringsAsFactors = FALSE)
  cells$Person <- persons$Person[cells$Id]
  cells$Rater <- as.character(persons$Pool[cells$Id])
  cells$Score <- 1L + (cells$Id + (cells$Criterion == "C2")) %% 4L
  one_way <- cells[cells$Id %in% 1:3, ]
  one_way$Rater <- "B"
  one_way$Score <- 1L
  reverse <- cells[cells$Id %in% 13:15, ]
  reverse$Rater <- "A"
  reverse$Score <- 4L

  fit <- function(link, model = "RSM", home = cells, person_data = persons,
                  formula = ~ 0 + Pool, anchors = NULL, group_anchors = NULL) {
    fit_mfrm(rbind(home, link), "Person", c("Rater", "Criterion"), "Score",
      rating_min = 1, rating_max = 4, method = "MML", model = model,
      step_facet = if (model == "PCM") "Criterion" else NULL,
      population_formula = formula, person_data = person_data, quad_points = 21L,
      anchors = anchors, group_anchors = group_anchors)
  }
  separated <- rbind(one_way, reverse)
  expect_error(fit(one_way), "No finite maximum-likelihood.*Rater")
  expect_error(fit(separated, "PCM"), "No finite maximum-likelihood.*Rater")
  balanced <- separated
  balanced$Score[c(1L, nrow(balanced))] <- c(4L, 1L)
  mixed <- separated
  mixed$Score[c(1L, nrow(mixed))] <- c(2L, 3L)
  result_mixed <- fit(mixed, "PCM")
  fair_mixed <- fair_average_table(result_mixed, facets = "Person")$raw_by_facet$Person
  expect_equal(nrow(fair_mixed), nrow(persons))
  expect_true(all(is.finite(fair_mixed$FairM)))

  # Additional covariates, signs and anchor combinations run in the full suite.
  skip_on_cran()
  persons$X <- rep(c(-1, 1, 1, -1), 6)
  expect_error(fit(one_way, formula = ~ 0 + Pool + X),
    "No finite maximum-likelihood.*Rater")
  flipped <- separated
  flipped$Score <- 5L - flipped$Score
  expect_error(fit(flipped), "No finite maximum-likelihood.*Rater")

  third <- cells[cells$Id %in% 1:12, ]
  third$Person <- paste0("C", third$Person)
  three_home <- rbind(cells, third)
  three_people <- rbind(persons, data.frame(Person = unique(third$Person),
    Pool = factor("C"), X = persons$X[1:12]))
  expect_error(fit(one_way, "PCM", three_home, three_people),
    "No finite maximum-likelihood.*Rater")

  three_people$Trend <- as.integer(three_people$Pool) - 2L
  separate_home <- three_home
  separate_home$Rater[separate_home$Person %in% third$Person] <- "C"
  bridge <- cells[cells$Id %in% 13:15, ]
  bridge$Rater <- "C"
  bridge$Score <- 1L
  links <- rbind(one_way, bridge)
  anchor <- function(level) data.frame(Facet = "Rater", Level = level, Anchor = 0)
  for (formula in list(~ Trend, ~ 0 + Trend)) {
    expect_error(fit(links, "PCM", separate_home, three_people, formula),
      "No finite maximum-likelihood.*Rater")
  }
  expect_error(fit(links, "PCM", separate_home, three_people, anchors = anchor("B")),
    "No finite maximum-likelihood.*Rater")
  group <- data.frame(Facet = "Rater", Level = c("A", "C"), Group = "AC", GroupValue = 0)
  expect_error(fit(links, "PCM", separate_home, three_people, group_anchors = group),
    "No finite maximum-likelihood.*Rater")

  result_bal <- fit(balanced, "PCM")
  fair_bal <- fair_average_table(result_bal, facets = "Person")$raw_by_facet$Person
  expect_equal(nrow(fair_bal), nrow(persons))
  expect_true(all(is.finite(fair_bal$FairM)))

  result <- fit(mixed, "PCM", three_home, three_people)
  fair <- fair_average_table(result, facets = "Person")$raw_by_facet$Person
  expect_equal(nrow(fair), nrow(three_people))
  expect_true(all(is.finite(fair$FairM)))

  three_people$Trend <- as.integer(three_people$Pool) - 2L
  separate_home <- three_home
  separate_home$Rater[separate_home$Person %in% third$Person] <- "C"
  bridge <- cells[cells$Id %in% 13:15, ]
  bridge$Rater <- "C"
  bridge$Score <- 1L
  links <- rbind(one_way, bridge)
  anchor <- function(level) data.frame(Facet = "Rater", Level = level, Anchor = 0)
  finite_fair <- function(result) {
    fair <- fair_average_table(result, facets = "Person")$raw_by_facet$Person
    expect_equal(nrow(fair), nrow(three_people))
    expect_true(all(is.finite(fair$FairM)))
  }
  finite_fair(fit(links, "PCM", separate_home, three_people, anchors = anchor(LETTERS[1:3])))
  links$Score[match(c("B", "C"), links$Rater)] <- 2L
  finite_fair(fit(links, "PCM", separate_home, three_people, ~ Trend))
  mixed_covariate <- fit(mixed, person_data = persons, formula = ~ 0 + Pool + X)
  expect_true(all(is.finite(fair_average_table(mixed_covariate, facets = "Person")$raw_by_facet$Person$FairM)))

  # A single numeric trend cannot freely shift A and C together relative to B.
  three_people$Pool <- as.integer(three_people$Pool)
  expect_s3_class(fit(one_way, "PCM", three_home, three_people, formula = ~ Pool), "mfrm_fit")
  three_people$Pool <- factor(three_people$Pool)
})
