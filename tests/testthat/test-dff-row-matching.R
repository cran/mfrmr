dff_matching_fixture <- function() {
  obs <- data.frame(Person = c("p1", "p1", "p1", "p2", "p2"),
    Item = c("i1", "i1", "i1", "i2", "i2"), Score = c(0, 0, 1, 0, 1), Weight = 1,
    Expected = c(.2, .2, .2, .8, .8), Var = .16)
  fit <- list(prep = list(facet_names = "Item", source_columns = list(person = "Person", score = "Score"),
    score_map = data.frame(OriginalScore = 0:1, InternalScore = 0:1)), config = list())
  d <- obs[c("Person", "Item", "Score")]; d$Group <- c("B", "A", "B", "A", "B")
  list(obs = obs, fit = fit, data = d)
}

test_that("response multiplicities and scores determine group allocation, not row order", {
  z <- dff_matching_fixture()
  local_mocked_bindings(compute_obs_table = function(res) z$obs)
  a <- .dff_match_observations(z$fit, z$data, "Group")
  b <- .dff_match_observations(z$fit, z$data[5:1, ], "Group")
  expect_identical(a$data, b$data)
  expect_equal(table(a$data$.group_var), table(z$data$Group))
  expect_equal(a$data$.group_var, c("A", "B", "B", "A", "B"))
  expect_true(all(a$cases$Included))
  expect_equal(sort(a$cases$FitRow), 1:5)
  expect_equal(z$data$Score, a$data$Score[match(a$cases$FitRow, sort(a$cases$FitRow))])
  subset <- .dff_match_observations(z$fit, z$data[-1, ], "Group")
  expect_equal(nrow(subset$data), 4L)
  expect_equal(sum(subset$data$.group_var == "B"), 2L)
  expect_equal(sum(subset$data$Score), 2)
})

test_that("excess copies are rejected without silently choosing between groups", {
  z <- dff_matching_fixture()
  local_mocked_bindings(compute_obs_table = function(res) z$obs)
  expect_error(.dff_match_observations(z$fit, rbind(z$data, z$data[1, ]), "Group"), "more copies")
  z$obs$Expected[2] <- .3
  expect_error(.dff_match_observations(z$fit, z$data, "Group"), "different model moments")
})

test_that("the case audit distinguishes missing groups, invalid responses and absent events", {
  z <- dff_matching_fixture()
  local_mocked_bindings(compute_obs_table = function(res) z$obs)
  extra <- z$data[rep(1, 5), ]; extra$Group <- c(NA, "A", "A", "B", "A")
  extra$Score[2] <- NA; extra$Person[3] <- "new"; extra$Score[4] <- 9; extra$Score[5] <- .5
  expect_message(a <- .dff_match_observations(z$fit, rbind(z$data, extra), "Group"), "Excluded 2 row")
  expect_identical(a$cases$Reason[6:10], c("missing_group", "invalid_or_missing_response",
    "not_in_fitted_responses", "not_in_fitted_responses", "invalid_or_missing_response"))
  expect_equal(sum(a$cases$Included), 5L)
  expect_true(all(is.na(a$cases$FitRow[6:10])))
  expect_error(suppressMessages(.dff_match_observations(z$fit, extra, "Group")), "No requested response")
})

test_that("source column names, score maps, trimming and fitted weights are replayed", {
  z <- dff_matching_fixture()
  z$fit$prep$source_columns <- list(person = "candidate id", score = "raw score", weight = "frequency")
  z$fit$prep$score_map$OriginalScore <- c(2, 4)
  z$obs$Weight <- c(1, 2, 1, 1, 1)
  local_mocked_bindings(compute_obs_table = function(res) z$obs)
  d <- z$data; names(d)[c(1, 3)] <- c("candidate id", "raw score")
  d[["raw score"]] <- factor(c(2, 2, 4, 2, 4)); d$frequency <- z$obs$Weight
  d[["candidate id"]][1] <- " p1 "
  a <- .dff_match_observations(z$fit, d, "Group")
  expect_true(all(a$cases$Included))
  expect_equal(a$data$Score, z$obs$Score)
  expect_equal(a$data$.group_var, z$data$Group)
  expect_error(.dff_match_observations(z$fit, d[setdiff(names(d), "frequency")], "Group"), "weight")
  d$frequency[2] <- 0
  a <- .dff_match_observations(z$fit, d, "Group")
  expect_equal(a$cases$Reason[2], "invalid_or_missing_response")
  normalized <- z$obs; normalized$Group <- z$data$Group
  expect_true(all(.dff_match_observations(z$fit, normalized, "Group", normalized = TRUE)$cases$Included))
})

test_that("declared missing-code policy does not turn omitted responses into fitted scores", {
  z <- dff_matching_fixture()
  z$fit$config$replay_inputs <- list(missing_codes = TRUE)
  local_mocked_bindings(compute_obs_table = function(res) z$obs)
  z$data$Score[1] <- 99
  a <- .dff_match_observations(z$fit, z$data, "Group")
  expect_identical(a$cases$Reason[1], "invalid_or_missing_response")
  expect_equal(nrow(a$data), 4L)
})

test_that("both public residual routes and summaries preserve exact selected observations", {
  fit <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[["PCM-JML"]]
  d <- as.data.frame(fit$prep$data)
  d$Group <- ifelse(as.integer(d$Person) %% 2, "A", "B")
  d <- d[seq(1, nrow(d), by = 2), ]
  a <- analyze_dff(fit, facet = "Rater", group = "Group", data = d, min_obs = 1)
  b <- dif_interaction_table(fit, facet = "Rater", group = "Group", data = d, min_obs = 1)
  obs <- compute_obs_table(fit)
  for (i in seq_len(nrow(a$cell_table))) {
    cell <- a$cell_table[i, ]; input <- which(as.character(d$Rater) == cell$Level & d$Group == cell$GroupValue)
    selected <- a$cases$FitRow[input]
    expect_equal(cell$N, length(input))
    expect_equal(cell$NPersons, length(unique(d$Person[input])))
    expect_equal(cell$ObsExpAvg, mean(obs$Observed[selected] - obs$Expected[selected]))
  }
  expect_identical(a$cases, b$cases)
  expect_identical(summary(a)$cases, a$cases)
  expect_identical(summary(b)$cases, b$cases)
  expect_equal(sum(a$cell_table$N), nrow(d))
  expect_identical(a$cell_table$NPersons, b$table$NPersons)
  expect_true(all(is.na(a$dif_table$p_value)))
  expect_equal(a$dif_table$NPersons_Group1, a$cell_table$NPersons[a$cell_table$GroupValue == "A"])
  old <- a; old$config$residual_output_version <- 2L
  expect_error(summary(old), "Recompute it")
  expect_error(dif_report(old), "Recompute it")
  expect_error(plot_dif_summary(old, draw = FALSE), "Recompute it")
})
