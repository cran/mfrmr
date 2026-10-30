comparison_review_data <- function(missing = FALSE) {
  d <- expand.grid(Person = paste0("P", 1:8), Rater = paste0("R", 1:4),
    Criterion = c("C1", "C2"), stringsAsFactors = FALSE)
  d$School <- ifelse(as.integer(sub("P", "", d$Person)) <= 4, "A", "B")
  d$Score <- (as.integer(factor(d$Person)) + as.integer(factor(d$Rater)) +
    as.integer(factor(d$Criterion))) %% 3
  if (missing) d <- d[!(d$Rater == "R4" & d$School == "B"), ]
  d
}
comparison_review_spec <- function(...) {
  c(list(facet = "Rater", group = "School", groups = c("A", "B"),
    contrasts = rbind(R1_panel = c(R1=.75, R2=-.25, R3=-.25, R4=-.25),
      R1_R2 = c(R1=1, R2=-1, R3=0, R4=0))), list(...))
}
comparison_review_run <- function(data, specification = comparison_review_spec()) {
  describe_mfrm_data(data, "Person", c("Rater", "Criterion", "School"), "Score",
    rating_min=0, rating_max=2, category_policy="preserve", include_agreement=FALSE,
    comparisons=specification)
}

test_that("prospective review separates full, partial and unidentified interaction targets", {
  local_mocked_bindings(fit_mfrm = function(...) stop("Review must not fit"), .package="mfrmr")
  for (model in c("RSM", "PCM")) {
    s <- comparison_review_spec(model=model)
    if (model == "PCM") s$step_facet <- "Criterion"
    full <- comparison_review_run(comparison_review_data(), s)$comparison_review
    expect_true(all(full$table$State == "all_parameters_identified"))
    d <- comparison_review_data(TRUE)
    review <- comparison_review_run(d, s)
    x <- review$comparison_review
    expect_identical(x$table$State, c("target_unidentified", "target_identified_only"))
    expect_true(all(is.na(x$table$Estimate)))
    expect_false(any(x$table$FormalInferenceEligible))
    expect_identical(summary(review)$comparison_review, x)
    expect_output(print(review), "Prospective interaction-comparison review")
    expect_output(print(summary(review)), "target_identified_only")
    expect_equal(x$cells$N[x$cells$Level == "R4" & x$cells$Group == "B"], 0L)
    expect_equal(x$table$N_Group2, c(24L, 16L))
    expect_equal(x$table$NPersons_Group2, c(4L, 4L))
    d <- rbind(d, subset(comparison_review_data(), Person == "P5" & Rater == "R4" & Criterion == "C1"))
    expect_true(all(comparison_review_run(d, s)$comparison_review$table$State == "all_parameters_identified"))
  }
})

test_that("coding preserves the target while origin direction reverses its map", {
  d <- comparison_review_data(TRUE); s <- comparison_review_spec()
  x <- comparison_review_run(d, s)$comparison_review
  s$groups <- rev(s$groups)
  y <- comparison_review_run(d[rev(seq_len(nrow(d))), ], s)$comparison_review
  expect_equal(y$target_map, -x$target_map)
  expect_identical(y$table$State, x$table$State)
  expect_identical(y$source$rows, x$source$rows)
  s$contrasts <- s$contrasts[, 4:1, drop=FALSE]
  expect_equal(comparison_review_run(d, s)$comparison_review$target_map, y$target_map)
  s$positive_facets <- "Rater"
  expect_identical(comparison_review_run(d, s)$comparison_review$table$State, x$table$State)
  expect_identical(x$fit_arguments$dummy_facets, "School")
  expect_identical(x$fit_arguments$facet_interactions, "Rater:School")
})

test_that("MML conditional ranks and unsupported GPCM stay outside global positive claims", {
  d <- comparison_review_data(TRUE)
  fixed <- comparison_review_run(d, comparison_review_spec(method="MML"))$comparison_review
  expect_true(all(fixed$table$State == "not_evaluated"))
  expect_true(all(fixed$table$ConditionalState == "all_parameters_identified"))
  expect_identical(fixed$table$FreePersonState, c("target_unidentified", "target_identified_only"))
  school <- comparison_review_run(d, comparison_review_spec(method="MML", population_formula=~School))$comparison_review
  expect_identical(school$table$State, c("target_unidentified", "not_evaluated"))
  expect_match(school$table$EvidenceScope[1], "survives_mixing")
  expect_equal(nrow(school$fit_arguments$person_data), 8L)
  testthat::with_mocked_bindings(
    expect_error(do.call(fit_mfrm,c(list(data=d),school$fit_arguments)), "population prepared"),
    mfrm_estimate=function(...) stop("population prepared"), .package="mfrmr")
  gpcm <- comparison_review_run(d, comparison_review_spec(model="GPCM"))$comparison_review
  expect_true(all(gpcm$table$State == "unsupported"))
  expect_true(all(is.na(gpcm$table$Rank)))
  expect_true(all(is.na(gpcm$table$StepFacet)))
})

test_that("a review rejects ambiguous targets and changed recording origins", {
  d <- comparison_review_data(); s <- comparison_review_spec()
  bad <- s; bad$contrasts <- bad$contrasts[, 1:3]
  expect_error(comparison_review_run(d, bad), "every facet level")
  bad <- s; bad$contrasts[1, 1] <- 1
  expect_error(comparison_review_run(d, bad), "summing to zero")
  bad <- s; bad$groups <- c("A", "A")
  expect_error(comparison_review_run(d, bad), "two distinct")
  bad <- s; bad$anchors <- data.frame()
  expect_error(comparison_review_run(d, bad), "named specification")
  d$School[1] <- "B"
  expect_error(comparison_review_run(d, s), "retain one origin")
  expect_error(comparison_review_run(comparison_review_data(), comparison_review_spec(model="PCM")), "explicit")
  expect_error(comparison_review_run(comparison_review_data(), comparison_review_spec(method="MML",population_formula=~Criterion)), "origin group")
  expect_null(comparison_review_run(comparison_review_data(), NULL)$comparison_review)
})

test_that("comparison bases distinguish population restrictions from response-model identification", {
  local_mocked_bindings(fit_mfrm=function(...)stop("Review must not fit"), .package="mfrmr")
  for (model in c("RSM", "PCM")) {
    d <- comparison_review_data(TRUE)
    s <- comparison_review_spec(model=model)
    if (model == "PCM") s$step_facet <- "Criterion"
    for (method in c("JML", "CML", "MML")) {
      s$method <- method
      x <- comparison_review_run(d,s)$comparison_review
      expect_identical(x$table$IdentificationBasis, if (method == "MML")
        c("population_restriction", "free_person_response_model") else
        c("unidentified", if (method == "CML") "conditional_support" else "free_person_response_model"))
      expect_true(all(is.na(x$table$Estimate)))
      expect_false(any(x$table$FormalInferenceEligible))
      if (method == "MML") {
        expect_true(all(x$table$State == "not_evaluated"))
        expect_match(x$table$IdentificationDetail[1], "global marginal identification is not established")
        expect_match(paste(x$notes,collapse=" "), "substantive restriction, not scale coding")
        expect_output(mfrmr_print_comparison_review(x), "Population-restriction rows")
      }
    }
    s$method <- "MML"; s$population_formula <- ~1
    common <- comparison_review_run(d,s)$comparison_review
    expect_identical(common$table$IdentificationBasis,c("population_restriction","free_person_response_model"))
    expect_match(paste(common$notes,collapse=" "),"estimates one common normal ability mean and SD")
    s$population_formula <- ~School
    regression <- comparison_review_run(d,s)$comparison_review
    expect_identical(regression$table$IdentificationBasis,c("unidentified","free_person_response_model"))
    expect_identical(regression$table$State,c("target_unidentified","not_evaluated"))
    expect_match(paste(regression$notes,collapse=" "),"common residual SD")
  }
  gpcm <- comparison_review_run(comparison_review_data(),comparison_review_spec(model="GPCM"))$comparison_review
  expect_true(all(gpcm$table$IdentificationBasis=="unsupported"))
  x <- comparison_review_run(comparison_review_data(TRUE),comparison_review_spec(method="MML"))$comparison_review
  # Historical reviews still print the interpretation from their retained map checks.
  x$table[c("IdentificationBasis","IdentificationDetail")] <- NULL
  expect_output(mfrmr_print_comparison_review(x),"population_restriction")
  x$table$FreePersonState <- "not_evaluated"
  expect_true(all(mfrmr_comparison_basis(x$table)$IdentificationBasis=="not_evaluated"))
  x$table$ConditionalState <- "not_evaluated"
  x$table$FreePersonState <- "all_parameters_identified"
  expect_true(all(mfrmr_comparison_basis(x$table)$IdentificationBasis=="not_evaluated"))
})

test_that("complete null spaces and numerical limits cannot become positive claims", {
  X <- cbind(diag(3), matrix(0,3,10))
  C <- rbind(observed=c(1,rep(0,12)),last_missing=c(rep(0,12),1))
  result <- mfrmr:::mfrmr_comparison_rowspace(X,C)
  expect_identical(result$state,c("target_identified_only","target_unidentified"))
  sensitive <- mfrmr:::mfrmr_comparison_rowspace(rbind(c(1,1),c(0,1e-9)),diag(2))
  expect_true(all(sensitive$state == "not_evaluated"))
  large <- Matrix::sparseMatrix(i=1,j=1,x=1,dims=c(1000001,2))
  expect_true(all(mfrmr:::mfrmr_comparison_rowspace(large,diag(2))$state == "not_evaluated"))
})

test_that("DFF tables, plots and reports retain a separately bound prospective review", {
  # Extended local/CI regression; see CONTRIBUTING.md.
  skip_on_cran()
  fit <- readRDS(test_path("fixtures", "ordinary-score-integration.rds"))$fits[["PCM-JML"]]
  d <- as.data.frame(fit$prep$data)
  d$School <- ifelse(as.integer(d$Person) %% 2, "A", "B")
  levels <- as.character(fit$prep$levels$Rater)
  W <- matrix(0, 1, length(levels), dimnames=list("first_minus_second",levels))
  W[1, 1:2] <- c(1, -1)
  review <- describe_mfrm_data(d, "Person", c(fit$prep$facet_names,"School"), "Score",
    rating_min=fit$prep$rating_min,rating_max=fit$prep$rating_max,category_policy="preserve",
    include_agreement=FALSE,comparisons=list(facet="Rater",group="School",groups=c("A","B"),contrasts=W))
  local_mocked_bindings(fit_mfrm=function(...) stop("No refit"), .package="mfrmr")
  plain <- analyze_dff(fit, facet="Rater",group="School",data=d,min_obs=1)
  x <- analyze_dff(fit, facet="Rater",group="School",data=d,min_obs=1,comparison_review=review)
  expect_identical(x$dif_table, plain$dif_table)
  expect_identical(x$comparison_review, review$comparison_review)
  expect_identical(summary(x)$comparison_review, x$comparison_review)
  for (type in c("summary", "heatmap")) {
    p <- plot(x,type=type,draw=FALSE)
    expect_identical(plot_data(p,"comparison_review"), x$comparison_review)
    expect_match(paste(plot_data(p,"notes")$Text, collapse=" "), "do not classify the displayed")
    if (requireNamespace("ggplot2", quietly=TRUE)) {
      expect_identical(plot_data(as_ggplot(p),"comparison_review"), x$comparison_review)
    }
  }
  report <- dif_report(x)
  expect_identical(report$comparison_review, x$comparison_review)
  expect_identical(summary(report)$comparison_review, x$comparison_review)
  expect_match(report$narrative, "not to the residual or subgroup-refit")
  restored <- unserialize(serialize(x, NULL))
  expect_identical(summary(restored)$comparison_review, x$comparison_review)
  expect_error(analyze_dff(fit,facet="Rater",group="School",data=d[-1,],min_obs=1,comparison_review=review), "same retained responses")
  expect_error(analyze_dff(fit,facet="Rater",group="School",data=d,min_obs=1,focal="B",comparison_review=review), "direction")
  shuffled <- analyze_dff(fit,facet="Rater",group="School",data=d[nrow(d):1,],min_obs=1,comparison_review=review)
  expect_identical(shuffled$comparison_review, x$comparison_review)
})
