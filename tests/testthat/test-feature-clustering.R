test_that("feature review preserves entity identity, types, and missing reasons", {
  input <- data.frame(
    Rater = c("R|1", "R2", "R3"), Experience = c(NA, 2, 3),
    Specialty = c("Language", NA, "Science"),
    Training = ordered(c("low", "high", NA), levels = c("low", "medium", "high")),
    Certified = c(TRUE, FALSE, NA), Unselected = c(NA, NA, NA)
  )
  reasons <- data.frame(ID = c("R|1", "R2"),
    Feature = c("Experience", "Specialty"), Reason = c("Not collected", "Not applicable"))
  x <- mfrm_features(input, "Rater", c("Experience", "Specialty", "Training", "Certified"), reasons)
  expect_identical(x$data$Rater, input$Rater)
  expect_identical(x$data$Experience, input$Experience)
  expect_identical(levels(x$data$Training), levels(input$Training))
  expect_identical(x$feature_summary$Type, c("Numeric", "Nominal", "Ordinal", "Symmetric binary"))
  expect_equal(x$feature_summary$Missing, rep(1, 4))
  expect_equal(x$row_summary$MissingFeatures, c(1, 1, 2))
  expect_equal(nrow(x$missing), 4L)
  expect_identical(x$missing$Reason, c("Not collected", "Not applicable", "Not supplied", "Not supplied"))
  expect_equal(summary(x), x$feature_summary)
  expect_output(print(x), "no imputation is performed")
  expect_error(mfrm_features(input[c(1, 1), ], "Rater", "Experience"), "unique")
  expect_error(mfrm_features(input, "Rater", "Rater"), "excluding")
  expect_error(mfrm_features(input, "Rater", c("Experience", "Experience")), "distinct")
  reasons$Feature[1] <- "Training"
  expect_error(mfrm_features(input, "Rater", c("Experience", "Training"), reasons), "missing cells")
  reasons <- data.frame(ID = "R|1", Feature = "Experience", Reason = "Not collected")
  expect_error(mfrm_features(input, "Rater", "Experience", rbind(reasons, reasons)), "unique per")
  input$Experience[1] <- Inf
  expect_error(mfrm_features(input, "Rater", "Experience"), "infinite")
  input$Specialty[1] <- " "
  expect_error(mfrm_features(input, "Rater", "Specialty"), "blank")
  input$Date <- as.Date("2026-01-01") + 0:2
  expect_error(mfrm_features(input, "Rater", "Date"), "unsupported type")
  input$Complex <- 1:3 + 1i
  expect_error(mfrm_features(input, "Rater", "Complex"), "unsupported type")
})

test_that("mixed-feature groups retain ID mapping, separation, and profiles", {
  skip_if_not_installed("cluster")
  input <- data.frame(Rater = c("z", "a", "é", "R|4", "R5", "R6"),
    Experience = c(1, 2, 3, 12, 13, 14),
    Specialty = rep(c("Language", "Science"), each = 3),
    Certified = rep(c(FALSE, TRUE), each = 3))
  x <- mfrm_features(input, "Rater", c("Experience", "Specialty", "Certified"))
  groups <- mfrm_cluster(x, 2, weights = c(Certified = 1, Specialty = 2, Experience = 3))
  expect_identical(groups$membership$ID, input$Rater)
  expect_equal(outer(groups$membership$Cluster, groups$membership$Cluster, "=="),
               outer(rep(1:2, each = 3), rep(1:2, each = 3), "=="))
  # Hand-derived distances: within-group difference / 13 * 3/6;
  # across-group differences additionally have nominal/binary contributions 3/6.
  distances <- abs(outer(input$Experience, input$Experience, "-")) / 26 +
    outer(input$Specialty, input$Specialty, "!=") / 2
  expected_silhouette <- cluster::silhouette(groups$membership$Cluster, stats::as.dist(distances))[, "sil_width"]
  expect_equal(groups$membership$Silhouette, unname(expected_silhouette), tolerance = 1e-12)
  expect_setequal(groups$medoids$Rater, c("a", "R5"))
  expect_setequal(groups$membership$ID[groups$membership$Medoid], groups$medoids$Rater)
  expect_equal(groups$cluster_summary$N, c(3L, 3L))
  expect_equal(sort(groups$profiles$numeric$Mean), c(2, 13))
  expect_true(all(tapply(groups$profiles$categorical$Proportion,
    list(groups$profiles$categorical$Cluster, groups$profiles$categorical$Feature), sum) == 1))
  expect_identical(groups$settings$weights, c(Experience = 3, Specialty = 2, Certified = 1))
  expect_equal(groups$settings$numeric_ranges$Minimum, 1)
  expect_equal(groups$settings$numeric_ranges$Maximum, 14)
  expect_equal(summary(groups), groups$cluster_summary)
  expect_output(print(groups), "stability is not assessed")
  scaled <- mfrm_cluster(x, 2, weights = c(Experience = 3e307, Specialty = 2e307, Certified = 1e307))
  expect_equal(scaled$membership, groups$membership)
  expect_identical(x$data$Rater, input$Rater)
})

test_that("missing features require an explicit policy and excluded IDs remain visible", {
  skip_if_not_installed("cluster")
  input <- data.frame(Person = 1:6, A = c(NA, 2, 3, 12, 13, 14), B = c(NA, 1, 1, 2, 2, 2))
  x <- mfrm_features(input, "Person", c("A", "B"),
    data.frame(ID = 1, Feature = "A", Reason = "Not measured"))
  expect_error(mfrm_cluster(x, 2), "missing values")
  out <- mfrm_cluster(x, 2, missing = "omit")
  expect_identical(out$membership$ID, as.character(input$Person))
  expect_true(all(is.na(out$membership[1, c("Cluster", "Medoid", "Silhouette")])))
  expect_equal(out$settings$included, 5)
  expect_equal(out$settings$excluded, 1)
  expect_identical(out$feature_data$missing, x$missing)
  expect_equal(out$settings$numeric_ranges$Minimum, c(2, 1))
  expect_equal(sum(out$cluster_summary$N), 5)
  expect_equal(names(out$profiles$categorical), c("Cluster", "Feature", "Level", "N", "Proportion"))
  expect_equal(nrow(out$profiles$categorical), 0L)
  # A stale cached review must not override an edited missing cell.
  x$data$A[2] <- NA_real_
  expect_error(mfrm_cluster(x, 2), "missing values")
})

test_that("unusable features and clustering requests are refused before fitting", {
  skip_if_not_installed("cluster")
  x <- mfrm_features(data.frame(ID = letters[1:4], A = c(0, 0, 1, 1)), "ID", "A")
  expect_error(mfrm_cluster(x, 1), "integer from 2")
  expect_error(mfrm_cluster(x, 2.5), "integer from 2")
  expect_error(mfrm_cluster(x, 2 + 1i), "integer from 2")
  expect_error(mfrm_cluster(x, 4), "integer from 2")
  expect_error(mfrm_cluster(x, 3), "distinct included")
  expect_error(mfrm_cluster(x, 2, weights = 1), "named")
  expect_error(mfrm_cluster(x, 2, weights = c(A = 0)), "positive finite")
  expect_error(mfrm_cluster(x, 2, weights = c(B = 1)), "covering")
  expect_error(mfrm_cluster(x, 2, weights = c(A = 1 + 1i)), "positive finite")
  constant <- mfrm_features(data.frame(ID = 1:4, A = rep(1, 4)), "ID", "A")
  expect_error(mfrm_cluster(constant, 2), "must vary")
  empty <- mfrm_features(data.frame(ID = 1:4, A = rep(NA_real_, 4)), "ID", "A")
  expect_equal(empty$feature_summary$Distinct, 0)
  expect_error(mfrm_cluster(empty, 2, missing = "omit"), "included entities")
  large <- mfrm_features(data.frame(ID = 1:5001, A = 1:5001), "ID", "A")
  expect_error(mfrm_cluster(large, 2), "5,000")
  nominal <- mfrm_features(data.frame(ID = letters[1:4], A = c("a", "a", "b", "b")), "ID", "A")
  out <- mfrm_cluster(nominal, 2)
  expect_equal(nrow(out$profiles$numeric), 0)
  expect_equal(names(out$profiles$numeric), c("Cluster", "Feature", "N", "Mean", "Median", "SD"))
})

test_that("numeric overflow and weight underflow cannot silently discard features", {
  skip_if_not_installed("cluster")
  input <- data.frame(ID = letters[1:6],
    Wide = c(-1e308, -1e307, 0, 1e307, 5e307, 1e308), Group = rep(c("a", "b"), each = 3))
  x <- mfrm_features(input, "ID", c("Wide", "Group"))
  expect_error(mfrm_cluster(x, 2), "ranges overflow.*rescale")
  input$Wide <- input$Wide / 1e308
  x <- mfrm_features(input, "ID", c("Wide", "Group"))
  groups <- mfrm_cluster(x, 2)
  expect_true(all(is.finite(groups$membership$Silhouette)))
  expect_true(all(groups$membership$Silhouette < 1))
  expect_error(mfrm_cluster(x, 2, weights = c(Wide = 1e-308, Group = 1e308)), "underflow")
})

test_that("original-unit profile SDs survive squared overflow and underflow", {
  skip_if_not_installed("cluster")
  values <- c(1, 2, 3, 11, 12, 13)
  for (factor in c(1e200, 1e-200, 1e-160)) {
    x <- mfrm_features(data.frame(ID = letters[1:6], X = values * factor), "ID", "X")
    for (fun in list(mfrm_cluster_pam, mfrm_cluster_hierarchical)) {
      out <- fun(x, 2)
      expect_equal(out$profiles$numeric$SD / factor, c(1, 1), tolerance = 1e-12)
      expect_equal(out$profiles$numeric$Mean / factor, c(2, 12), tolerance = 1e-12)
      expect_equal(out$membership$Cluster, c(1L, 1L, 1L, 2L, 2L, 2L))
    }
  }
  # Preserve a true zero and an undefined singleton SD, independently of scale.
  profiles <- cluster_feature_profiles(data.frame(X = c(1, 1, 3) * 1e200), c(1, 1, 2))
  expect_identical(profiles$numeric$SD, c(0, NA_real_))
  # The shared helper also serves k-means; ordinary results remain exact.
  x <- mfrm_features(data.frame(ID = letters[1:6], X = values), "ID", "X")
  out <- mfrm_cluster_kmeans(x, 2, silhouette = FALSE)
  expected <- vapply(seq_len(2), function(g) sd(values[out$membership$Cluster == g]), numeric(1))
  expect_identical(out$profiles$numeric$SD, expected)
})

test_that("Gower ordinal gaps and singleton silhouettes retain their definitions", {
  skip_if_not_installed("cluster")
  ordinal <- ordered(c("b", "b", "d", "e"), levels = letters[1:6])
  x <- mfrm_features(data.frame(ID = letters[1:4], Grade = ordinal), "ID", "Grade")
  out <- mfrm_cluster_hierarchical(x, 3)
  # Observed codes 2, 2, 4, 5 use range 3; unused middle level c is not dropped.
  distance <- abs(outer(c(2, 2, 4, 5), c(2, 2, 4, 5), "-")) / 3
  expected <- hclust(as.dist(distance), method = "average")
  expect_equal(out$tree$height, expected$height, tolerance = 1e-14)
  expect_identical(out$tree$merge, expected$merge)
  expect_identical(out$membership$Silhouette, c(1, 1, 0, 0))
  expect_identical(levels(out$feature_data$data$Grade), letters[1:6])
})

test_that("mice completions preserve observed mixed features and model diagnostics", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  attributes <- mice::nhanes2
  attributes$Person <- paste0("P", seq_len(nrow(attributes)))
  attributes$Certified <- seq_len(nrow(attributes)) %% 2 == 0
  attributes$Training <- ordered(rep(c("low", "medium", "high"), length.out = nrow(attributes)),
                                 levels = c("low", "medium", "high"))
  x <- mfrm_features(attributes, "Person", c("bmi", "chl", "age", "hyp", "Certified", "Training"))
  method <- mice::make.method(attributes)
  method["Person"] <- ""
  predictors <- mice::make.predictorMatrix(attributes)
  predictors[, "Person"] <- 0
  predictors["Person", ] <- 0
  model <- mice::mice(attributes, m = 3, maxit = 2, method = method,
    predictorMatrix = predictors, seed = 42, printFlag = FALSE)
  out <- mfrm_cluster_imputed(x, model, x$missing, k = 2)
  expect_identical(out$imputation_model, model)
  expect_identical(out$feature_data, x)
  expect_equal(out$imputed_cells, x$missing)
  expect_equal(length(out$analyses), 3L)
  expect_true(all(out$analysis_summary$Included == nrow(attributes)))
  expect_true(all(out$analysis_summary$Excluded == 0L))
  expect_equal(out$co_membership, t(out$co_membership))
  expect_equal(unname(diag(out$co_membership)), rep(1, nrow(attributes)))
  expect_equal(summary(out), out$analysis_summary)
  expect_output(print(out), "not membership probabilities or sampling stability")
  for (a in out$analyses) {
    expect_identical(a$membership$ID, attributes$Person)
    observed <- !is.na(attributes$bmi)
    expect_identical(a$feature_data$data$bmi[observed], attributes$bmi[observed])
    expect_identical(a$feature_data$data$Training, attributes$Training)
    expect_identical(a$feature_data$data$Certified, attributes$Certified)
    expect_identical(levels(a$feature_data$data$hyp), levels(attributes$hyp))
  }
  # Auxiliary variables may belong to the imputation model without entering distances.
  narrow <- mfrm_features(attributes, "Person", c("bmi", "chl"))
  fewer <- mfrm_cluster_imputed(narrow, model, narrow$missing, k = 2)
  expect_identical(names(fewer$settings$weights), c("bmi", "chl"))
  expect_identical(fewer$imputation_model$data$age, attributes$age)
  changed_levels <- model
  changed_levels$data$age <- factor(as.character(attributes$age), levels = rev(levels(attributes$age)))
  expect_error(mfrm_cluster_imputed(x, changed_levels, x$missing, 2), "factor levels/order")
  changed_type <- model
  changed_type$data$Certified <- as.integer(attributes$Certified)
  expect_error(mfrm_cluster_imputed(x, changed_type, x$missing, 2), "Feature types")
})

test_that("mice imputed logical and ordered features retain their meaning", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  attributes <- mice::nhanes2
  attributes$hyp <- as.logical(as.integer(attributes$hyp) - 1L)
  attributes$age <- ordered(attributes$age)
  attributes$age[c(2, 5)] <- NA
  attributes$Person <- paste0("P", seq_len(nrow(attributes)))
  method <- mice::make.method(attributes)
  method["Person"] <- ""
  predictors <- mice::make.predictorMatrix(attributes)
  predictors[, "Person"] <- 0
  predictors["Person", ] <- 0
  model <- mice::mice(attributes, m = 2, maxit = 2, method = method,
    predictorMatrix = predictors, seed = 5, printFlag = FALSE)
  x <- mfrm_features(attributes, "Person", c("age", "bmi", "hyp", "chl"))
  out <- mfrm_cluster_imputed(x, model, x$missing, 2)
  for (a in out$analyses) {
    expect_type(a$feature_data$data$hyp, "logical")
    expect_true(is.ordered(a$feature_data$data$age))
    expect_identical(levels(a$feature_data$data$age), levels(attributes$age))
  }
  bad <- model
  bad$imp$hyp[1, 2] <- 2
  expect_error(mfrm_cluster_imputed(x, bad, x$missing, 2), "Imputation 2: Feature types")
})

test_that("imputation co-membership uses all completions and retains undefined cells", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  original <- data.frame(ID = 1:6, A = c(0, 1, NA, 9, 10, NA))
  first <- second <- original
  first$A[3] <- 2
  second$A[3] <- 8
  long <- rbind(cbind(.imp = 0L, original), cbind(.imp = 1L, first), cbind(.imp = 2L, second))
  where <- is.na(original)
  where[6, "A"] <- FALSE
  model <- mice::as.mids(long, where = where)
  reasons <- data.frame(ID = c(3, 6), Feature = "A", Reason = c("Not recorded", "Not applicable"))
  # The reviewed table is deliberately reordered; rows must be matched by ID.
  x <- mfrm_features(original[c(6, 3, 1, 5, 2, 4), ], "ID", "A", reasons)
  cells <- subset(x$missing, Reason == "Not recorded")
  expect_error(mfrm_cluster_imputed(x, model, cells, 2), "not selected for imputation")
  out <- mfrm_cluster_imputed(x, model, cells, 2, missing = "omit")
  # The unknown entity joins {1,2} in the first completion and {4,5} in the second.
  expected <- matrix(c(1, 1, .5, 0, 0,
                       1, 1, .5, 0, 0,
                       .5, .5, 1, .5, .5,
                       0, 0, .5, 1, 1,
                       0, 0, .5, 1, 1), 5, 5)
  expect_equal(unname(out$co_membership[as.character(1:5), as.character(1:5)]), expected)
  expect_true(all(is.na(out$co_membership["6", ])))
  expect_true(all(is.na(out$co_membership[, "6"])))
  expect_identical(rownames(out$co_membership), x$row_summary$ID)
  expect_true(all(out$analysis_summary$Included == 5L))
  expect_true(all(out$analysis_summary$Excluded == 1L))
  expect_equal(out$settings$imputations, 2L)
  expect_identical(out$imputed_cells$Reason, "Not recorded")
  expect_identical(out$feature_data$missing$Reason, c("Not applicable", "Not recorded"))
  for (a in out$analyses) {
    expect_true(is.na(a$membership$Cluster[1]))
    expect_identical(a$feature_data$missing$Reason, "Not applicable")
  }
  # Relabelling a partition cannot change pairwise co-membership.
  relabelled <- lapply(out$analyses, function(a) 10L - a$membership$Cluster)
  independent <- (outer(relabelled[[1]], relabelled[[1]], "==") +
                    outer(relabelled[[2]], relabelled[[2]], "==")) / 2
  expect_equal(unname(out$co_membership), independent)

  bad <- model
  bad$imp$A[1, 2] <- NA_real_
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2, missing = "omit"), "Imputation 2: Missing cells")
  bad <- model
  bad$imp$A[1, 2] <- Inf
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2, missing = "omit"), "Imputation 2:.*infinite")
  bad <- model
  bad$data$A[1] <- .01
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2, missing = "omit"), "Observed feature values")
  bad <- model
  bad$data$ID[1] <- 100
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2, missing = "omit"), "Entity IDs")
  bad <- model
  bad$where[1, "A"] <- TRUE
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2, missing = "omit"), "must match.*impute")
  expect_error(mfrm_cluster_imputed(x, model, x$missing, 2), "must match.*impute")
  expect_error(mfrm_cluster_imputed(x, model, rbind(cells, cells), 2), "unique missing cells")
  expect_error(mfrm_cluster_imputed(x, model, data.frame(ID = 1, Feature = "A"), 2), "unique missing cells")
  expect_error(mfrm_cluster_imputed(x, model, cells[FALSE, ], 2), "explicitly list")
  bad <- model
  bad$m <- 1L
  expect_error(mfrm_cluster_imputed(x, bad, cells, 2), "at least two")
})

test_that("one unusable completed partition stops the entire comparison", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  original <- data.frame(ID = 1:4, A = c(0, 0, 1, NA))
  first <- second <- original
  first$A[4] <- 2
  second$A[4] <- 1
  long <- rbind(cbind(.imp = 0L, original), cbind(.imp = 1L, first), cbind(.imp = 2L, second))
  model <- mice::as.mids(long)
  x <- mfrm_features(original, "ID", "A")
  expect_error(mfrm_cluster_imputed(x, model, x$missing, 3), "Imputation 2:.*distinct included")
  large <- mfrm_features(data.frame(ID = 1:5001, A = c(NA, 2:5001)), "ID", "A")
  expect_error(mfrm_cluster_imputed(large, model, large$missing, 2), "5,000 total")
})
