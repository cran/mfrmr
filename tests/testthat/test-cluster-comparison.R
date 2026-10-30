test_that("group-count and weight comparisons use label-invariant pair counts", {
  skip_if_not_installed("cluster")
  input <- data.frame(ID = c("a", "b", "c", "d"),
                      X = c(0, 0, 1, 1), Y = c(0, 1, 0, 1))
  x <- mfrm_features(input, "ID", c("X", "Y"))
  horizontal <- mfrm_cluster(x, 2, weights = c(X = 3, Y = 1))
  vertical <- mfrm_cluster(x, 2, weights = c(X = 1, Y = 3))
  three <- mfrm_cluster(x, 3)
  out <- mfrm_cluster_compare(list(XWeighted = horizontal, YWeighted = vertical, Three = three))
  expect_equal(nrow(out$comparisons), 3L)
  first <- out$comparisons[1, ]
  expect_equal(first$Pairs, 6)
  expect_equal(first$SplitPairs, 2)
  expect_equal(first$JoinedPairs, 2)
  expect_equal(first$ChangedFraction, 2/3)
  expect_equal(first$AdjustedRand, -0.5)
  expect_equal(out$analysis_summary$K, c(2, 2, 3))
  expect_equal(out$analysis_summary$MinGroupSize, c(2, 2, 1))
  expect_equal(out$analysis_summary$MaxGroupSize, c(2, 2, 2))
  expect_true(all(is.na(out$comparisons$Imputation)))
  expect_equal(out$weights$Weight, c(3, 1, 1, 3, 1, 1))
  expect_identical(out$analyses$XWeighted, horizontal)
  expect_equal(summary(out), out$comparison_summary)
  expect_output(print(out), "not pooled inference")
  # Enumerate entity pairs independently, including the different-k comparison.
  ij <- utils::combn(1:4, 2)
  for (r in seq_len(nrow(out$comparisons))) {
    row <- out$comparisons[r, ]
    a <- out$analyses[[row$First]]$membership$Cluster
    b <- out$analyses[[row$Second]]$membership$Cluster
    same_a <- a[ij[1, ]] == a[ij[2, ]]
    same_b <- b[ij[1, ]] == b[ij[2, ]]
    expect_equal(row$SplitPairs, sum(same_a & !same_b))
    expect_equal(row$JoinedPairs, sum(!same_a & same_b))
    expect_equal(row$ChangedFraction, mean(same_a != same_b))
  }
  relabelled <- horizontal
  relabelled$membership$Cluster <- 10L - relabelled$membership$Cluster
  relabelled$membership <- relabelled$membership[c(4, 2, 1, 3), ]
  relabelled$feature_data <- mfrm_features(input[4:1, ], "ID", c("Y", "X"))
  same <- mfrm_cluster_compare(list(Original = horizontal, Reordered = relabelled))
  expect_equal(same$comparisons$ChangedFraction, 0)
  expect_equal(same$comparisons$AdjustedRand, 1)
  reverse <- mfrm_cluster_compare(list(Three = three, XWeighted = horizontal))
  forward <- subset(out$comparisons, First == "XWeighted" & Second == "Three")
  expect_equal(reverse$comparisons$SplitPairs, forward$JoinedPairs)
  expect_equal(reverse$comparisons$JoinedPairs, forward$SplitPairs)
})

test_that("comparison refuses different data and preserves excluded IDs", {
  skip_if_not_installed("cluster")
  input <- data.frame(ID = 1:6, X = c(0, 1, 2, 8, 9, NA))
  x <- mfrm_features(input, "ID", "X")
  two <- mfrm_cluster(x, 2, missing = "omit")
  three <- mfrm_cluster(x, 3, missing = "omit")
  out <- mfrm_cluster_compare(list(Two = two, Three = three))
  expect_equal(out$comparisons$Pairs, 10)
  expect_equal(out$analysis_summary$Excluded, c(1, 1))
  expect_identical(out$analyses$Three$membership$ID, as.character(1:6))
  expect_true(is.na(out$analyses$Three$membership$Cluster[6]))
  expect_error(mfrm_cluster_compare(list(two, three)), "nonblank names")
  expect_error(mfrm_cluster_compare(list(A = two)), "at least two")
  expect_error(mfrm_cluster_compare(list(A = two, A = three)), "unique")
  altered <- input
  altered$X[1] <- -1
  changed <- mfrm_cluster(mfrm_features(altered, "ID", "X"), 2, missing = "omit")
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "original feature values")
  altered$ID[1] <- 99
  changed <- mfrm_cluster(mfrm_features(altered, "ID", "X"), 2, missing = "omit")
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "same entity IDs")
  changed <- three
  changed$membership$Cluster[1] <- NA_integer_
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "same entities")
  changed <- three
  changed$membership$ID[1] <- "99"
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "Membership IDs")
  changed <- three
  changed$membership$Cluster[!is.na(changed$membership$Cluster)] <- 1L
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "valid group memberships")
})

test_that("imputed comparisons pair identical completions without refitting", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  original <- data.frame(ID = 1:6, X = c(0, 1, NA, 9, 10, NA))
  first <- second <- original
  first$X[3] <- 2
  second$X[3] <- 8
  long <- rbind(cbind(.imp = 0L, original), cbind(.imp = 1L, first), cbind(.imp = 2L, second))
  where <- is.na(original)
  where[6, "X"] <- FALSE
  model <- mice::as.mids(long, where = where)
  x <- mfrm_features(original, "ID", "X")
  cells <- subset(x$missing, ID == "3")
  two <- mfrm_cluster_imputed(x, model, cells, k = 2, missing = "omit")
  three <- mfrm_cluster_imputed(x, model, cells, k = 3, missing = "omit")
  # Existing results suffice; comparison must never invoke clustering again.
  local_mocked_bindings(mfrm_cluster = function(...) stop("Unexpected refit"), .package = "mfrmr")
  out <- mfrm_cluster_compare(list(Two = two, Three = three))
  expect_equal(out$comparisons$Imputation, 1:2)
  expect_equal(out$comparison_summary$Partitions, 2L)
  expect_equal(out$comparisons$Pairs, c(10, 10))
  expect_equal(out$comparison_summary$MeanChangedFraction, mean(out$comparisons$ChangedFraction))
  expect_equal(out$comparison_summary$MinChangedFraction, min(out$comparisons$ChangedFraction))
  expect_equal(out$comparison_summary$MaxChangedFraction, max(out$comparisons$ChangedFraction))
  expect_identical(out$analyses$Two, two)
  expect_equal(nrow(out$analysis_summary), 4L)
  for (i in 1:2) {
    separate <- mfrm_cluster_compare(list(Two = two$analyses[[i]], Three = three$analyses[[i]]))
    expect_equal(out$comparisons[i, c("SplitPairs", "JoinedPairs", "AdjustedRand")],
                 separate$comparisons[1, c("SplitPairs", "JoinedPairs", "AdjustedRand")],
                 ignore_attr = TRUE)
  }
  expect_error(mfrm_cluster_compare(list(A = two, B = two$analyses[[1]])), "only mfrm_cluster")
  changed <- three
  changed$analyses <- rev(changed$analyses)
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "Completed feature values.*1")
  changed <- three
  changed$analyses[[2]] <- NULL
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "same number of imputations")
  changed <- three
  changed$analyses[[2]]$feature_data$data$X[3] <- 7
  expect_error(mfrm_cluster_compare(list(A = two, B = changed)), "Completed feature values.*2")
})

test_that("feature selections compare the same entities and check every shared column", {
  skip_if_not_installed("cluster")
  input <- data.frame(ID = letters[1:4], X = c(0, 0, 1, 1), Y = c(0, 1, 0, 1))
  x <- mfrm_cluster(mfrm_features(input, "ID", "X"), 2)
  y <- mfrm_cluster_hierarchical(mfrm_features(input[4:1, ], "ID", "Y"), 2)
  xy <- mfrm_cluster(mfrm_features(input, "ID", c("X", "Y")), 2)
  out <- mfrm_cluster_compare(list(XOnly = x, YOnly = y, Both = xy))
  expect_equal(out$analysis_summary$Features, c(1L, 1L, 2L))
  expect_equal(out$weights$Feature, c("X", "Y", "X", "Y"))
  expect_equal(out$comparisons$ChangedFraction[1], 2/3)
  expect_equal(out$comparisons$AdjustedRand[1], -0.5)
  expect_equal(out$comparisons$Pairs, rep(6, 3))
  expect_identical(out$analyses$YOnly, y)
  reverse <- mfrm_cluster_compare(list(Both = xy, YOnly = y, XOnly = x))
  expect_equal(tail(reverse$comparisons$ChangedFraction, 1), 2/3)

  # XOnly shares no columns with either later analysis: still detect Y conflicts.
  changed <- input
  changed$Y[1] <- 2
  other <- mfrm_cluster(mfrm_features(changed, "ID", "Y"), 2)
  expect_error(mfrm_cluster_compare(list(X = x, Y = y, Changed = other)),
    "original feature values and types")
  changed$Y <- factor(input$Y)
  other <- mfrm_cluster(mfrm_features(changed, "ID", "Y"), 2)
  expect_error(mfrm_cluster_compare(list(X = x, Y = y, Changed = other)),
    "original feature values and types")

  changed <- input
  changed$X[4] <- NA_real_
  incomplete <- mfrm_cluster(mfrm_features(changed, "ID", c("X", "Y")), 2,
                             missing = "omit")
  expect_error(mfrm_cluster_compare(list(All = incomplete, YOnly = y)), "same entities")
  same_omission <- mfrm_cluster(mfrm_features(changed, "ID", "X"), 2, missing = "omit")
  omitted <- mfrm_cluster_compare(list(All = incomplete, XOnly = same_omission))
  expect_equal(omitted$analysis_summary$Excluded, c(1L, 1L))
  expect_equal(omitted$comparisons$Pairs, 3)
  expect_identical(omitted$analyses$XOnly$membership$ID, letters[1:4])
})

test_that("feature selection preserves imputation pairing even without shared features", {
  skip_if_not_installed("mice")
  skip_if_not_installed("cluster")
  original <- data.frame(ID = letters[1:6],
    X = c(0, 1, NA, 9, 10, 11), Y = c(0, 9, 1, NA, 2, 10),
    Complete = c(0, 0, 0, 1, 1, 1),
    Flag = c(FALSE, NA, FALSE, TRUE, TRUE, FALSE))
  first <- second <- original
  first$X[3] <- 2; second$X[3] <- 8
  first$Y[4] <- 8; second$Y[4] <- 3
  first$Flag[2] <- 0; second$Flag[2] <- 1
  # Numeric 0/1 storage for completed logical features matches mice's imputer.
  first$Flag <- as.numeric(first$Flag); second$Flag <- as.numeric(second$Flag)
  model <- mice::as.mids(rbind(cbind(.imp = 0L, original),
    cbind(.imp = 1L, first), cbind(.imp = 2L, second)))
  model$data$Flag <- original$Flag
  fit <- function(names) {
    review <- mfrm_features(original, "ID", names)
    mfrm_cluster_imputed(review, model, review$missing, k = 2)
  }
  x <- fit("X"); y <- fit("Y"); xy <- fit(c("X", "Y"))
  fixed <- fit("Complete"); flag <- fit("Flag")
  local_mocked_bindings(mfrm_cluster = function(...) stop("Unexpected refit"), .package = "mfrmr")
  out <- mfrm_cluster_compare(list(X = x, Y = y, Both = xy, Complete = fixed, Flag = flag))
  expect_equal(nrow(out$comparisons), 20L)
  expect_equal(out$analysis_summary$Features, rep(c(1L, 1L, 2L, 1L, 1L), each = 2))
  expect_true(all(out$comparison_summary$Partitions == 2L))
  expect_equal(nrow(fixed$imputed_cells), 0L)
  expect_identical(fixed$analyses[[1]]$membership, fixed$analyses[[2]]$membership)
  # Separate RDS reads recreate formula environments, not different imputations.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(y, path)
  restored <- readRDS(path)
  saved_environment <- environment(restored$imputation_model$formulas[[1]])
  reopened <- mfrm_cluster_compare(list(X = x, Y = restored))
  expect_equal(reopened$comparisons,
    mfrm_cluster_compare(list(X = x, Y = y))$comparisons)
  expect_identical(reopened$analyses$Y, restored)
  expect_identical(environment(restored$imputation_model$formulas[[1]]), saved_environment)
  changed_record <- restored
  changed_record$imputation_model$imp$Y[1, 1] <-
    changed_record$imputation_model$imp$Y[1, 1] + 1e-12
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed_record)), "same fitted mids")
  changed_record <- restored
  changed_record$imputation_model$formulas$Y <- Y ~ Complete
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed_record)), "same fitted mids")
  # Independent pair enumeration for the disjoint selections in each completion.
  ij <- utils::combn(1:6, 2)
  for (j in 1:2) {
    a <- x$analyses[[j]]$membership$Cluster
    b <- y$analyses[[j]]$membership$Cluster
    expect_equal(out$comparisons$ChangedFraction[j],
      mean((a[ij[1, ]] == a[ij[2, ]]) != (b[ij[1, ]] == b[ij[2, ]])))
  }
  changed <- y
  changed$analyses <- rev(changed$analyses)
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed)), "Completed feature values.*1")
  changed <- y
  changed$imputation_model$method["Y"] <- "norm"
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed)), "same fitted mids")
  changed <- y
  changed$feature_data$data$Y[1] <- -1
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed)), "original feature values")
  changed <- y
  changed$analyses[[1]]$feature_data <- x$analyses[[1]]$feature_data
  expect_error(mfrm_cluster_compare(list(X = x, Y = changed)), "selected features")
})
