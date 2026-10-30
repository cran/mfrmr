# Development-only aggregation; no numerical fitting or random generation.
path <- test_path("..", "..", "inst", "validation", "mfrm-wide-map-output-summary-20261001.R")
skip_if_not(file.exists(path), "Repository-only study aggregation helper")
accounting <- new.env(parent = globalenv())
sys.source(path, envir = accounting)

accounting_plan <- function(arms = "MML", output = "interval") {
  x <- expand.grid(Replicate = 1:50, Arm = arms, stringsAsFactors = FALSE)
  x$ConditionId <- "synthetic-accounting-fixture"
  x$Target <- "Rater1-Rater2"; x$Output <- output
  x$InputId <- paste0("input-", x$Replicate)
  x$Eligibility <- "eligible"; x$Reference <- "generating_parameter"; x$ReferenceValue <- 0
  x
}

accounting_results <- function(plan) {
  x <- plan[c("ConditionId", "Arm", "Replicate", "Target", "Output")]
  x$Status <- "returned"; x$Available <- FALSE; x$Estimate <- 0
  x$SE <- x$Lower <- x$Upper <- NA_real_
  x$SourceStage <- "q31:saved-stage"; x$Reason <- "Numerical interval review failed"
  x
}

test_that("interval availability and coverage retain all planned attempts", {
  p <- accounting_plan(); r <- accounting_results(p)
  r$Available[1:20] <- TRUE; r$SE[1:20] <- 1
  r$Lower[1:20] <- -1; r$Upper[1:20] <- 1; r$Lower[20] <- .1
  r$Status[41:50] <- "error"; r$Estimate[41:50] <- NA_real_; r$Reason[41:50] <- "Fit failed"
  s <- accounting$wide_output_summary(p, r)
  expect_true(s$Complete)
  expect_equal(s[c("Planned", "Eligible", "Available", "ReturnedUnavailable", "Errors")],
    data.frame(Planned=50, Eligible=50, Available=20, ReturnedUnavailable=20, Errors=10), ignore_attr=TRUE)
  expect_equal(s$ConditionalCoverage, 19/20)
  expect_equal(s$ReturnedAndCovered, 19/50)
  expect_equal(s$CoverageMCSE, sqrt(.95*.05/20))
  expect_equal(s$Availability, .4)
  none <- accounting_results(p)
  s0 <- accounting$wide_output_summary(p, none)
  expect_equal(s0$ReturnedAndCovered, 0)
  expect_equal(s0$Covered, 0)
  expect_true(is.na(s0$ConditionalCoverage))
  no_reference <- p; no_reference$Reference <- "none"; no_reference$ReferenceValue <- NA_real_
  sn <- accounting$wide_output_summary(no_reference, r)
  expect_true(is.na(sn$ConditionalCoverage) && is.na(sn$ReturnedAndCovered))
  expect_true(is.finite(sn$MeanWidth))
  expect_error(accounting$wide_output_summary(p, rbind(r,r[1,])), "Duplicate selected output")
  leak <- none; leak$Lower[1] <- -1
  expect_error(accounting$wide_output_summary(p, leak), "must not leak")
  untraceable <- none; untraceable$SourceStage[1] <- ""
  expect_error(accounting$wide_output_summary(p, untraceable), "need a source stage")
})

test_that("unstarted and interrupted work cannot masquerade as completed sampling", {
  p <- accounting_plan(); r <- accounting_results(p)
  r$Status[1:2] <- c("interrupted", "running")
  s <- accounting$wide_output_summary(p,r[-3,])
  expect_false(s$Complete)
  expect_equal(c(s$NotStarted,s$Running,s$Interrupted,s$Settled),c(1,1,1,47))
  expect_true(all(is.na(s[c("Availability","RMSE","ConditionalCoverage","ReturnedAndCovered")])))
  fixed <- p; fixed$Eligibility <- "fixed"
  sf <- accounting$wide_output_summary(fixed,r[FALSE,])
  expect_equal(c(sf$Planned,sf$Eligible,sf$Excluded,sf$Started),c(50,0,50,0))
  expect_true(sf$Complete)
  expect_error(accounting$wide_output_summary(fixed,r), "no sampling attempt")
})

test_that("paired errors keep both marginal delivery counts and input identity", {
  p <- accounting_plan(c("MML","JML"),"point"); r <- accounting_results(p)
  r$Available <- r$Arm=="MML" & r$Replicate<=40 | r$Arm=="JML" & r$Replicate %in% 16:45
  r$Estimate[r$Available] <- ifelse(r$Arm[r$Available]=="MML",1,2)
  r$Estimate[!r$Available] <- NA_real_;r$Status[!r$Available] <- "error"
  r$Reason[!r$Available] <- "Fit failed"
  s <- accounting$wide_paired_point_summary(p,r,"MML","JML")
  expect_equal(c(s$EligiblePairs,s$AAvailable,s$BAvailable,s$BothAvailable),c(50,40,30,25))
  expect_equal(s$MeanSquaredErrorDifference,-3)
  expect_equal(s$PairedMCSE,0)
  expect_match(s$Conditioning,"Both selected point outputs")
  mismatch <- p;mismatch$InputId[1] <- "different-input"
  expect_error(accounting$wide_paired_point_summary(mismatch,r,"MML","JML"),"identical inputs")
  root <- p;root$Reference[root$Arm=="JML"] <- "population_equation_root"
  expect_error(accounting$wide_paired_point_summary(root,r,"MML","JML"),"evaluation references")
  incomplete <- accounting$wide_paired_point_summary(p,r[-1,],"MML","JML")
  expect_false(incomplete$Complete)
  expect_true(is.na(incomplete$MeanSquaredErrorDifference))
})
