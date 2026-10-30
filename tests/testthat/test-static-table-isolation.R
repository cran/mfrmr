test_that("editing a guide cannot change subsequent guide or report output", {
  original <- mfrmr_output_guide("all")
  changed <- mfrmr_output_guide("all")
  changed$Question[1L] <- "caller annotation"
  changed$new_column <- TRUE
  expect_identical(mfrmr_output_guide("all"), original)
  plots <- mfrmr_output_guide("plots")
  modified <- plots
  modified[1L, 1L] <- "caller annotation"
  expect_identical(mfrmr_output_guide("plots"), plots)
})
