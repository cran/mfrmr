test_that("anchor_review accessor uses only the canonical field", {
  new_review <- structure(list(source = "new"), class = "mfrm_anchor_review")

  fit_like <- list(config = list(anchor_review = new_review))
  expect_identical(anchor_review(fit_like)$source, "new")
  expect_identical(anchor_review(new_review), new_review)

  expect_null(anchor_review(list(config = list()), required = FALSE))
  expect_error(
    anchor_review(list(config = list(anchor_audit = list(source = "old")))),
    "No anchor-review component"
  )
})

test_that("precision_review accessor uses only the canonical table", {
  new_review <- data.frame(Check = "new", Status = "ok", Detail = "canonical", stringsAsFactors = FALSE)

  diagnostics <- list(precision_review = new_review)
  expect_identical(precision_review(diagnostics), new_review)

  expect_null(precision_review(list(), required = FALSE))
  expect_error(
    precision_review(list(precision_audit = data.frame())),
    "No precision-review component"
  )
})
