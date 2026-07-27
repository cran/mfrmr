test_that("as_kable appends a pipe-table note exactly once", {
  skip_if_not_installed("knitr")

  tbl <- structure(
    list(
      table = data.frame(Term = c("A", "B"), Estimate = c(-0.1, 0.2)),
      caption = "Toy caption",
      note = "Toy note"
    ),
    class = "apa_table"
  )

  rendered <- as_kable(tbl, format = "pipe")
  rendered_text <- paste(as.character(rendered), collapse = "\n")
  note_locations <- gregexpr("Note. Toy note", rendered_text, fixed = TRUE)[[1]]

  expect_s3_class(rendered, "knitr_kable")
  expect_identical(attr(rendered, "format"), "pipe")
  expect_identical(sum(note_locations > 0L), 1L)
  expect_match(rendered_text, "Table: Toy caption", fixed = TRUE)
  expect_true(endsWith(rendered_text, "Note. Toy note"))
})

test_that("as_kable leaves note-free pipe tables as native kable vectors", {
  skip_if_not_installed("knitr")

  tbl <- structure(
    list(table = data.frame(A = 1:2), caption = "", note = ""),
    class = "apa_table"
  )
  rendered <- as_kable(tbl, format = "pipe")

  expect_s3_class(rendered, "knitr_kable")
  expect_gt(length(rendered), 1L)
  expect_false(any(grepl("Note.", as.character(rendered), fixed = TRUE)))
})
