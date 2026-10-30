test_that("result CSVs preserve UTF-8 names, quoting and missing values in the C locale", {
  withr::local_locale(c(LC_CTYPE="C"))
  x <- data.frame(label=c("\u89b3\u70b9, \"A\"", "line1\nline2", NA_character_),
    value=c(1/3, NA_real_, -1.25e100), group=factor(c("\u5224\u5b9a", "x", NA)),
    date=as.Date(c("2026-10-23", NA, "2026-01-01")), check.names=FALSE)
  names(x)[1] <- "\u89b3\u70b9 \"name\""
  path <- tempfile(fileext=".csv"); withr::defer(unlink(path))
  expect_warning(mfrm_results_export_write_csv(x,path), NA)
  y <- read.csv(path, encoding="UTF-8", check.names=FALSE, na.strings="")
  expect_identical(names(y),names(x))
  expect_identical(y[[1]],x[[1]])
  expect_identical(y$group,as.character(x$group))
  expect_identical(y$date,as.character(x$date))
  expect_equal(y$value,x$value,tolerance=1e-14)
  expect_false(any(grepl("<U+",readLines(path,encoding="UTF-8"),fixed=TRUE)))
  mfrm_results_export_write_csv(x[FALSE,],path)
  expect_identical(names(read.csv(path,encoding="UTF-8",check.names=FALSE)),names(x))
})
