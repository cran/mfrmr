test_that("FACETS display maps join by facet and exact ID without changing results", {
  x <- read_facets_results(data.frame(Facet = c("Person", "Rater", "Person", "Person"),
    ElementID = c("01", "01", "02", "03"), Label = c("same", "rater", "same", "missing"),
    Measure = c(.2, -.4, 1, NA), Infit = c(.9, 1.1, 1, NA),
    Outfit = c(1.2, 1, 1, NA), Status = c(-1, 1, -2, 0)), common_scale = TRUE)
  saved <- x
  map <- data.frame(Facet = c("Person", "Rater", "Person"), ElementID = c("02", "01", "01"),
    DisplayLabel = c("\u53d7\u691c\u8005\u4e59", "\u8a55\u5b9a\u8005\u7532", "\u53d7\u691c\u8005\u7532"))
  for (type in c("measure", "wright", "fit", "infit_outfit")) {
    base <- plot(x, type = type, draw = FALSE)
    p <- plot(x, type = type, label_map = map, draw = FALSE)
    fields <- setdiff(names(base$data$plot_data), "DisplayLabel")
    expect_identical(p$data$plot_data[fields], base$data$plot_data[fields])
    expect_identical(p$data$data, base$data$data)
    expect_identical(p$data$sources, base$data$sources)
    expect_identical(p$data$plot_data$DisplayLabel, c(map$DisplayLabel[c(3, 2, 1)], "missing"))
    expect_identical(p$data$label_map, map)
  }
  filtered <- plot(x, facets = "Person", label_map = map, draw = FALSE)
  expect_identical(filtered$data$plot_data$DisplayLabel, c(map$DisplayLabel[c(3, 1)], "missing"))
  expect_identical(x, saved)
  expect_identical(plot(x, label_map = map[0, ], draw = FALSE)$data$plot_data,
    plot(x, draw = FALSE)$data$plot_data)
  file <- tempfile(fileext = ".rds")
  withr::defer(unlink(file))
  saveRDS(list(x = x, map = map, plotted = filtered), file)
  replay <- readRDS(file)
  expect_identical(plot(replay$x, facets = "Person", label_map = replay$map, draw = FALSE), filtered)
})

test_that("display maps reject ambiguous, unknown and missing exported IDs", {
  x <- read_facets_results(data.frame(Facet = "Person", ElementID = c("01", "02"),
    Label = c("A", "B"), Measure = c(0, 1)))
  map <- data.frame(Facet = "Person", ElementID = "01", DisplayLabel = "Local name")
  expect_error(plot(x, label_map = list(), draw = FALSE), "data frame")
  expect_error(plot(x, label_map = map[-3], draw = FALSE), "columns")
  duplicate_column <- map; duplicate_column[[4]] <- "extra"
  names(duplicate_column)[4] <- "DisplayLabel"
  expect_error(plot(x, label_map = duplicate_column, draw = FALSE), "columns")
  expect_error(plot(x, label_map = map[c(1, 1), ], draw = FALSE), "duplicate")
  for (column in names(map)) {
    bad <- map; bad[[column]] <- NA_character_
    expect_error(plot(x, label_map = bad, draw = FALSE), "nonempty character")
    bad[[column]] <- "  "
    expect_error(plot(x, label_map = bad, draw = FALSE), "nonempty character")
  }
  bad <- map; bad$ElementID <- 1L
  expect_error(plot(x, label_map = bad, draw = FALSE), "preserve IDs as text")
  bad$ElementID <- "1"
  expect_error(plot(x, label_map = bad, draw = FALSE), "match an exported")
  bad <- map; bad$Facet <- "Rater"
  expect_error(plot(x, label_map = bad, draw = FALSE), "match an exported")
  no_id <- read_facets_results(data.frame(Facet = "Person", Label = "01", Measure = 0))
  expect_error(plot(no_id, label_map = map, draw = FALSE), "match an exported")
  expect_error(plot(x, type = "ccc", label_map = map, draw = FALSE), "only to element")
})

test_that("mapped duplicate names remain distinct even when a name resembles an ID suffix", {
  x <- read_facets_results(data.frame(Facet = "Person", ElementID = c("01", "02", "03"),
    Label = c("A", "B", "C"), Measure = 1:3))
  map <- data.frame(Facet = "Person", ElementID = c("01", "02", "03"),
    DisplayLabel = c("same", "same", "same [01]"))
  p <- plot(x, label_map = map, draw = FALSE)$data$plot_data
  expect_false(anyDuplicated(p$DisplayLabel) > 0)
  expect_identical(p$ElementID, c("01", "02", "03"))
  expect_identical(p$Label, c("A", "B", "C"))
  expect_identical(p$DisplayLabel, c("same [01] [01]", "same [02]", "same [01] [03]"))
})
