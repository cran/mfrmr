# Graph contracts use synthetic links; no estimation or inference is implied.
graph_chain <- function(waves = c("A", "B", "C")) {
  structure(list(
    config = list(waves = waves),
    links = data.frame(Link = 1:2, FromID = 1:2, ToID = 2:3,
      From = waves[1:2], To = waves[2:3], N_Common = c(2L, 1L),
      N_Retained = c(1L, 0L), LinkSupportAdequate = FALSE),
    cumulative = data.frame(Wave = waves, Cumulative_Offset = c(0, 0.1, NA)),
    element_detail = data.frame(LinkID = c(1L, 1L, 2L),
      Facet = "Item", Level = c("I1", "I2", "I1"),
      Retained = c(TRUE, FALSE, FALSE), Flag = c(FALSE, TRUE, TRUE))
  ), class = c("mfrm_equating_chain", "list"))
}

test_that("graph preserves isolated waves, pair screening and shared identities", {
  x <- graph_chain()
  p <- plot(x, type = "graph", draw = FALSE)
  d <- p$data$data
  expect_equal(d$n_waves, 3L)
  expect_equal(d$n_anchors, 2L)
  expect_equal(d$n_element_nodes, 3L)
  expect_equal(d$nodes$Label[d$nodes$Kind == "wave"], c("A", "B", "C"))
  expect_equal(d$nodes$RetainedDegree[d$nodes$Kind == "wave"], c(1L, 1L, 0L))
  expect_equal(sum(d$edges$Retained), 2L)
  expect_equal(sum(d$edges$Status == "excluded"), 4L)
  expect_equal(length(unique(d$elements$AnchorId[d$elements$Level == "I1"])), 1L)
  expect_equal(length(unique(d$elements$NodeId[d$elements$Level == "I1"])), 2L)
  expect_true(all(c("LinkID", "Retained", "Flag", "Reason") %in% names(d$edges)))
  expect_true(any(grepl("No retained.*C", p$data$notes$Text)))
  expect_true(any(grepl("Thin or absent", p$data$notes$Text)))
  expect_output(print(p), "Connectedness and retention")
})

test_that("zero common elements and zero links still produce wave vertices", {
  x <- graph_chain()
  x$element_detail <- data.frame()
  x$links$N_Common <- x$links$N_Retained <- 0L
  for (empty_links in c(FALSE, TRUE)) {
    if (empty_links) x$links <- data.frame()
    p <- plot(x, type = "graph", draw = FALSE)
    expect_equal(p$data$data$n_waves, 3L)
    expect_equal(p$data$data$n_anchors, 0L)
    expect_equal(nrow(p$data$data$edges), 0L)
    expect_equal(nrow(p$data$data$nodes), 3L)
    grDevices::pdf(NULL, width = 5, height = 4)
    tryCatch(expect_no_warning(plot(x, type = "graph")), finally = grDevices::dev.off())
  }
})

test_that("unknown retention and retained-but-flagged are not called unflagged retention", {
  x <- graph_chain()
  x$element_detail$Retained <- c(TRUE, TRUE, NA)
  x$element_detail$Flag <- c(TRUE, NA, FALSE)
  p <- plot(x, type = "graph", draw = FALSE)
  expect_equal(p$data$data$elements$Status, c("retained_review", "retained_review", "candidate"))
  x$element_detail$Retained <- x$element_detail$Flag <- NULL
  p <- plot(x, type = "graph", draw = FALSE)
  expect_true(all(p$data$data$edges$Status == "candidate"))
  expect_true(all(p$data$data$nodes$RetainedDegree == 0L))
})

test_that("wave and element identity survives duplicate names and delimiter punctuation", {
  x <- graph_chain(c("Item: I1 -> A", "Item: I1 -> A", "B | C"))
  x$element_detail$Facet <- c("a: b", "a", "a: b")
  x$element_detail$Level <- c("c", "b: c", "c")
  p <- plot(x, type = "graph", draw = FALSE)
  expect_equal(nrow(p$data$data$nodes), 6L)
  expect_equal(anyDuplicated(p$data$data$nodes$NodeId), 0L)
  expect_equal(p$data$data$n_anchors, 2L)
  expect_equal(p$data$data$links$from, c("wave_1", "wave_2"))
  expect_equal(p$data$data$links$to, c("wave_2", "wave_3"))
  # Old native chains are matched against whole labels, without splitting arrows.
  x <- graph_chain(c("A -> B", "C | D", "E"))
  x$element_detail$Link <- c("A -> B -> C | D", "A -> B -> C | D", "C | D -> E")
  x$element_detail$LinkID <- x$links$FromID <- x$links$ToID <- NULL
  expect_equal(plot(x, type = "graph", draw = FALSE)$data$data$n_waves, 3L)
  x$element_detail$Link[1] <- "unknown"
  expect_error(plot(x, type = "graph", draw = FALSE), "do not match")
})

test_that("malformed graph references and flags fail explicitly", {
  x <- graph_chain()
  x$links$ToID[2] <- 4L
  expect_error(plot(x, type = "graph", draw = FALSE), "Invalid wave IDs")
  x <- graph_chain()
  x$element_detail$LinkID[1] <- 10L
  expect_error(plot(x, type = "graph", draw = FALSE), "do not match")
  x <- graph_chain()
  x$element_detail$Retained <- c("TRUE", "FALSE", "FALSE")
  expect_error(plot(x, type = "graph", draw = FALSE), "must be logical")
  expect_error(plot(graph_chain(), type = "graph", show_notes = NA), "TRUE or FALSE")
})

test_that("monochrome and clean drawing preserve data and the caller's panel sequence", {
  x <- graph_chain()
  expected <- plot(x, type = "graph", draw = FALSE, preset = "monochrome",
    show_title = FALSE, show_notes = FALSE)
  expect_equal(length(unique(expected$data$styles$Linetype)), 4L)
  rgb <- grDevices::col2rgb(expected$data$styles$Color)
  expect_equal(rgb[1, ], rgb[2, ])
  expect_equal(rgb[2, ], rgb[3, ])
  expect_equal(mfrmr:::.mfrm_plot_display_title(expected$data), "")
  expect_true(nrow(expected$data$notes) >= 3L)
  grDevices::pdf(NULL, width = 10, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  old <- graphics::par(c("mar", "bg", "xpd"))
  expect_no_warning(actual <- plot(x, type = "graph", preset = "monochrome",
    show_title = FALSE, show_notes = FALSE))
  expect_identical(actual, expected)
  expect_equal(graphics::par(c("mar", "bg", "xpd")), old)
  expect_equal(graphics::par("mfg")[1:2], c(1L, 1L))
  graphics::plot(1:3)
  expect_equal(graphics::par("mfg")[1:2], c(1L, 2L))
})

test_that("supported paths do not bypass excluded pair-specific links", {
  skip_if_not_installed("igraph")
  x <- graph_chain()
  x$element_detail$Retained <- c(TRUE, FALSE, TRUE)
  d <- plot(x, type = "graph", draw = FALSE)$data$data
  g <- igraph::graph_from_data_frame(d$edges[d$edges$Retained, 1:2],
    directed = FALSE, vertices = d$nodes$NodeId)
  expect_true(is.finite(igraph::distances(g, "wave_1", "wave_3")[1, 1]))
  # The same element is excluded only in B-C; A-C must now be disconnected.
  x$element_detail$Retained[3] <- FALSE
  d <- plot(x, type = "graph", draw = FALSE)$data$data
  g <- igraph::graph_from_data_frame(d$edges[d$edges$Retained, 1:2],
    directed = FALSE, vertices = d$nodes$NodeId)
  expect_equal(igraph::distances(g, "wave_1", "wave_3")[1, 1], Inf)
})

test_that("removal deletes an element across all links and counts only new losses", {
  x <- graph_chain()
  x$element_detail$Retained <- c(TRUE, TRUE, TRUE)
  x$element_detail$Flag <- FALSE
  x$config$waves <- c(x$config$waves, "D") # already isolated at baseline
  p <- plot(x, type = "anchor_removal", draw = FALSE)
  d <- p$data$data
  expect_equal(d$summary$Components, 2L)
  expect_equal(d$summary$ConnectedPairs, 3L)
  i1 <- d$removal[d$removal$Level == "I1", ]
  expect_equal(i1$RetainedOccurrences, 2L)
  expect_equal(i1$LostLinks, 1L) # I2 still supports A-B
  expect_equal(i1$ComponentsAfter, 3L)
  expect_equal(i1$NewlyDisconnectedPairs, 2L) # A-C and B-C; D excluded
  expect_equal(d$lost_pairs$From, c("A", "B"))
  expect_equal(d$lost_pairs$To, c("C", "C"))
  expect_equal(d$lost_links$LinkID, "2")
  expect_equal(d$removal$NewlyDisconnectedPairs[d$removal$Level == "I2"], 0L)
  expect_equal(nrow(d$components_after), 8L)
  expect_match(p$data$notes$Text[nrow(p$data$notes)], "No offsets, estimates or SEs")
})

test_that("alternative paths preserve connectivity after a direct link disappears", {
  x <- graph_chain()
  x$links <- data.frame(Link = 1:3, FromID = c(1L, 2L, 1L), ToID = c(2L, 3L, 3L),
    From = c("A", "B", "A"), To = c("B", "C", "C"), LinkSupportAdequate = TRUE)
  x$element_detail <- data.frame(LinkID = 1:3, Facet = "Item", Level = c("I1", "I2", "I3"),
    Retained = TRUE, Flag = FALSE)
  d <- plot(x, type = "anchor_removal", draw = FALSE)$data$data
  expect_equal(d$removal$LostLinks, rep(1L, 3))
  expect_equal(d$removal$NewlyDisconnectedPairs, rep(0, 3))
  expect_equal(nrow(d$lost_links), 3L)
  expect_equal(nrow(d$lost_pairs), 0L)
  expect_equal(d$removal$ComponentsAfter, rep(1L, 3))
})

test_that("one shared element can disconnect every baseline pair", {
  x <- graph_chain()
  x$element_detail <- x$element_detail[c(1, 3), ]
  x$element_detail$Retained <- TRUE
  d <- plot(x, type = "anchor_removal", draw = FALSE)$data$data
  expect_equal(nrow(d$removal), 1L)
  expect_equal(d$removal$LostLinks, 2L)
  expect_equal(d$removal$ComponentsAfter, 3L)
  expect_equal(d$removal$NewlyDisconnectedPairs, 3)
  expect_equal(nrow(d$lost_pairs), 3L)
})

test_that("topology retains unknown and flagged provenance without asserting support", {
  x <- graph_chain()
  x$element_detail$Retained <- c(TRUE, FALSE, NA)
  x$element_detail$Flag <- c(TRUE, TRUE, NA)
  p <- plot(x, type = "links", draw = FALSE)
  expect_equal(p$data$data$links$RetainedElements, c(1L, 0L))
  expect_equal(p$data$data$links$UnknownRetention, c(0L, 1L))
  expect_equal(p$data$data$links$Status, c("retained_review", "no_retained"))
  expect_equal(p$data$data$summary$ConnectedPairs, 1L)
  expect_true(any(grepl("Retention is missing", p$data$notes$Text)))
  d <- plot(x, type = "anchor_removal", draw = FALSE)$data$data
  expect_equal(d$removal$NewlyDisconnectedPairs, c(1, 0))
  expect_equal(d$removal$RetainedOccurrences, c(1L, 0L))
})

test_that("empty chains and duplicate labels remain valid topology inputs", {
  x <- graph_chain(rep("same -> name", 3))
  x$element_detail$Retained <- TRUE
  d <- plot(x, type = "anchor_removal", draw = FALSE)$data$data
  expect_equal(d$summary$Waves, 3L)
  expect_equal(d$lost_pairs$FromID, c(1L, 2L))
  expect_equal(d$lost_pairs$ToID, c(3L, 3L))
  x$links <- x$element_detail <- data.frame()
  p <- plot(x, type = "anchor_removal", draw = FALSE)
  expect_equal(p$data$data$summary$Components, 3L)
  expect_equal(p$data$data$summary$ConnectedPairs, 0)
  expect_equal(nrow(p$data$data$removal), 0L)
  expect_equal(nrow(p$data$data$lost_pairs), 0L)
  grDevices::pdf(NULL, width = 5, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_warning(plot(x, type = "anchor_removal"))
  expect_no_warning(plot(x, type = "links"))
})

test_that("new topology views preserve clean payloads and graphics state", {
  grDevices::pdf(NULL, width = 10, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  x <- graph_chain()
  for (type in c("links", "anchor_removal")) {
    p <- plot(x, type = type, draw = FALSE, preset = "monochrome",
      show_title = FALSE, show_notes = FALSE)
    old <- graphics::par(c("mar", "bg", "xpd"))
    expect_no_warning(q <- plot(x, type = type, preset = "monochrome",
      show_title = FALSE, show_notes = FALSE))
    expect_identical(q, p)
    expect_equal(graphics::par(c("mar", "bg", "xpd")), old)
    expect_true(nrow(p$data$notes) > 0L)
    expect_identical(mfrmr:::.mfrm_plot_display_title(p$data), "")
    rgb <- grDevices::col2rgb(p$data$styles$Color)
    expect_equal(rgb[1, ], rgb[2, ])
    expect_equal(rgb[2, ], rgb[3, ])
  }
  expect_equal(graphics::par("mfg")[1:2], c(1L, 2L))
})

test_that("narrow topology plots warn before replacing wide labels with IDs", {
  x <- graph_chain(c(strrep("W", 80), strrep("M", 80), "C"))
  x$element_detail$Facet <- strrep("W", 40)
  grDevices::pdf(NULL, width = 3, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(p <- plot(x, type = "links", show_title = FALSE, show_notes = FALSE),
    "wave IDs are shown")
  expect_identical(p$data$data$nodes$Label, x$config$waves)
  expect_warning(p <- plot(x, type = "anchor_removal", show_title = FALSE, show_notes = FALSE),
    "anchor IDs are shown")
  expect_true(all(p$data$data$removal$Facet == strrep("W", 40)))
})
