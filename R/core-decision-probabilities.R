# Cut scores are supplied on the scoring scale, never inferred from the data.
mfrm_validate_cut_scores <- function(cut_scores) {
  if (is.null(cut_scores)) return(NULL)
  if (!is.numeric(cut_scores) || is.complex(cut_scores) || is.object(cut_scores) ||
      !is.null(dim(cut_scores)) || !length(cut_scores) || any(!is.finite(cut_scores)) ||
      any(diff(cut_scores) <= 0))
    stop("`cut_scores` must be a finite, strictly increasing numeric vector.", call. = FALSE)
  labels <- names(cut_scores)
  if (!is.null(labels) && (anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)))
    stop("Names of `cut_scores` must be distinct, nonempty labels.", call. = FALSE)
  cut_scores
}

# Integrate bands directly: subtracting two CDFs loses narrow upper-tail mass.
# Reuse the caller's normalization and split wide intervals at its density mode.
mfrm_cut_probabilities <- function(integral, cut_points, mass,
                                   lower = -Inf, upper = Inf, mode = 0) {
  segment <- function(a, b) {
    if (a == b) return(0)
    if (a < mode && b > mode) integral(a, mode) + integral(mode, b) else integral(a, b)
  }
  edges <- pmax(lower, pmin(upper, c(lower,cut_points,upper)))
  # Integrate toward the nearest tail, avoiding quadrature over enormous finite
  # intervals. Retain both tails so upper-tail bands avoid CDF cancellation.
  tails <- vapply(edges, function(x) {
    if (x==lower) return(c(0,mass))
    if (x==upper) return(c(mass,0))
    if (x<=mode) {
      left <- integral(lower,x); c(left,mass-left)
    } else {
      right <- integral(x,upper); c(mass-right,right)
    }
  },numeric(2))
  above <- tails[2,-c(1L,length(edges))]/mass
  bands <- vapply(seq_len(length(edges)-1L),function(i) {
    a <- edges[i]; b <- edges[i+1L]
    if (a==b) return(0)
    if (b<=mode) {
      value <- tails[1,i+1L]-tails[1,i]; size <- tails[1,i+1L]
    } else if (a>=mode) {
      value <- tails[2,i]-tails[2,i+1L]; size <- tails[2,i]
    } else {
      value <- mass-tails[1,i]-tails[2,i+1L]; size <- mass
    }
    # Direct integration preserves narrow bands when subtraction would lose
    # significant digits. Wide bands use the already computed tail masses.
    if (value < size/2) value <- segment(a,b)
    value/mass
  },0)
  # Callers integrate at relative tolerance <= 1e-9; this checks accumulated
  # absolute error, not relative accuracy of probabilities below that scale.
  if (any(!is.finite(c(above, bands))) || any(c(above, bands) < 0) ||
      any(c(above, bands) > 1 + 1e-9) || abs(sum(bands)-1) > 1e-9 ||
      any(diff(above) > 1e-9) ||
      max(abs(above - rev(cumsum(rev(bands)))[-1L])) > 1e-9)
    stop("Conditional decision probability integration is unresolved.", call. = FALSE)
  stats::setNames(pmin(1, c(above, bands)),
    c(paste0("Above", seq_along(above)), paste0("Band", seq_along(bands))))
}

mfrm_decision_table <- function(probabilities, persons, cut_scores, prior_mean,
                                prior_sd, status, reason = "") {
  k <- length(cut_scores); n <- length(persons); columns <- 2L*k+1L
  labels <- names(cut_scores) %||% paste0("Cut", seq_len(k))
  if (!identical(dim(probabilities), c(n, columns)))
    stop("Decision probabilities do not match the scoring roster and cuts.", call. = FALSE)
  data.frame(Person = rep(persons, each=columns),
    Type = rep(c(rep("above_cut", k), rep("ability_band", k+1L)), n),
    Label = rep(c(labels, paste0("Band", seq_len(k+1L))), n),
    Lower = rep(c(unname(cut_scores), -Inf, unname(cut_scores)), n),
    Upper = rep(c(rep(Inf, k), unname(cut_scores), Inf), n),
    Probability = as.vector(t(probabilities)),
    PriorMean = rep(rep_len(prior_mean,n), each=columns),
    PriorSD = rep(rep_len(prior_sd,n), each=columns),
    Status = rep(rep_len(status,n), each=columns),
    Reason = rep(rep_len(reason,n), each=columns),
    Uncertainty = rep("conditional_on_point_calibration_and_prior", n*columns),
    row.names = NULL)
}

# Saved probabilities must remain attached to their cuts, priors and roster.
# This checks the saved result; it never repeats posterior scoring.
mfrm_validate_decisions <- function(x) {
  if (!is.list(x)) return(invisible(x))
  cuts <- x$settings$cut_scores; d <- x$decisions
  if (is.null(cuts) && is.null(d)) return(invisible(x))
  cuts <- mfrm_validate_cut_scores(cuts)
  fail <- function() stop("Saved decision probabilities must match their cuts, prior and scoring roster.", call. = FALSE)
  columns <- 2L*length(cuts)+1L
  required <- c("Person","Type","Label","Lower","Upper","Probability","PriorMean",
    "PriorSD","Status","Reason","Uncertainty")
  if (is.null(cuts) || !is.data.frame(d) || !all(required %in% names(d)) ||
      nrow(d) %% columns || anyNA(d[setdiff(required,"Probability")])) fail()
  persons <- unique(d$Person)
  if (nrow(d) != length(persons)*columns) fail()
  first <- match(persons,d$Person)
  expected <- mfrm_decision_table(matrix(NA_real_,length(persons),columns),persons,cuts,
    d$PriorMean[first],d$PriorSD[first],d$Status[first],d$Reason[first])
  for (name in setdiff(required,"Probability")) if (!identical(d[[name]],expected[[name]])) fail()
  if (any(!is.finite(d$PriorMean)) || any(!is.finite(d$PriorSD)) || any(d$PriorSD < 0) ||
      any(!d$Status %in% c("available_conditional","prior_only","review_only","unavailable"))) fail()
  scores <- x$estimates %||% x$table %||% x$scores
  if (is.data.frame(scores) && "Person" %in% names(scores)) {
    if (!all(scores$Person %in% persons)) fail()
    if (all(c("PriorMean","PriorSD") %in% names(scores))) {
      at <- match(scores$Person,d$Person)
      if (!identical(as.numeric(scores$PriorMean),d$PriorMean[at]) ||
          !identical(as.numeric(scores$PriorSD),d$PriorSD[at])) fail()
    }
  }
  for (person in persons) {
    rows <- d[d$Person==person,]; p <- rows$Probability
    if (rows$Status[1L]=="unavailable") {
      if (!all(is.na(p))) fail()
    } else {
      if (rows$PriorSD[1L] <= 0 || any(!is.finite(p)) || any(p<0 | p>1)) fail()
      bands <- p[-seq_along(cuts)]
      if (abs(sum(bands)-1)>1e-9 ||
          max(abs(p[seq_along(cuts)]-rev(cumsum(rev(bands)))[-1L]))>1e-9) fail()
    }
  }
  invisible(x)
}

mfrm_decision_note <- function() {
  paste("Decision probabilities describe Lower < ability <= Upper on the calibration scale.",
    "Calibration, population estimates and supplied cut scores are held fixed.",
    "Probability accuracy depends on the scoring prior matching the target population; numerical checks do not validate that match.",
    "Uncertainty in an estimated prior is excluded, and a larger calibration sample does not correct a mismatched scoring prior.",
    "They are not Person-difference tests or probabilities of passing a future administration.",
    "Cut labels do not establish CEFR alignment. Very small tails have absolute, not relative, numerical accuracy.")
}
