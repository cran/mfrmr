#' G-study for one or several observed scores
#'
#' Investigate how observed assessment scores vary across people, raters and
#' tasks. With several score columns, such as fluency and accuracy, also study
#' how their sources of variation relate. Use [mfrm_d_study()] next
#' to compare plans with different numbers of raters or tasks.
#'
#' Estimate observed-score variance-covariance components for fixed score
#' components with one or two common random measurement facets, such as
#' raters, tasks, or occasions. Use balanced ANOVA for a complete design or
#' explicitly select MINQUE(0) for incomplete or unequal observed designs. Conditions
#' have common identities across persons and scores. This is not an MFRM fit.
#'
#' @param data A data frame with at most one row per observed combination of
#'   person and selected facet identifiers. Do not add rows for unassigned
#'   cells or code absent scores as zero.
#' @param scores Names of finite numeric score columns, in the desired order.
#'   A single column is allowed as the univariate special case. Scores are
#'   neither standardized nor converted from category labels.
#' @param person,rater,task Distinct columns identifying persons and facets.
#'   By default both rater and task facets are included. Set `rater = NULL`
#'   for Person-by-Task, or `task = NULL` for Person-by-Rater. At least one
#'   facet is required. Each included factor must have at least two observed levels. Labels
#'   may be character, factor, finite numeric, or logical. Blank or infinite
#'   labels fail; missing labels follow `missing`.
#' @param method `"anova"` (default) requires a complete balanced design.
#'   `"minque0"` estimates the same covariance components from a complete or
#'   incomplete or unequal design using identity-working-covariance MINQUE. Neither
#'   method constrains covariance estimates to be positive semidefinite.
#' @param missing `"error"` (default) refuses missing selected scores or
#'   factor identifiers. `"omit"` explicitly excludes any such row from
#'   every score's analysis and retains exclusion accounting. It does not
#'   impute values or correct missing-data bias. Infinite/nonnumeric scores
#'   are always refused.
#' @param facets Optional character vector selecting one or two common random
#'   facets instead of `rater`/`task`, for example
#'   `c(Rater = "Assessor", Occasion = "Session")`. Values identify data
#'   columns; names label covariance components, D-study count columns, and
#'   plots. An unnamed vector uses the column names as labels. Labels must be
#'   unique and nonblank, contain no `:`, and not use reserved result names
#'   such as `Person`, `Residual`, `Score`, `Scenario`, or plotting fields.
#'   Do not supply `rater` or `task` together with `facets`. `person` and
#'   `scores` remain explicit. `NULL` preserves the task/rater interface.
#' @param nesting Optional named character vector specifying one measurement
#'   facet nested within the other, for example `c(Rater = "Task")`.
#'   The name is the child facet label and the value is its parent label;
#'   with `facets`, use its labels rather than input column names. Persons
#'   remain crossed with these conditions. Child identities are local to each
#'   parent: R1 within Task 1 differs from R1 within Task 2. `NULL` (default)
#'   specifies crossed facets. See "Nested measurement facets" below.
#'
#' @details Every retained cell must have every selected score. Duplicate
#'   cells are not supported. The same identifier (parent/child pair for a
#'   nested child) must denote the same condition across persons and scores.
#'   Equal level counts alone cannot establish that identity or random sampling
#'   from the intended universe.
#'
#'   In 0.2.4, `mfrm_g_study()` replaces the longer public name
#'   `mfrm_multivariate_gstudy()` for both one and several scores. The arguments,
#'   default estimator and result class are unchanged. Use [mfrm_d_study()]
#'   for future plans. Existing saved G-study objects remain usable without
#'   re-estimation; update constructor calls in old scripts. See
#'   [mfrmr_compatibility_layer] for migration from both former names.
#'
#'   Choose facets from the intended use of scores: common raters for judging
#'   performances, common tasks for sampling content, or common occasions for
#'   repeat assessments. All selected facets are random and, unless `nesting`
#'   is supplied, fully crossed; `method = "minque0"` permits incomplete observations
#'   of that model. Naming a column Occasion does not model growth, practice,
#'   trends, or serial correlation. Its levels must be defensibly treated as
#'   exchangeable conditions under the stated random-effects assumptions.
#'   Score columns are fixed components of the assessment, not another random
#'   facet. Three or more facets, fixed measurement facets, nesting within
#'   persons, and partial sharing across score components are not supported.
#'
#'   MINQUE(0) models a common mean for each score and independent, zero-mean
#'   random effects with a common covariance matrix for each component. With
#'   intercept-removal matrix `H` and shared-level covariance kernels `K_s`, it
#'   solves `S_st = tr(H K_s H K_t)` against `Q_s = Y' H K_s H Y` for every
#'   score pair. Group-count calculations avoid an observation-by-observation
#'   matrix or a full Cartesian grid. On balanced data it agrees with ANOVA.
#'   It does not optimize a likelihood, fit covariate-dependent means, or
#'   iteratively estimate working covariance weights.
#'
#'   The diagonally scaled moment system must have all eigenvalues above
#'   `sqrt(.Machine$double.eps)` times its largest eigenvalue. Otherwise the
#'   function stops because the components cannot be separated reliably by
#'   these equations. `estimation` retains the rank, scaled eigenvalues,
#'   condition number, and per-component replication counts. These are design
#'   and numerical diagnostics, not precision estimates or model-fit tests.
#'   Graph connectedness alone does not establish component identifiability.
#'   For example, four tasks scored by two distinct raters per performance
#'   and eight tasks scored once have the same per-person workload. Only the
#'   former provides repeated ratings within a Person/Task cell; with one
#'   rating per cell, Person-by-Task and residual variation cannot be separated
#'   in the seven-component model. Other overlaps are still needed for the
#'   remaining components. Uneven workloads across raters can also affect
#'   precision even when the moment system has full rank.
#'   [vcov.mfrm_multivariate_gstudy()] supplies normal-model sampling
#'   covariance and SEs for the original component estimates, including
#'   incomplete sources. Inspect its separate negative/nonidentifiable
#'   sampling-variance diagnostics; it does not provide normal-reference intervals.
#'
#'   Conditioning on the observed assignments must preserve the stated
#'   random-effect means and covariances. Outcome-dependent assignment or
#'   missingness may violate this assumption; calling missingness MAR does
#'   not correct omitted covariates or selection. Review planned assignments
#'   separately from recorded scores, for example with [describe_mfrm_data()]
#'   and its `expected_design` argument. `design$observed_fraction` uses the
#'   cross-product of retained observed levels for crossed facets, or persons
#'   times the number of observed parent/child pairs for nested facets. It is not an assignment
#'   completion rate and cannot reveal entirely unobserved persons or facets.
#'
#'   With crossed rater and task facets, the seven components are `Person`, `Rater`, `Task`, `Person:Rater`,
#'   `Person:Task`, `Rater:Task`, and `Residual`. With one observation per cell,
#'   `Residual` combines the three-way interaction and within-cell error;
#'   they cannot be separated. With one facet, the three components are
#'   `Person`, the facet, and `Residual`; the last combines Person-by-facet
#'   interaction and within-cell error. Labels supplied in `facets` replace
#'   Rater/Task in these component names. Omitting a facet does not remove its
#'   effects from scores or support generalization to new conditions of that facet. No
#'   scores are averaged automatically. Each component is a matrix whose
#'   diagonal contains variances and whose off-diagonal contains covariances between
#'   scores. For ANOVA, mean products replace the mean squares used for a
#'   single score.
#'
#'   Raw component estimates, including negative variances and indefinite
#'   matrices, are retained without clipping or nearest-PSD repair. The
#'   `component_diagnostics` assess eigenvalues after scaling by the observed
#'   score SDs, using a relative numerical tolerance. This is a numerical
#'   admissibility check, not a significance test or precision assessment.
#'   Rank-deficient components are reported. [mfrm_d_study()]
#'   retains raw projected variances and evaluates each score/composite and
#'   metric separately. A non-PSD component does not automatically withhold
#'   coefficients: for example, Rater-by-Task does not enter relative error.
#'   D-study tables and plots flag non-PSD components even when a requested
#'   projection can be calculated. Inspect these diagnostics before using it.
#'   Non-PSD estimates can arise from sampling variation even under a correctly
#'   specified model, particularly when a true component is small or zero.
#'   They do not by themselves prove bad data, informative missingness, or
#'   model misspecification. Automatic clipping or deleting a component would
#'   change the estimation procedure and is not performed.
#'   Distinguish `estimation$rank`, the number of separable components, from
#'   `component_diagnostics$Rank`, the rank of each between-score covariance
#'   matrix. A true zero covariance component can have matrix rank zero
#'   while the design still separates that component from the others.
#'
#'   The model concerns numeric observed scores. Treating ordered categories
#'   as numeric does not estimate latent ordinal or MFRM reliability. The
#'   The FACETS G-theory web example gives different coefficients for
#'   main-effects and interaction models (about 0.72930 versus 0.55144 in its
#'   nested example). Match the specified decomposition before comparing.
#'   Its Rater and Task:Rater terms share the same nested incidence and their
#'   summed variance corresponds to Rater(Task) here. A matching coefficient
#'   does not identify those redundant terms separately.
#'   GENOVA 3.1 and urGENOVA 2.1 balanced nested examples agree with the
#'   corresponding ANOVA variance components; mGENOVA 2.1 also checks
#'   multivariate components and D-study coefficients. urGENOVA uses
#'   Henderson's Method 1 for unbalanced designs, whereas this function offers
#'   MINQUE0 there; equal results are not required for different estimators.
#'   urGENOVA has no D-study output. These comparisons do not establish
#'   unbalanced interval coverage or general interchangeability. The
#'   mGENOVA Appendix E example checks balanced numerical calculations;
#'   MINQUE(0) additionally agrees with direct covariance-kernel calculations
#'   for incomplete designs. These checks do not establish population recovery
#'   for arbitrary sparse assignments or missingness mechanisms. This function
#'   returns point estimates, not sampling intervals.
#'
#' @section Component uncertainty versus measurement error:
#'   For a balanced ANOVA result, [confint.mfrm_multivariate_gstudy()] adds
#'   sampling SEs and pointwise intervals for variance components of one score
#'   or a prespecified weighted composite, under an explicit normal-effects
#'   assumption. `plot(g, intervals = TRUE, assumption = "normal")` adds
#'   these intervals to a component plot; its returned data retain SEs,
#'   raw endpoints and reasons for unavailable intervals.
#'   For balanced or identifiable incomplete sources, explicit
#'   `method = "reml_profile"` supplies approximate variance or two-score covariance intervals
#'   from a separate normal REML model. Original estimates remain unchanged,
#'   with the REML estimate in `IntervalCenter`. It does not supply a component
#'   sampling SE, and numerical availability is not a coverage guarantee.
#'   For between-score covariance, pass exactly two score names to `confint()`
#'   and explicitly select `method = "gci"` or `"regression_region"` for balanced ANOVA,
#'   or experimental `"reml_profile"` for the joint normal model.
#'   These retain signed bounds and distinguish approximate generalized coverage
#'   from conservative normal-model enclosures. `plot(g, type = "covariance",
#'   score = c("V", "W"), intervals = TRUE, method = "gci",
#'   assumption = "normal")` shows their intervals across sources.
#'   Explicit `method = "residual_contrast"` also accepts retained ANOVA/MINQUE
#'   ratings, but supplies only Residual intervals: exact normal-model
#'   chi-square for variance, or a conservative covariance enclosure. It
#'   removes all nonresidual incidence columns, requires residual degrees of
#'   freedom, retains original points and reports a separate `IntervalCenter`.
#'   It does not supply the original component estimator's SE or intervals
#'   for G/Phi or other sources. See [confint.mfrm_multivariate_gstudy()].
#'   `score_scale` contains observed score SDs used for
#'   numerical diagnostics; these are not component SEs. Neither matrix rank
#'   nor a nonnegative estimate establishes that a component is precisely known.
#'   A D-study's SEM instead describes projected measurement error in score
#'   units. The separate [mfrm_multivariate_d_compare()] provides approximate
#'   intervals for **plan differences** in its supported two-facet crossed
#'   scope; these are not component intervals or marginal intervals for each
#'   plan's coefficient. See [mfrmr_methods_reference] for the target map.
#'
#' @section Fixed tasks as score components:
#'   To plan ratings of the same fixed interview, presentation and discussion,
#'   put the three task scores in separate columns and use `task = NULL`.
#'   Each row is one Person/Rater pair; the same rater and person identities
#'   must apply to every task column. The tasks and prespecified score weights
#'   define the fixed composite. Only raters are sampled measurement conditions.
#'   Different task-specific rater teams do not satisfy this representation.
#'
#'   The Person covariance includes stable Person-by-fixed-task differences.
#'   For weight vector `w` and `n_r` planned raters, universe variance is
#'   `w' P w`, relative error is `w' E w / n_r`, and absolute error is
#'   `w' (R + E) w / n_r`. These match a univariate Person-by-Rater analysis
#'   of the directly weighted task score. Do not average task-specific
#'   reliability coefficients or additionally divide error by the task count.
#'   Each planned rater scores every fixed task. Two raters for three fixed
#'   tasks require six ratings per person, not two or three.
#'
#'   MINQUE(0) permits identifiable incomplete Person/Rater source designs;
#'   every retained row still needs all task scores. Explicit `missing = "omit"`
#'   removes a whole incomplete score vector. Do not impute unassigned tasks
#'   to manufacture common score identities. The D-study projects a future
#'   complete common-rater design; it does not estimate sparse-roster reliability.
#'   Adding or replacing tasks, partial rater sharing and an arbitrary mixture
#'   of fixed/random facets are not implemented by this representation. Holding
#'   the task count constant in a random-task model is a different assumption.
#'
#' @section Nested measurement facets:
#'   Suppose each task has its own rater team, and each team rates the same
#'   persons on both Content and Organization. Use `nesting = c(Rater = "Task")`
#'   for Person crossed with Rater-within-Task. The five components are
#'   `Person`, `Task`, `Rater(Task)`, `Person:Task`, and `Residual`. The last
#'   combines Person-by-Rater-within-Task interaction and within-cell error.
#'   There is no separately estimated common Rater or Rater-by-Task component.
#'   A rater who actually works across tasks is not an independent nested
#'   rater; changing their identifier cannot establish independence.
#'
#'   ANOVA requires all persons at every observed parent/child combination
#'   and the same number of children per parent. MINQUE(0) also permits missing
#'   cells and unequal child counts when its moment equations separate the
#'   five components. Having only one child per parent confounds parent and
#'   child components. Neither method corrects informative assignments.
#'
#'   `design$counts["Rater"]` is the total number of distinct task/rater
#'   pairs; `design$child_counts` gives the rater count within each task.
#'   `design$levels` retains the original labels, interpreted locally for the
#'   child facet. A design is `complete` when every person has every observed
#'   nested condition, and `balanced` when it also has equal child counts.
#'   In a future D-study, `Raters` instead means raters per task: two raters
#'   for each of six tasks use twelve distinct raters and twelve ratings per
#'   person. Future designs preserve this nesting and have equal child counts.
#'   Other two-facet names and the opposite nesting direction follow the same
#'   rules; nesting within persons and score-specific child identities do not.
#'
#' @section Reviewing incomplete designs:
#'   Begin with the planned roster and `data_usage` to distinguish unassigned
#'   cells from missing assigned ratings. Then inspect `estimation` for
#'   component separation and replication, `component_diagnostics` for
#'   covariance admissibility, and the D-study metric-specific status columns
#'   (`GStatus`, `PhiStatus`, `RelativeSEMStatus`, `AbsoluteSEMStatus`) before
#'   reading the corresponding result. None of these checks estimates sampling
#'   precision or identifies the missingness mechanism. No fixed percentage of observed cells,
#'   condition-number cutoff, or passed matrix check establishes adequate
#'   precision. D-studies project future complete balanced designs. The separate
#'   [mfrm_multivariate_d_compare()] offers approximate normal-theory intervals
#'   for prespecified differences in two-facet crossed designs; nested-design
#'   intervals are not currently available.
#'
#'   For the same rating workload, changing rater overlap or concentrating
#'   assignments can change estimation precision. Omitting ratings selected
#'   by their scores can introduce bias even when every requested coefficient
#'   is calculable. Neither a returned coefficient nor agreement with another
#'   program determines whether these assumptions suit the user's assessment.
#'   Observed pool sizes concern estimation in the G-study; they differ from
#'   the per-person counts in a future D-study scenario. Review which tasks
#'   and raters are shared and whether recorded ratings represent the intended
#'   population and conditions before interpreting a projected improvement.
#'
#' @return An `mfrm_multivariate_gstudy` list with `components` (three, five, or seven
#'   named covariance matrices), ANOVA `mean_products` and `degrees_of_freedom`
#'   (`NULL` for MINQUE), `estimation` (MINQUE moment matrices and diagnostics,
#'   `NULL` for ANOVA), `data_usage` (input/used/excluded counts, excluded input
#'   row positions and missing cells),
#'   `component_diagnostics`, `score_scale` (observed SDs used only for matrix
#'   diagnostics), `design` (`facets` maps facet labels to input columns and
#'   `count_columns` maps those labels to D-study columns; also levels, counts, method, and
#'   score convention, `nesting`, `child_counts`, completeness, balance,
#'   potential cells and observed cell fraction), and `data`
#'   (the retained selected input columns).
#'   `print()` displays the design and diagnostics and returns the same object
#'   invisibly. `summary()` returns the component-diagnostics data frame;
#'   use `$components` for the estimated covariance matrices.
#' @references Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Chapters 9--11. \doi{10.1007/978-1-4757-3456-0}.
#'
#'   Brennan, R. L. (2001). *Manual for mGENOVA, Version 2.1*.
#'   Iowa Testing Programs Occasional Papers, No. 50. Pages 7--8 and 19--22;
#'   Table 12 (page 32) and Appendix E (pages 74--77).
#'   Manual available from the author's
#'   [GENOVA Suite](https://brennancrickgenova.org/genova-suite/).
#'
#'   Brennan, R. L. (1992). Generalizability theory.
#'   *Educational Measurement: Issues and Practice*, 11(4), 27--34.
#'   Equations 13--16 and Table 3.
#'
#'   Rao, C. R. (1971). Estimation of variance and covariance components--MINQUE
#'   theory. *Journal of Multivariate Analysis*, 1, 257--275.
#'   \doi{10.1016/0047-259X(71)90001-7}.
#' @seealso [mfrm_d_study()], [mfrm_generalizability()]. For a
#'   guided data-to-plan-to-saved-result workflow, use
#'   `vignette("mfrmr-gtheory", package = "mfrmr")`.
#' @examples
#' # Common tasks, without a rater facet: Brennan's published synthetic data.
#' # mGENOVA manual Table 12: 10 persons, 6 common items, two scores V and W.
#' # The item facet is named Task in the supplied long-format data.
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g_task <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#' g_task$components
#' g_task$component_diagnostics
#' # To compare future designs, see ?mfrm_d_study.
#'
#' # Optional planning, incomplete designs, renamed facets and nesting examples.
#' # example("mfrm_g_study", run.dontrun = TRUE)
#' \dontrun{
#' d_task <- mfrm_d_study(g_task,
#'   design_grid = data.frame(Tasks = c(6, 12)), weights = c(V = -1, W = 1))
#' d_task$coefficients
#' # At six tasks, W - V has G = 0.30000 and Phi = 0.24116 (Appendix E).
#'
#' # An illustrative incomplete assignment roster, not a missingness model.
#' sparse <- tasks[(tasks$Person + tasks$Task) %% 3 != 0, ]
#' sparse$V[1] <- NA_real_ # One assigned score is additionally unrecorded.
#' g_sparse <- mfrm_g_study(sparse, c("V", "W"), rater = NULL,
#'   method = "minque0", missing = "omit")
#' g_sparse$data_usage
#' g_sparse$estimation$component_support
#' g_sparse$component_diagnostics
#' # Explicit future COMPLETE designs, not reliability of the sparse roster.
#' d_sparse <- mfrm_d_study(g_sparse,
#'   data.frame(Tasks = c(6, 12)), weights = c(V = -1, W = 1))
#' d_sparse$coefficients # Read metric status and ComponentPSD separately.
#'
#' # Existing continuous-score example; rename columns to show role mapping.
#' data("mfrmr_example_gtheory_crossed", package = "mfrmr")
#' ratings <- mfrmr_example_gtheory_crossed
#' names(ratings)[names(ratings) == "Rater"] <- "Assessor"
#' names(ratings)[names(ratings) == "Task"] <- "Session"
#' head(ratings) # Two score columns for each Person-Assessor-Session cell.
#' # Here sessions represent exchangeable assessments, without a time trend.
#' g <- mfrm_g_study(ratings, c("Content", "Organization"),
#'   facets = c(Rater = "Assessor", Occasion = "Session"))
#' g$components$Person
#' g$component_diagnostics
#' d <- mfrm_d_study(g,
#'   design_grid = data.frame(Rater = c(2, 4), Occasion = c(3, 3)),
#'   weights = c(Content = 0.6, Organization = 0.4))
#' d$coefficients
#' # Difference-score dependability, when subtraction is meaningful on these scales.
#' difference <- mfrm_d_study(g,
#'   weights = c(Content = 1, Organization = -1))
#' difference$coefficients
#'
#' # Rater-only planning for one occasion; no generalization across occasions.
#' one_session <- ratings[ratings$Session == 1, ]
#' g_rater <- mfrm_g_study(one_session, c("Content", "Organization"),
#'   rater = "Assessor", task = NULL)
#' d_rater <- mfrm_d_study(g_rater, data.frame(Raters = c(1, 2, 4)))
#' plot(d_rater, score = "Content")
#'
#' # Different rater teams for different tasks; the same Persons and scores.
#' data("mfrmr_example_gtheory_nested", package = "mfrmr")
#' nested <- mfrmr_example_gtheory_nested
#' head(nested) # R1 in Task 1 is a different rater from R1 in Task 2.
#' g_nested <- mfrm_g_study(nested, c("Content", "Organization"),
#'   nesting = c(Rater = "Task"))
#' g_nested$design$child_counts # Three raters per task; eighteen in total.
#' g_nested$component_diagnostics
#' d_nested <- mfrm_d_study(g_nested,
#'   expand.grid(Raters = c(2, 3), Tasks = c(4, 6)),
#'   weights = c(Content = 0.6, Organization = 0.4))
#' summary(d_nested)
#' plot(d_nested, x_var = "Raters") # Raters on the axis means raters PER TASK.
#' plot(d_nested, x_var = "Tasks", type = "sem")
#' }
#' @aliases mfrm_multivariate_gstudy
#' @export
mfrm_g_study <- function(data, scores, person = "Person",
                                     rater = "Rater", task = "Task",
                                     method = c("anova", "minque0"),
                                     missing = c("error", "omit"), facets = NULL,
                                     nesting = NULL) {
  method <- match.arg(method)
  missing <- match.arg(missing)
  valid_names <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
  if (!is.data.frame(data) || nrow(data) == 0L || !valid_names(names(data))) {
    stop("`data` must be a nonempty data frame with unique, nonblank column names.", call. = FALSE)
  }
  if (is.null(facets)) {
    facet_args <- list(Rater = rater, Task = task)
    facet_args <- facet_args[!vapply(facet_args, is.null, logical(1))]
    facet_names <- names(facet_args)
    count_columns <- c(Rater = "Raters", Task = "Tasks")[facet_names]
  } else {
    if (!missing(rater) || !missing(task)) {
      stop("Use either `facets` or `rater`/`task`, not both.", call. = FALSE)
    }
    if (!is.character(facets) || !is.null(dim(facets)) || !valid_names(facets)) {
      stop("`facets` must specify one or two distinct identifier columns.", call. = FALSE)
    }
    facet_names <- if (is.null(names(facets))) facets else names(facets)
    if (!.mfrm_mvgt_valid_facets(facet_names)) {
      stop("Supply one or two unique, nonblank facet labels without ':' or reserved result names (such as Person, Residual, Score, or Scenario).", call. = FALSE)
    }
    facet_args <- as.list(unname(facets))
    count_columns <- setNames(facet_names, facet_names)
    rater <- task <- NULL
  }
  one_facet <- length(facet_args) == 1L
  nesting <- .mfrm_mvgt_nesting(nesting, facet_names)
  id_args <- c(list(person), facet_args)
  ids <- unlist(id_args, use.names = FALSE)
  if (!length(facet_args) || any(lengths(id_args) != 1L) ||
      !all(vapply(id_args, is.character, logical(1))) ||
      !valid_names(ids) || !valid_names(scores) || any(scores %in% ids) ||
      !all(c(ids, scores) %in% names(data))) {
    stop("Supply distinct person, facet, and score column names present in `data` (at least one facet).", call. = FALSE)
  }
  data <- as.data.frame(data[c(ids, scores)])
  groups <- lapply(data[ids], function(x) {
    if (!is.null(dim(x)) || !(is.factor(x) || (!is.object(x) &&
        (is.character(x) || is.logical(x) || (is.numeric(x) && !is.complex(x))))) ||
        (is.numeric(x) && any(!is.na(x) & !is.finite(x))) ||
        any(!is.na(x) & !nzchar(trimws(as.character(x))))) {
      stop("Factor identifiers must be finite, nonmissing, nonblank labels.", call. = FALSE)
    }
    factor(x)
  })
  if (!all(vapply(data[scores], function(x) is.numeric(x) && !is.object(x) &&
      !is.complex(x) && is.null(dim(x)) && all(is.na(x) | is.finite(x)), logical(1)))) {
    stop("All selected scores must be finite numeric values or NA; correct invalid scores explicitly.", call. = FALSE)
  }
  absent <- is.na(data)
  keep <- rowSums(absent) == 0L
  cells <- which(absent, arr.ind = TRUE)
  data_usage <- list(source = "data", missing = missing,
    counts = c(InputRows = nrow(data), UsedRows = sum(keep), ExcludedRows = sum(!keep)),
    excluded_rows = which(!keep),
    missing_cells = data.frame(InputRow = cells[, 1L], Column = names(data)[cells[, 2L]],
      row.names = NULL))
  if (any(!keep) && missing == "error") {
    stop("Scores must be finite numeric values and identifiers nonmissing, nonblank; choose missing = 'omit' explicitly to exclude incomplete rows.", call. = FALSE)
  }
  if (!any(keep)) stop("No complete rows remain for the G-study.", call. = FALSE)
  # Refuse duplicate identifiable cells before omission can hide a replicate.
  identified <- rowSums(absent[, ids, drop = FALSE]) == 0L
  if (anyDuplicated(data[identified, ids, drop = FALSE])) {
    stop("Duplicate cells in the selected design are not supported; do not average replicates silently.", call. = FALSE)
  }
  data <- data[keep, , drop = FALSE]
  groups <- lapply(groups, function(g) droplevels(g[keep]))
  factor_names <- c("Person", facet_names)
  counts <- setNames(vapply(groups, nlevels, integer(1)), factor_names)
  if (any(counts < 2L)) stop("Each included factor needs at least two observed levels.", call. = FALSE)
  child_counts <- NULL
  potential_cells <- prod(as.double(counts))
  if (!is.null(nesting)) {
    child <- match(names(nesting), factor_names)
    parent <- match(unname(nesting), factor_names)
    pairs <- unique(data.frame(Parent = groups[[parent]], Child = groups[[child]]))
    child_counts <- setNames(as.integer(table(pairs$Parent)), levels(groups[[parent]]))
    counts[child] <- nrow(pairs)
    potential_cells <- as.double(counts[1L]) * sum(child_counts)
  }
  complete <- nrow(data) == potential_cells
  balanced <- complete && (is.null(nesting) || length(unique(child_counts)) == 1L)
  if (!balanced && method == "anova") {
    stop("ANOVA requires a complete, balanced design; consider method = 'minque0' for an incomplete or unequal observed design.", call. = FALSE)
  }
  if (!is.null(nesting) && method == "anova" && any(child_counts < 2L)) {
    stop("Cannot separate parent and nested-child components with one child per parent.", call. = FALSE)
  }
  y <- as.matrix(data[scores])
  y <- sweep(y, 2L, colMeans(y), "-")
  if (any(!is.finite(y))) stop("Score centering overflowed; rescale the scores.", call. = FALSE)
  score_scale <- sqrt(colSums(y^2) / (nrow(y) - 1))
  if (any(!is.finite(score_scale)) || any(score_scale == 0 & colSums(abs(y)) > 0)) {
    stop("Score variation exceeds numeric range; rescale the scores.", call. = FALSE)
  }
  names(score_scale) <- scores
  # Shared groupings for the crossed components; unique highest-order cells
  # represent the combined interaction/residual kernel.
  subsets <- if (one_facet) list(1L, 2L, 1:2) else
    list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
  if (!is.null(nesting)) subsets <- list(1L, parent, c(child, parent), c(1L, parent), 1:3)
  sources <- .mfrm_mvgt_sources(factor_names, nesting)
  # Integer level codes prevent collisions between literal IDs containing dots.
  groupings <- setNames(lapply(subsets, function(s) .mfrm_mvgt_group(groups[s])), sources)
  estimation <- NULL
  if (method == "minque0") {
    fitted <- .mfrm_mvgt_minque0(y, groupings)
    components <- fitted$components
    estimation <- fitted$estimation
    mean_products <- degrees_of_freedom <- NULL
  } else {
    effects <- mean_products <- setNames(vector("list", length(sources)), sources)
    df <- setNames(numeric(length(sources)), sources)
    for (i in seq_along(subsets)) {
      subset <- subsets[[i]]
      group <- groupings[[i]]
      effect <- (rowsum(y, group) / tabulate(group))[group, , drop = FALSE]
      for (j in seq_len(i - 1L)) {
        if (all(subsets[[j]] %in% subset)) effect <- effect - effects[[j]]
      }
      effects[[i]] <- effect
      df[i] <- if (is.null(nesting)) prod(counts[subset] - 1) else {
        P <- counts[1L]; T <- counts[parent]; R <- child_counts[1L]
        c(P - 1, T - 1, T * (R - 1), (P - 1) * (T - 1), T * (P - 1) * (R - 1))[i]
      }
      mean_products[[i]] <- crossprod(effect) / df[i]
    }
    m <- mean_products
    if (!is.null(nesting)) {
      components <- list((m[[1L]] - m[[4L]]) / (R * T),
        (m[[2L]] - m[[3L]] - m[[4L]] + m[[5L]]) / (P * R),
        (m[[3L]] - m[[5L]]) / P, (m[[4L]] - m[[5L]]) / R, m[[5L]])
    } else if (one_facet) {
      components <- list((m[[1L]] - m[[3L]]) / counts[2L],
        (m[[2L]] - m[[3L]]) / counts[1L], m[[3L]])
    } else {
      # The first and second random facets have the same averaging rules
      # whatever their substantive names (e.g., raters, tasks, occasions).
      n1 <- counts[2L]
      n2 <- counts[3L]
      components <- list(
        (m[[1L]] - m[[4L]] - m[[5L]] + m[[7L]]) / (n1 * n2),
        (m[[2L]] - m[[4L]] - m[[6L]] + m[[7L]]) / (counts[1L] * n2),
        (m[[3L]] - m[[5L]] - m[[6L]] + m[[7L]]) / (counts[1L] * n1),
        (m[[4L]] - m[[7L]]) / n2,
        (m[[5L]] - m[[7L]]) / n1,
        (m[[6L]] - m[[7L]]) / counts[1L], m[[7L]]
      )
    }
    names(components) <- sources
    degrees_of_freedom <- data.frame(Source = sources, DF = unname(df))
  }
  if (any(!is.finite(unlist(components)))) {
    stop("Covariance estimation exceeded numeric range; rescale the scores.", call. = FALSE)
  }
  diagnostics <- .mfrm_mvgt_diagnostics(components, score_scale)
  structure(list(components = components, mean_products = mean_products,
    degrees_of_freedom = degrees_of_freedom, estimation = estimation, data_usage = data_usage,
    component_diagnostics = diagnostics, score_scale = score_scale,
    design = list(person = person, rater = rater, task = task, scores = scores,
      facets = setNames(ids[-1L], facet_names), count_columns = count_columns,
      counts = counts, levels = lapply(groups, levels), rows = nrow(data),
      nesting = nesting, child_counts = child_counts, balanced = balanced,
      complete = complete, potential_cells = potential_cells,
      observed_fraction = nrow(data) / potential_cells,
      method = if (method == "anova") "Balanced multivariate ANOVA" else "Multivariate MINQUE(0)",
      score_convention = if (!is.null(nesting)) paste0("Means over random ",
        unname(nesting), " conditions and their distinct random ", names(nesting),
        " conditions; ", count_columns[[names(nesting)]], " is the count per ", unname(nesting)) else
        if (is.null(facets)) paste("Means over common random",
        paste(tolower(count_columns), collapse = " and ")) else
        paste("Means over common random conditions:", paste(facet_names, collapse = ", ")),
      calculation_version = if (!is.null(nesting)) 3L else if (method == "minque0") 2L else 1L), data = data), class = "mfrm_multivariate_gstudy")
}

.mfrm_mvgt_valid_facets <- function(x) {
  is.character(x) && length(x) %in% 1:2 && !anyNA(x) &&
    all(nzchar(trimws(x))) && !anyDuplicated(x) && !any(grepl(":", x, fixed = TRUE)) &&
    !any(x %in% c("Person", "Residual", "Scenario", "Kind", "Score", "UniverseVariance",
      "RelativeErrorVariance", "AbsoluteErrorVariance", "G", "Phi", "RelativeSEM",
      "GStatus", "PhiStatus", "RelativeSEMStatus", "AbsoluteSEMStatus", "ComponentPSD",
      "AbsoluteSEM", "Status", "Metric", "X", "Group", "Value", "Panel"))
}

.mfrm_mvgt_nesting <- function(nesting, facets) {
  if (is.null(nesting)) return(NULL)
  if (length(facets) != 2L || !is.character(nesting) || !is.null(dim(nesting)) ||
      length(nesting) != 1L || is.na(nesting) || is.null(names(nesting)) ||
      is.na(names(nesting)) || !names(nesting) %in% facets || !nesting %in% facets ||
      names(nesting) == unname(nesting)) {
    stop("`nesting` must name one child facet and its other parent facet, for example c(Rater = 'Task'); nesting within Person is not supported.", call. = FALSE)
  }
  nesting
}

.mfrm_mvgt_sources <- function(factors, nesting = NULL) {
  if (!is.null(nesting)) return(c("Person", unname(nesting),
    paste0(names(nesting), "(", nesting, ")"), paste0("Person:", nesting), "Residual"))
  c(factors, if (length(factors) == 3L) utils::combn(factors, 2L, paste, collapse = ":"),
    "Residual")
}

.mfrm_mvgt_group <- function(groups) {
  as.integer(interaction(lapply(groups, as.integer), drop = TRUE))
}

.mfrm_gt_kernel_gram <- function(groupings, centered = TRUE) {
  # Group counts compute the traces without any n-by-n matrix or full grid.
  n <- as.double(length(groupings[[1L]]))
  k <- length(groupings)
  sizes <- lapply(groupings, tabulate)
  row_sizes <- Map(function(g, size) as.double(size[g]), groupings, sizes)
  totals <- vapply(row_sizes, sum, numeric(1))
  gram <- matrix(0, k, k, dimnames = list(names(groupings), names(groupings)))
  for (i in seq_len(k)) for (j in seq_len(i)) {
    joint_sizes <- tabulate(.mfrm_mvgt_group(list(groupings[[i]], groupings[[j]])))
    value <- sum(as.double(joint_sizes)^2)
    if (centered) value <- value -
      2 * sum(row_sizes[[i]] * row_sizes[[j]]) / n + totals[i] * totals[j] / n^2
    gram[i, j] <- gram[j, i] <- value
  }
  gram
}

.mfrm_mvgt_minque0 <- function(y, groupings) {
  # With H = I - 11'/n and K_s the shared-level covariance kernels, solve
  # S_st = tr(H K_s H K_t), Q_s = Y' H K_s H Y. y is already centered.
  gram <- .mfrm_gt_kernel_gram(groupings)
  k <- length(groupings)
  sizes <- lapply(groupings, tabulate)
  if (any(!is.finite(gram)) || any(diag(gram) <= 0)) {
    stop("The observed design cannot support the requested covariance components.", call. = FALSE)
  }
  scale <- sqrt(diag(gram))
  normalized <- gram / outer(scale, scale)
  eigenvalues <- eigen(normalized, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps)
  rank <- sum(eigenvalues > tolerance * max(eigenvalues))
  if (rank < k) {
    stop("Cannot separate covariance components in this observed design (numerical rank ",
      rank, " of ", k, "). Review facet overlap and confounding; no estimates were returned.", call. = FALSE)
  }
  products <- lapply(groupings, function(g) crossprod(rowsum(y, g)))
  rhs <- do.call(rbind, lapply(products, as.vector))
  estimates <- solve(normalized, rhs / scale) / scale
  components <- setNames(lapply(seq_len(k), function(i) {
    matrix(estimates[i, ], ncol(y), dimnames = list(colnames(y), colnames(y)))
  }), names(groupings))
  list(components = components, estimation = list(
    kernel_gram = gram, quadratic_products = products,
    scaled_eigenvalues = eigenvalues, rank = rank, relative_tolerance = tolerance,
    condition_number = max(eigenvalues) / min(eigenvalues),
    component_support = data.frame(Source = names(groupings),
      Groups = lengths(sizes), MinimumRows = vapply(sizes, min, integer(1)),
      MaximumRows = vapply(sizes, max, integer(1)), row.names = NULL)))
}

.mfrm_mvgt_diagnostics <- function(components, score_scale) {
  scale <- score_scale
  scale[scale == 0] <- 1
  do.call(rbind, lapply(names(components), function(source) {
    standardized <- sweep(sweep(components[[source]], 1L, scale, "/"), 2L, scale, "/")
    values <- eigen(standardized, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(values))
    data.frame(Source = source, MinimumScaledEigenvalue = min(values),
      Tolerance = tolerance, PositiveSemidefinite = min(values) >= -tolerance,
      Rank = sum(values > tolerance), Dimension = length(values), row.names = NULL)
  }))
}

#' Plan measurement conditions with a D-study
#'
#' Would adding tasks or raters make assessment scores more dependable?
#' Use an estimated G-study to compare specified future designs for individual
#' score components and optional weighted composites. Each scenario reports
#' dependability and measurement error for mean scores over common random
#' measurement conditions, such as tasks, raters, or occasions. Projection from
#' an estimated G-study does not refit. Supplying an MFRM first estimates the
#' separate main-effects G-study described below.
#'
#' @param x A result from [mfrm_g_study()]. A main-effects result from
#'   [mfrm_generalizability()] or a fitted MFRM uses the separate
#'   [main-effects method][mfrm_d_study_main_effects], with its collapsed-residual
#'   assumptions. Dispatch preserves the source model; it does not select a new estimator.
#' @param design_grid A nonempty data frame or named list with one positive
#'   integer count column per G-study facet. Both G-study routes accept the
#'   exact facet label or `n_<facet>`, such as `n_Rater` and `n_Task`.
#'   Existing saved count-column names are also accepted: with the `rater`/`task` interface, use `Raters`
#'   and/or `Tasks` for the included facets. With `facets`, use its exact
#'   labels, for example `Rater` and `Occasion` for
#'   `facets = c(Rater = "Assessor", Occasion = "Session")`. The mapping is
#'   available in `x$design$count_columns`.
#'   Supply exactly one column per facet; extra columns, duplicate aliases and
#'   ambiguous matches are refused. Output keeps the saved count-column names,
#'   so existing plot and interval selections retain their meaning.
#'   Each row is one future complete balanced design preserving the G-study's
#'   crossed or nested structure. For `nesting = c(Rater = "Task")`, `Raters`
#'   means raters per task, not the total rater pool. `NULL` uses the G-study
#'   counts (children per parent for a nested facet) only when its retained
#'   data are complete and balanced; incomplete or unequal designs require an
#'   explicit grid. Counts need not match or exceed the G-study counts.
#'   Conditions are shared across persons and scores in every scenario.
#' @param weights Optional named numeric vector for one composite, or a numeric
#'   matrix for several named composites. A vector must name every original
#'   score exactly once; matrix rows must name every score exactly once, and
#'   columns must have unique, nonempty composite names. Rows are matched by
#'   name, regardless of order. Every weight must be finite and each composite
#'   must have at least one nonzero weight. Signed weights allow difference
#'   scores, for example `c(Content = 1, Organization = -1)`. Weights are used
#'   as supplied, without normalization. `NULL` reports original scores only.
#'   Otherwise, each scenario reports the original scores once, followed by
#'   each composite. A vector's composite is named `"Composite"`; a matrix
#'   preserves its column names and order.
#' @param negative_components Component treatment for D-study point projections.
#'   `"raw"` (default) retains the original estimates. `"genova_algorithm"`
#'   replaces negative variance components with zero; `"genova_ems"` first
#'   recalculates components by source-design expected-mean-square equations,
#'   substituting zero as each negative component is encountered. The latter
#'   choices require a complete balanced univariate ANOVA G-study with one or
#'   two supported crossed or nested random facets. Original G-study estimates
#'   remain unchanged. These options do not supply adjusted-estimator intervals.
#'
#' @details Start with a concrete comparison, such as six versus twelve tasks
#'   for every examinee. Each row of `design_grid` is one scenario. With both
#'   facets, `data.frame(Raters = c(2, 4), Tasks = c(3, 3))` compares two versus
#'   four raters at three tasks; `expand.grid(Raters = c(2, 4),
#'   Tasks = c(3, 6))` requests all four combinations. These counts apply to
#'   every person and score, not to the total number of observed ratings.
#'   For `facets = c(Rater = "Assessor", Occasion = "Session")`, use
#'   `expand.grid(Rater = c(1, 2), Occasion = c(2, 4))`. The same workflow
#'   applies to any supported one- or two-facet model; extra columns cannot
#'   introduce conditions absent from the G-study.
#'
#'   This is the primary D-study entry for one or several scores. It replaces
#'   the public `mfrm_multivariate_d_study()` name without changing its result
#'   class, component estimates, projections or uncertainty methods. For either
#'   source model, `summary(d)` returns its full-precision coefficient table.
#'   Optional interval summaries remain specific to the interaction-component
#'   route. Main-effects projections do not acquire those interval methods.
#'
#'   Read `summary(d)` together with `plot(d)`:
#'   * `G` concerns consistency of relative ordering, such as ranking examinees.
#'     It is not the probability of identifying the best person or the winner
#'     of a competition. Close leading scores can remain hard to distinguish
#'     even when G is high; this function provides no Person-rank intervals.
#'   * `Phi` concerns absolute score levels and also counts shifts from easier
#'     tasks or more lenient raters as error. It is not pass/fail accuracy.
#'   * `RelativeSEM` and `AbsoluteSEM` express these errors in the units of the
#'     mean score or composite. Smaller SEMs mean less error. They are not
#'     confidence intervals for G or Phi.
#'
#'   Larger G/Phi and smaller SEM indicate greater dependability under the
#'   model. There is no universally acceptable coefficient: consider the use
#'   of scores, consequences of error, and workload before choosing a design.
#'   Compare a chosen score/composite and metric across feasible plans. For
#'   example, `data.frame(Raters = c(2, 3, 4), Tasks = c(6, 4, 3))` holds the
#'   number of ratings per person at twelve. It does not hold examinee task
#'   burden or total cost constant. The largest projected coefficient is a
#'   point-estimate ranking, not evidence that one plan is reliably better.
#'   Consider the size of the projected difference and which differences
#'   would matter for the assessment. Small rank changes need not imply large
#'   practical losses. If a metric is unavailable for any candidate, comparing
#'   only the remaining values does not resolve the full planning question.
#'   `plot(d, type = "sem")` shows the SEMs. Plots select a sole composite
#'   by default, or the first score when there are no composites. With several
#'   composites, specify one with `plot(d, composite = "Equal")`; use
#'   `plot(d, score = "V")` to inspect an original score named V. The title
#'   metadata identifies the plotted score and any composite weights; it is
#'   returned with the plot, with titles and notes omitted from the figure
#'   by default.
#'   Check `GStatus`, `PhiStatus`, `RelativeSEMStatus`, and `AbsoluteSEMStatus`:
#'   each explains the corresponding metric's availability. `NA` means
#'   unavailable, not zero dependability. `ComponentPSD` separately flags
#'   covariance components needing review. See [plot.mfrm_multivariate_d_study()] for
#'   interpreting curves and unavailable estimates.
#'
#'   For a balanced one-facet ANOVA source, [confint.mfrm_multivariate_d_study()]
#'   returns exact F intervals for G and generalized confidence intervals for
#'   Phi (`parm = "Phi"`) under independent normal random effects, using the
#'   source study's degrees of freedom at each requested future count. Phi
#'   includes Monte Carlo endpoint diagnostics and is not an exact-coverage
#'   method. Add `reference = 1` for intervals on changes from the first plan:
#'   G maps the F set; Phi maps the generalized set and retains approximate
#'   coverage. Both include interior extrema; Phi differences retain joint
#'   numerical endpoint bands. These methods also apply to balanced nested
#'   plans when future child counts per parent equal the source count, using
#'   the parent count in the one-facet reduction. Complete balanced crossed
#'   two-facet sources default to approximate MLS null inversion for individual
#'   G/Phi intervals. Undefined MLS widths are conservatively retained and
#'   flagged by `FormulaComplete`; original estimates and numerical failures
#'   remain explicit. Small-source intervals may span almost all of `[0,1]`.
#'   These sources and general balanced nested plans also have explicitly selected
#'   experimental `method = "approximate_f"` and `"mean_square_region"`
#'   options for G/Phi; a vector requests both and `plot()` compares them.
#'   The confidence-interval help explains their assumptions, limitations and
#'   unavailable results. Use [mfrm_multivariate_d_compare()]
#'   separately for its supported paired-plan differences.
#'   For balanced or identifiable incomplete sources, explicit
#'   `method = "reml_profile"` fits approximate G/Phi or SEM intervals from the used
#'   ratings under normal effects. It retains original ANOVA/MINQUE estimates
#'   and reports the separate REML estimates as `IntervalCenter`. Numerical
#'   diagnostics and experimental finite-sample/boundary qualifications remain
#'   explicit; it is not a default. With `reference`, it supplies shared-profile
#'   interval images for one-facet or common-child nested plans, and direct
#'   difference profiles for other crossed/nested comparisons. The latter need
#'   optional nloptr and retain experimental calibration, especially at zero.
#'   SEM uses a direct error-variance profile followed by square roots; the
#'   numerical interval surrounds the REML center without enumerating remote
#'   disconnected accepted regions. Checked boundary populations showed
#'   below-nominal inclusion despite high numerical availability.
#'   Select `parm = "RelativeSEM"` or `"AbsoluteSEM"` for intervals on SEM
#'   itself from supported balanced crossed/nested sources. A single contributing
#'   mean square uses chi-square limits; otherwise MLS limits for the combined
#'   error variance are transformed by the square root. MLS coverage is
#'   approximate. One-facet SEM differences also accept `reference = 1`;
#'   nested MLS comparisons require a common future child count across plans.
#'   For one-facet RelativeSEM, explicit `method = "residual_contrast"` also
#'   accepts retained ANOVA/MINQUE ratings, including incomplete MINQUE sources.
#'   It transforms one residual
#'   chi-square interval and retains original points, a separate projected
#'   `IntervalCenter`, unavailable states and shared uncertainty across plans.
#'   Normal independent residuals and noninformative omissions are required;
#'   this does not supply AbsoluteSEM, G/Phi or multifacet intervals.
#'   Explicit `method = "mean_square_region"` supplies conservative individual
#'   SEM intervals and G/Phi/SEM difference enclosures for general balanced
#'   crossed/nested two-facet plans. The enclosure retains common source
#'   uncertainty but can be much wider than the tight difference image.
#'   `plot(d, type = "sem", metric = "AbsoluteSEM", intervals = TRUE,
#'   assumption = "normal")` adds the corresponding ribbon without refitting.
#'
#' @section Calculation and interpretation limits:
#'   For crossed covariance component matrices `P`, `R`, `T`, `PR`, `PT`, `RT`,
#'   and `E`, universe-score covariance is `P`. Relative-error covariance is
#'   `PR/n_r + PT/n_t + E/(n_r*n_t)`. Absolute-error covariance adds
#'   `R/n_r + T/n_t + RT/(n_r*n_t)`. The G-study number of persons does not
#'   divide individual-score universe variance. Holding a count constant does
#'   not turn its random facet into a fixed facet.
#'   Composite weights define linear score combinations. They do not implement
#'   score-dependent panel trimming or an entire competition scoring rule.
#'   For a Person-by-Task G-study, `E` combines Person-by-Task interaction and
#'   within-cell error: relative-error covariance is `E/n_t`, absolute-error
#'   covariance is `(T + E)/n_t`, and universe-score covariance remains `P`.
#'   The same formulas apply to a Person-by-Rater design or to facets selected
#'   by `facets`, with names replaced in the declared order. For a single
#'   facet F and its planned count n, relative error is `E/n` and absolute
#'   error is `(F + E)/n`. A D-study cannot introduce an absent facet.
#'   For fixed tasks represented by score columns in a Person-by-Rater model,
#'   vary only `Raters`. The task set remains fixed, and each planned rater
#'   scores all tasks. See "Fixed tasks as score components" in
#'   [mfrm_g_study()] for the target and allocation requirements.
#'
#'   For raters nested within tasks, use the five components `P`, `T`,
#'   `R(T)`, `PT`, and `E`, where `E` includes the Person-by-Rater-within-Task
#'   interaction. With `n_r` raters per task and `n_t` tasks, relative error is
#'   `PT/n_t + E/(n_r*n_t)`; absolute error adds `T/n_t + R(T)/(n_r*n_t)`.
#'   Universe covariance remains `P`. Other nested facet labels follow the
#'   same rules. All persons receive the same tasks and their nested raters;
#'   two raters per task at six tasks means twelve distinct raters, not two
#'   raters each judging six tasks. Future unequal child counts, changes of
#'   nesting are not supported. Nested G/Phi defaults and their differences
#'   require unchanged source child counts. Nested MLS SEM differences need
#'   a common future child count across plans, which may differ from the source.
#'   General nested coefficient intervals require explicit experimental methods;
#'   general difference enclosures use `method = "mean_square_region"`.
#'   Individual nested SEM intervals allow both counts to change.
#'   See `confint()` above for the derivations and limitations.
#'
#'   For incomplete source data, counts of distinct observed levels are not
#'   per-person replication counts. This function projects a future complete
#'   design; it does not estimate the reliability of the observed sparse
#'   roster, heterogeneous person-specific assignments, or averages of them.
#'
#'   Each score uses diagonal entries; a composite uses `w' Sigma w` for each
#'   of the three covariance matrices. Thus between-score covariance affects
#'   composite dependability. `G = UniverseVariance / (UniverseVariance +
#'   RelativeErrorVariance)` and `Phi` uses absolute error instead. SEMs are
#'   square roots of the corresponding error variances, in the units of the
#'   score or specified composite of mean scores. Multiplying all weights by
#'   a nonzero constant multiplies SEMs by its absolute value but leaves G/Phi
#'   unchanged. Weights express a substantive choice; this function neither
#'   chooses them nor identifies an optimal design. Interpret a difference
#'   only when subtraction is meaningful on the supplied score scales.
#'   Comparing several composites holds the estimated components and planned
#'   designs constant, but different weights may change the construct being
#'   assessed. A higher coefficient alone does not justify that change.
#'   Compare SEMs across composites only when their scales are comparable.
#'
#'   The same weights define the universe-score composite and its observed
#'   mean-score estimate (the equal WWTS/AWTS case in mGENOVA). Different
#'   target and estimation weights, profile reliability, and accuracy at a
#'   cut score are not provided. Unlike mGENOVA's default D-study procedure,
#'   negative variance estimates are not replaced with zero by default. mGENOVA replaces
#'   negative variance components with zero for the D-study unless its
#'   `NEGATIVE` option is selected (manual, pages 15 and 21). That rule concerns
#'   variances, not all negative covariances, and does not guarantee PSD
#'   component matrices.
#'   Do not transfer that option's meaning to GENOVA 3.1 or jGENOVA: their
#'   `NEGATIVE` option controls printing of negative estimates. Their D-studies zero negative
#'   variances with either `ALGORITHM` or its default `EMS` convention; `EMS`
#'   also recalculates other components using the zero replacements. Matching
#'   a raw G-study estimate therefore need not give matching default D-studies.
#'   A matched negative-interaction example gives Phi = 0.73096 with the raw
#'   components retained here, 0.72727 in GENOVA 3.1 `ALGORITHM`, and 0.73532
#'   with its default `EMS`. This is a component-treatment difference, not a
#'   coefficient-formula discrepancy. Select `negative_components` explicitly
#'   to reproduce either GENOVA D-study convention in supported univariate
#'   designs. EMS recalculation uses the source counts, before division by the
#'   future-plan counts. It is not REML or a positive-semidefinite matrix repair.
#'   Multivariate and MINQUE sources require `"raw"`; negative cross-score
#'   covariance is not a negative variance. GENOVA's `NEGATIVE` printing option
#'   is separate from this component-treatment choice.
#'
#'   `gstudy` retains raw estimates and diagnostics. `effective_components`
#'   contains the matrices actually used in projections; their diagnostics are
#'   in `effective_component_diagnostics`. `ComponentTreatment` identifies the
#'   choice in coefficient tables. `ComponentPSD` still describes the original
#'   G-study, so zero replacement does not erase evidence of problematic raw
#'   components. Save the complete result to preserve both sets.
#'   For GENOVA-adjusted results, `confint()` retains population-parameter
#'   intervals based on the original source study alongside adjusted points.
#'   `IntervalEstimate` records the original projection; bounds are not
#'   recentered on the adjusted `Estimate`. Any `IntervalSE` and other sampling
#'   diagnostics retain the original interval-method basis. Summary and plot
#'   methods accept intervals with the same recorded treatment. Supported
#'   population plan differences use `confint(d, reference = ...)`; the separate
#'   delta-estimator comparison API still requires raw components. These
#'   intervals do not characterize the adjusted estimator's sampling law.
#'   See [confint.mfrm_multivariate_d_study()] for assumptions and distinctions.
#'
#'   Original covariance matrices and the selected projected variances remain available. G
#'   requires nonnegative universe and relative-error variances and a positive
#'   sum; Phi uses absolute-error variance instead. Each SEM requires only its
#'   corresponding nonnegative error variance. A failed requirement makes
#'   only that metric `NA`. Zero universe variance yields a coefficient of
#'   zero when its error variance is positive; a zero total variance leaves
#'   the coefficient undefined. Under `"raw"`, no estimates are clipped or replaced.
#'
#'   `Status` summarizes a row as `"Available"`, `"Partially available"`, or
#'   `"Unavailable"`; the four metric-specific status columns give the reasons.
#'   `ComponentPSD` is TRUE only when all G-study component matrices pass
#'   their numerical PSD check. A FALSE value flags components for review
#'   without suppressing otherwise calculable projections. Print and plot
#'   methods display a note for non-PSD components. Calculability does not
#'   validate the joint covariance model or establish precise estimation;
#'   inspect `component_diagnostics` before interpreting estimates. Raw estimates
#'   can yield `Phi > G` when projected absolute error is less than relative
#'   error. Under the stated model, the extra absolute-error contributions
#'   are nonnegative; such a reversal calls for review of the estimated
#'   components, not a conclusion that absolute decisions are more dependable.
#'   Passing the matrix checks does not establish model fit. Outputs are
#'   point projections conditional on the estimated components, without
#'   confidence intervals, missing-data correction, or MFRM latent inference.
#'
#' @return An `mfrm_multivariate_d_study` list containing `coefficients`
#'   (scenario, counts of included facets, score/composite identity,
#'   universe/relative/absolute variances under the selected treatment,
#'   G/Phi, relative/absolute SEM, `ComponentTreatment`, and
#'   row and metric-specific status, and `ComponentPSD`), `covariances`
#'   (three matrices per scenario), `design_grid`,
#'   `weights` (the supplied vector or matrix, reordered to the G-study score
#'   order), `component_diagnostics`,
#'   `negative_components`, `effective_components`,
#'   `effective_component_diagnostics`, and `gstudy` (the unchanged source result).
#'   `summary()` returns `coefficients`.
#'   `print()` returns the same D-study invisibly without computing intervals.
#'   [summary.mfrm_multivariate_d_study()] with `intervals = TRUE` explicitly
#'   calculates and collects supported intervals, methods and source context.
#'   `print()` also shows source level/row counts, future-plan conventions and
#'   the supplied composite weights. Save the full result with `saveRDS(d)`
#'   and reopen with `readRDS()` to preserve that context. A CSV of
#'   `summary(d)` contains the coefficient table only; accompany it with the
#'   weights, source design and row usage, future design and score convention.
#'   Main-effects inputs instead retain their `mfrm_d_study` data-frame result;
#'   see [mfrm_d_study_main_effects] for its fields and limitations.
#' @references Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'   Chapters 9--11. \doi{10.1007/978-1-4757-3456-0}.
#'
#'   Brennan, R. L. (2001). *Manual for mGENOVA, Version 2.1*.
#'   Iowa Testing Programs Occasional Papers, No. 50. Pages 16 and 20--22;
#'   Appendices E and F, pages 74--81.
#'   Crick, J. E., and Brennan, R. L. (2021). *Manual for jGENOVA, Version 1.0*.
#'   Pages 2-10--2-11, 3-5, and A-22--A-23.
#'   Both manuals are available from the author's
#'   [GENOVA Suite](https://brennancrickgenova.org/genova-suite/).
#' @seealso [plot.mfrm_multivariate_d_study()], [mfrm_g_study()],
#'   [mfrm_multivariate_d_compare()], [mfrm_d_study()]
#' @examples
#' # Published synthetic two-score data: the same six tasks for ten persons.
#' tasks <- read.csv(system.file("extdata", "mgenova-table12.csv", package = "mfrmr"))
#' g <- mfrm_g_study(tasks, c("V", "W"), rater = NULL)
#'
#' # Question: would doubling tasks improve the dependability of W minus V?
#' d <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)),
#'   weights = c(V = -1, W = 1))
#' summary(d)
#' plot(d)
#' # G rises from 0.300 to 0.462, conditional on these estimated components.
#' # This is a planning projection, not evidence from twelve observed tasks.
#'
#' # Optional alternative views, GENOVA conventions and composite weights.
#' # example("mfrm_d_study", run.dontrun = TRUE)
#' \dontrun{
#' plot(d, type = "sem")
#' # Without weights, report the original scores without creating a composite.
#' individual <- mfrm_d_study(g, data.frame(Tasks = c(6, 12)))
#' plot(individual, score = "V")
#'
#' # Explicit GENOVA convention for a single observed score.
#' g_one <- mfrm_g_study(tasks, "V", rater = NULL)
#' d_ems <- mfrm_d_study(g_one, data.frame(Tasks = c(6, 12)),
#'   negative_components = "genova_ems")
#' summary(d_ems)
#' ci_ems <- confint(d_ems, assumption = "normal")
#' # Adjusted point and original-source population interval have separate bases.
#' ci_ems[, c("Estimate", "IntervalEstimate", "Lower", "Upper", "IntervalBasis")]
#' d_ems$effective_components # Compare with d_ems$gstudy$components.
#' plot(d_ems, metric = "G", intervals = ci_ems) # Original-source population bounds.
#'
#' # Compare weight choices on the same task-count scenarios.
#' weights <- cbind(Equal = c(V = 0.5, W = 0.5),
#'   W_focus = c(V = 0.2, W = 0.8), Difference = c(V = -1, W = 1))
#' alternatives <- mfrm_d_study(g,
#'   data.frame(Tasks = c(6, 12)), weights = weights)
#' summary(alternatives) # Original scores plus all three named composites.
#' plot(alternatives, composite = "Equal")
#' plot(alternatives, composite = "Difference", type = "sem")
#' # Choose weights for the intended interpretation, not solely for higher G.
#' }
#' @name mfrm_d_study
#' @aliases mfrm_multivariate_d_study
#' @method mfrm_d_study mfrm_multivariate_gstudy
#' @export
mfrm_d_study.mfrm_multivariate_gstudy <- function(x, design_grid = NULL, weights = NULL,
    negative_components = c("raw", "genova_algorithm", "genova_ems"), ...) {
  rlang::check_dots_empty()
  negative_components <- match.arg(negative_components)
  if (!inherits(x, "mfrm_multivariate_gstudy") ||
      !(identical(x$design$calculation_version, 1L) || identical(x$design$calculation_version, 2L) ||
        identical(x$design$calculation_version, 3L))) {
    stop("`x` must be a current mfrm_multivariate_gstudy result.", call. = FALSE)
  }
  scores <- x$design$scores
  counts <- x$design$counts
  facet_names <- names(counts)[-1L]
  count_columns <- x$design$count_columns
  # Previously saved task/rater results did not store an explicit count map.
  if (is.null(count_columns) && all(facet_names %in% c("Rater", "Task"))) {
    count_columns <- c(Rater = "Raters", Task = "Tasks")[facet_names]
  }
  if (!is.numeric(counts) || !identical(names(counts)[1L], "Person") ||
      !.mfrm_mvgt_valid_facets(facet_names) || any(!is.finite(counts)) ||
      any(counts < 2 | counts != floor(counts)) ||
      !identical(names(count_columns), facet_names) ||
      !.mfrm_mvgt_valid_facets(unname(count_columns))) {
    stop("The G-study factor or count identities are incomplete or altered; recreate the result.", call. = FALSE)
  }
  one_facet <- length(facet_names) == 1L
  nesting <- .mfrm_mvgt_nesting(x$design$nesting, facet_names)
  if (identical(x$design$calculation_version, 3L) != !is.null(nesting)) {
    stop("The G-study nesting specification is incomplete or altered; recreate the result.", call. = FALSE)
  }
  if (!is.null(nesting)) {
    child_counts <- x$design$child_counts
    if (!is.numeric(child_counts) || length(child_counts) != counts[[unname(nesting)]] ||
        any(!is.finite(child_counts)) || any(child_counts < 1 | child_counts != floor(child_counts)) ||
        sum(child_counts) != counts[[names(nesting)]]) {
      stop("The G-study child counts are incomplete or altered; recreate the result.", call. = FALSE)
    }
  }
  sources <- .mfrm_mvgt_sources(names(counts), nesting)
  if (!is.character(scores) || !length(scores) || anyNA(scores) || anyDuplicated(scores) ||
      !identical(names(x$components), sources) ||
      !is.numeric(x$score_scale) || !identical(names(x$score_scale), scores) ||
      any(!is.finite(x$score_scale)) || any(x$score_scale < 0) ||
      !all(vapply(x$components, function(a) is.matrix(a) && is.numeric(a) &&
        !is.complex(a) && identical(dim(a), rep(length(scores), 2L)) &&
        identical(dimnames(a), list(scores, scores)) && all(is.finite(a)) &&
        identical(a, t(a)), logical(1)))) {
    stop("The G-study covariance matrices or score identities are incomplete or altered; recreate the result.", call. = FALSE)
  }
  diagnostics <- .mfrm_mvgt_diagnostics(x$components, x$score_scale)
  grid_names <- unname(count_columns)
  if (is.null(design_grid)) {
    if (identical(x$design$complete, FALSE) ||
        (!is.null(nesting) && length(unique(child_counts)) != 1L)) {
      stop("An incomplete or unequal G-study requires an explicit `design_grid` for a future complete balanced design; observed pool counts are not planned replication.", call. = FALSE)
    }
    future_counts <- setNames(counts[-1L], grid_names)
    if (!is.null(nesting)) future_counts[count_columns[[names(nesting)]]] <- child_counts[1L]
    design_grid <- as.data.frame(as.list(future_counts), check.names = FALSE)
  }
  design_grid <- mfrm_d_study_count_columns(design_grid, facet_names, grid_names)
  if (!is.data.frame(design_grid) || !nrow(design_grid) ||
      !setequal(names(design_grid), grid_names) || ncol(design_grid) != length(grid_names) ||
      !all(vapply(design_grid, function(n) is.numeric(n) && !is.object(n) &&
        is.null(dim(n)) && !is.complex(n) && all(is.finite(n)) &&
        all(n >= 1 & n <= .Machine$integer.max & n == floor(n)), logical(1)))) {
    stop(paste0("`design_grid` must contain only positive integer ",
      paste(grid_names, collapse = " and "), " columns, one scenario per row."), call. = FALSE)
  }
  design_grid <- as.data.frame(design_grid[grid_names])
  rownames(design_grid) <- NULL
  vectors <- diag(length(scores))
  colnames(vectors) <- scores
  kind <- rep("Score", length(scores))
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (!is.numeric(weights) || is.object(weights) || is.complex(weights) ||
          nrow(weights) != length(scores) || !ncol(weights) ||
          is.null(rownames(weights)) || anyNA(rownames(weights)) ||
          anyDuplicated(rownames(weights)) || !setequal(rownames(weights), scores) ||
          is.null(colnames(weights)) || anyNA(colnames(weights)) ||
          anyDuplicated(colnames(weights)) || any(!nzchar(trimws(colnames(weights)))) ||
          any(!is.finite(weights)) || any(colSums(weights != 0) == 0)) {
        stop("A `weights` matrix must name every score once in its rows and have unique, nonempty composite names in its columns, with finite values and at least one nonzero weight per column.", call. = FALSE)
      }
      weights <- weights[scores, , drop = FALSE]
      vectors <- cbind(vectors, weights)
      kind <- c(kind, rep("Composite", ncol(weights)))
    } else {
      if (!is.numeric(weights) || is.object(weights) || is.complex(weights) ||
          !is.null(dim(weights)) || length(weights) != length(scores) ||
          is.null(names(weights)) || anyNA(names(weights)) || anyDuplicated(names(weights)) ||
          !setequal(names(weights), scores) || any(!is.finite(weights)) ||
          !any(weights != 0)) {
        stop("`weights` must name every score once with finite values and at least one nonzero value.", call. = FALSE)
      }
      weights <- weights[scores]
      vectors <- cbind(vectors, Composite = weights)
      kind <- c(kind, "Composite")
    }
  }
  projection_scale <- apply(abs(vectors), 2L, max)
  if (any(!is.finite(projection_scale^2)) || any(projection_scale^2 == 0)) {
    stop("Composite variances exceed numeric range; rescale the weights.", call. = FALSE)
  }
  vectors <- sweep(vectors, 2L, projection_scale, "/")
  c <- .mfrm_mvds_components(x, negative_components)
  effective_diagnostics <- .mfrm_mvgt_diagnostics(c, x$score_scale)
  covariances <- vector("list", nrow(design_grid))
  tables <- vector("list", nrow(design_grid))
  for (i in seq_len(nrow(design_grid))) {
    # Count columns may be integer storage, but their products need doubles.
    n1 <- as.double(design_grid[[1L]][i])
    if (!is.null(nesting)) {
      n_child <- as.double(design_grid[[count_columns[[names(nesting)]]]][i])
      n_parent <- as.double(design_grid[[count_columns[[unname(nesting)]]]][i])
      relative <- c[[4L]] / n_parent + c[[5L]] / (n_child * n_parent)
      absolute <- relative + c[[2L]] / n_parent + c[[3L]] / (n_child * n_parent)
    } else if (one_facet) {
      relative <- c[[3L]] / n1
      absolute <- relative + c[[2L]] / n1
    } else {
      n2 <- as.double(design_grid[[2L]][i])
      relative <- c[[4L]] / n1 + c[[5L]] / n2 + c[[7L]] / (n1 * n2)
      absolute <- relative + c[[2L]] / n1 + c[[3L]] / n2 + c[[6L]] / (n1 * n2)
    }
    covariances[[i]] <- list(Universe = c$Person, RelativeError = relative, AbsoluteError = absolute)
    unit_variances <- matrix(vapply(covariances[[i]],
      function(a) colSums(vectors * (a %*% vectors)), numeric(ncol(vectors))),
      nrow = ncol(vectors), ncol = 3L)
    variances <- unit_variances * projection_scale^2
    if (any(!is.finite(variances)) || any(variances == 0 & unit_variances != 0)) {
      stop("D-study projection exceeded numeric range; rescale scores or weights.", call. = FALSE)
    }
    u <- unit_variances[, 1L]
    r <- unit_variances[, 2L]
    a <- unit_variances[, 3L]
    metrics <- c("G", "Phi", "RelativeSEM", "AbsoluteSEM")
    values <- matrix(NA_real_, length(u), 4L, dimnames = list(NULL, metrics))
    metric_status <- matrix("Available", length(u), 4L,
      dimnames = list(NULL, paste0(metrics, "Status")))
    for (j in 1:2) {
      error <- if (j == 1L) r else a
      metric_status[error < 0, c(j, j + 2L)] <- "Negative projected error variance"
      metric_status[u < 0, j] <- "Negative universe variance"
      metric_status[u == 0 & error == 0, j] <- "Zero total variance"
      ok <- metric_status[, j] == "Available"
      # Scale each denominator before addition to avoid overflow.
      scale <- pmax(u[ok], error[ok])
      values[ok, j] <- (u[ok] / scale) / (u[ok] / scale + error[ok] / scale)
      ok_sem <- metric_status[, j + 2L] == "Available"
      values[ok_sem, j + 2L] <- sqrt(error[ok_sem]) * projection_scale[ok_sem]
    }
    available_count <- rowSums(metric_status == "Available")
    status <- ifelse(available_count == 4L, "Available",
      ifelse(available_count == 0L, "Unavailable", "Partially available"))
    tables[[i]] <- data.frame(Scenario = i, design_grid[i, , drop = FALSE],
      Kind = kind, Score = colnames(vectors), UniverseVariance = variances[, 1L],
      RelativeErrorVariance = variances[, 2L], AbsoluteErrorVariance = variances[, 3L],
      values, Status = status, metric_status,
      ComponentPSD = all(diagnostics$PositiveSemidefinite),
      ComponentTreatment = negative_components, row.names = NULL, check.names = FALSE)
  }
  structure(list(coefficients = do.call(rbind, tables), covariances = covariances,
    design_grid = design_grid, weights = weights, component_diagnostics = diagnostics,
    negative_components = negative_components, effective_components = c,
    effective_component_diagnostics = effective_diagnostics,
    gstudy = x), class = "mfrm_multivariate_d_study")
}

#' @rdname mfrm_g_study
#' @param x,object A result from the corresponding G-study or D-study function.
#' @param digits Display precision passed to the printed tables; `NULL` uses
#'   R's current setting. Stored estimates are unchanged.
#' @param row.names Logical; print table row names? Defaults to `FALSE`.
#' @param ... Further table-printing arguments for `print()`; unused arguments
#'   to `summary()` are rejected.
#' @method print mfrm_multivariate_gstudy
#' @export
print.mfrm_multivariate_gstudy <- function(x, digits = NULL, row.names = FALSE, ...) {
  cat("Multivariate observed-score G-study\n")
  label <- if (identical(x$design$complete, FALSE)) "Incomplete observed design:" else
    if (!is.null(x$design$nesting)) "Complete nested design:" else "Complete crossed design:"
  counts <- x$design$counts
  units <- c(Person = "persons", Rater = "raters", Task = "tasks")[names(counts)]
  units[is.na(units)] <- paste(names(counts)[is.na(units)], "levels")
  cat(label, paste(counts, units, collapse = ", "), "\n")
  if (!is.null(x$design$nesting)) {
    cat(names(x$design$nesting), "identities are local to", unname(x$design$nesting),
      "; the observed child count is the total across parents.\n")
    cat("Children per parent:", paste(names(x$design$child_counts), x$design$child_counts,
      sep = " = ", collapse = ", "), "\n")
  }
  cat(x$design$method, "; highest-order interaction and residual are combined.\n", sep = "")
  if (!is.null(x$data_usage)) {
    cat("Rows:", x$data_usage$counts[["UsedRows"]], "used of", x$data_usage$counts[["InputRows"]],
      ";", x$data_usage$counts[["ExcludedRows"]], "explicitly omitted.\n")
  }
  if (!is.null(x$estimation)) {
    cat("Covariance-component rank:", x$estimation$rank,
      "; scaled moment-system condition number:", format(x$estimation$condition_number, digits = 4), "\n")
  }
  print(x$component_diagnostics, digits = digits, row.names = row.names, ...)
  cat("Raw covariance components are retained; this table shows point estimates.\n")
  if (identical(x$design$method, "Balanced multivariate ANOVA"))
    cat("Use confint(x, assumption = 'normal') for eligible score-component SEs and intervals.\n")
  else cat("Use vcov(x, assumption = 'normal') for component sampling SEs; confint(x, method = 'reml_profile', assumption = 'normal') provides explicit experimental intervals.\n")
  invisible(x)
}

#' @rdname mfrm_g_study
#' @method summary mfrm_multivariate_gstudy
#' @export
summary.mfrm_multivariate_gstudy <- function(object, ...) {
  rlang::check_dots_empty()
  object$component_diagnostics
}

#' @rdname mfrm_d_study
#' @inheritParams print.mfrm_multivariate_gstudy
#' @param ... Arguments for the selected method. The fitted-MFRM method accepts
#'   arguments to [mfrm_generalizability()]; already estimated G-studies reject
#'   unused arguments. Printing accepts arguments for its coefficient table.
#' @method print mfrm_multivariate_d_study
#' @export
print.mfrm_multivariate_d_study <- function(x, digits = NULL, row.names = FALSE, ...) {
  if (!identical(.mfrm_mvds_policy(x), "raw")) .mfrm_mvds_replay(x)
  cat("Multivariate observed-score D-study\n")
  .mfrm_mvgt_print_source(x$gstudy$design, x$gstudy$data_usage$counts)
  cat("Future complete balanced plans; counts apply to every person and score.\n")
  .mfrm_mvgt_print_future_counts(x$gstudy$design)
  convention <- x$gstudy$design$score_convention
  if (is.null(convention)) convention <- paste("Means over common random",
    paste(tolower(names(x$design_grid)), collapse = " and "))
  cat(paste0(convention, "; weights are used as supplied.\n"))
  cat(.mfrm_mvds_component_note(x), "\n")
  if (!is.null(x$weights)) {
    cat("Composite weights (not normalized):\n")
    print(if (is.matrix(x$weights)) x$weights else
      matrix(x$weights, ncol = 1L, dimnames = list(names(x$weights), "Composite")), digits = digits)
  }
  if (any(!x$component_diagnostics$PositiveSemidefinite)) {
    cat("Warning: non-PSD covariance components in the original G-study; inspect component_diagnostics and the selected treatment.\n")
  }
  print(x$coefficients, digits = digits, row.names = row.names, ...)
  cat("Point projections conditional on estimated covariance components; no confidence intervals.\n")
  invisible(x)
}

.mfrm_mvgt_print_source <- function(design, rows) {
  if (is.null(design)) {
    cat("Source design not recorded in this saved comparison; consult the original G-study.\n")
    return(invisible(NULL))
  }
  cat("Source G-study:", design$method, "\n")
  cat("Observed level counts (not future replication):",
    paste(names(design$counts), design$counts, sep = " = ", collapse = ", "), "\n")
  if (!is.null(design$nesting)) cat("The observed child count is the total across parents.\n")
  if (!is.null(rows)) {
    cat("Source rows:", rows[["UsedRows"]], "used of", rows[["InputRows"]],
      ";", rows[["ExcludedRows"]], "explicitly omitted.\n")
  } else if (!is.null(design$rows)) cat("Source rows analyzed:", design$rows, "\n")
  if (identical(design$complete, FALSE) || identical(design$balanced, FALSE)) {
    cat("The observed source is incomplete or unequal; future counts do not reproduce its observed assignments.\n")
  }
  invisible(NULL)
}

.mfrm_mvgt_print_future_counts <- function(design) {
  if (!is.null(design$nesting)) {
    child <- names(design$nesting)
    cat("Future", design$count_columns[[child]], "counts levels per",
      unname(design$nesting), "; not the total across parents.\n")
  }
  invisible(NULL)
}
