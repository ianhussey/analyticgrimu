#' GRIM-U Forensic Consistency Check
#'
#' Checks the mathematical consistency of reported Mann-Whitney U test results.
#' Performs a "Triangulation" check across Sample Sizes (N), Test Statistic (U), and P-value (P).
#'
#' 1. N <-> U: Is the reported U physically possible? (Bounds & Granularity).
#' 2. N <-> P: Is the reported P mathematically possible? (Grid search of attainable p-values).
#' 3. U <-> P: Does the specific reported U generate the reported P?
#'
#' The candidate p-value set checked at each U covers:
#' * Exact p-value, R / SciPy doubled-tail convention (`p_exact`).
#' * Exact p-value, SPSS Exact / StatXact deviation-from-mean convention (`p_exact_dev`).
#' * Mid-p (`p_mid`).
#' * Asymptotic p-value (continuity-corrected and uncorrected) under every achievable
#'   tie-correction K = sum(t_i^3 - t_i) for tie-pattern partitions of size <= N.
#'   This covers no-ties (K = 0) and arbitrary multi-tie patterns, matching the
#'   tie-corrected variance used by R `wilcox.test()`, SPSS, Stata, and SciPy.
#'
#' For one-sided alternatives the search is run for *both* directions and the results
#' are OR'd together. This is because which group is labelled "1" determines whether
#' the reported one-sided p is the left- or right-tail value, and this labelling
#' convention differs between R `wilcox.test()`, SPSS, and SciPy. Running both
#' directions prevents a correctly-reported result from being rejected purely on
#' opposite-labelling grounds. The supplied `alternative` argument is still used to
#' shape the asymptotic deviation and continuity-correction signs.
#'
#' **Caveat: exact permutation tests with ties.** Tools such as
#' `coin::wilcox_test(distribution = "exact")`, SAS `PROC NPAR1WAY EXACT`, StatXact,
#' and the SPSS Exact Tests Add-on compute the exact permutation p-value
#' *conditional on the observed tie pattern*. That p-value is not reachable by any
#' enumeration indexed only by sample sizes, because the partition is not recoverable
#' from `(n1, n2)`. An `inconsistent` verdict against a report from one of those tools
#' may therefore reflect this method's coverage gap rather than a real inconsistency
#' in the report. Users investigating a flagged result should check the source
#' software before treating the flag as evidence of an error.
#'
#' @param n1 Integer. Sample size of group 1.
#' @param n2 Integer. Sample size of group 2.
#' @param u_reported Numeric, optional. The reported U statistic. If NA, checks P-value consistency only.
#' @param p_reported Numeric. The reported p-value.
#' @param comparison Character. "equal" (default) or "less_than" (for reported inequalities like "p < .05").
#' @param digits Integer. The number of decimal places the p-value was reported to. Mandatory if p_reported is supplied.
#' @param rounding Character vector, optional. Allowed rounding methods: "round", "trunc", "up".
#'   If NULL, defaults to `c("round")`. The `"round"` window is symmetric (`+/- 0.5 * 10^-digits`)
#'   around `p_reported` and so covers both round-half-up and banker's rounding for any
#'   non-midpoint value (the two algorithms differ only at exact midpoints, which have
#'   measure zero for continuous p-values). `"trunc"` (floor to precision) and `"up"`
#'   (ceiling to precision) are available as opt-in options.
#' @param p_min Numeric, optional. Manually override the lower bound of the p-value search window.
#' @param p_max Numeric, optional. Manually override the upper bound of the p-value search window.
#' @param alternative Character. One of "two.sided", "less", or "greater". Defaults to "two.sided".
#'   For one-sided alternatives, both directions are enumerated and the consistency check is OR'd
#'   across them (see Description).
#'
#' @return A list containing:
#'   \item{summary}{A single-row tibble with the overall consistency verdict (`consistent`), component checks
#'     (`u_bounds_consistent`, `u_matches_p`, `p_granularity_consistent`), and diagnostic flags.}
#'   \item{details}{A tibble of all candidate U values found within the p-value search window,
#'     showing which specific formula(s) matched the reported p-value.}
#' @export
grimu_check <- function(n1, n2,
                        u_reported = NA_real_,
                        comparison = "equal",
                        p_reported,
                        digits,
                        rounding = NULL,
                        p_min = NULL, p_max = NULL,
                        alternative = "two.sided") {

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

  # --- Fail-Fast Validation ---
  if (is.na(n1) || is.na(n2)) stop("n1 and n2 must be supplied")

  is_whole <- function(x) !is.na(x) && abs(x - round(x)) < 1e-10
  if (!is_whole(n1) || !is_whole(n2)) {
    return(list(summary = tibble(n1=n1, n2=n2, consistent=NA, error="Non-integer N"), details=tibble()))
  }

  if (!is.na(p_reported)) {
    if (is.null(digits) || is.na(digits) || !is_whole(digits)) {
      stop("If p_reported is supplied, digits must be supplied as an integer.")
    }
  }

  if (!is.na(p_reported) & (p_reported < 0 || p_reported > 1)) stop("p_reported must be [0,1]")
  if (!is.numeric(u_reported)) stop("u_reported must be numeric or NA")

  # --- Step 1: Check Global U Validity ---
  if (!is.na(u_reported)) {
    u_res <- check_u_consistency(n1, n2, u_reported)
  } else {
    u_res <- tibble::tibble(
      u_bounds_consistent = NA,
      u_granularity_consistent = NA,
      u_possible = NA
    )
  }

  # --- Step 2: Range Detection for P ---
  if (is.null(p_min) || is.null(p_max)) {
    # Default window for bounds search (wide enough to catch all rounding types)
    window <- 10^(-digits) * 5
    p_max_search <- p_reported + window
    # Calculate raw lower bound
    raw_p_min <- p_reported - window

    # If comparison is inequality, or window touches 0, anchor the "deep tail".
    if (comparison == "less_than" || raw_p_min <= 0) {
      p_min_search <- NA_real_
    } else {
      # Standard case: Window is strictly positive (e.g., 0.04 to 0.06)
      p_min_search <- raw_p_min
    }
  } else {
    p_min_search <- p_min
    p_max_search <- p_max
    # If user manually provides 0, treat it as NA for Z-score purposes
    if (p_min_search <= 0) p_min_search <- NA_real_
  }

  # --- Step 3: Bounds Calculation ---
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  sigma_est <- sqrt((n1 * n2 * (N + 1)) / 12)
  max_u <- n1 * n2

  # Helper: Convert P to Z (Signed)
  p_to_z <- function(p, alt) {
    p_safe <- min(1, max(0, p))
    if (is.na(p)) return(Inf)
    if (alt == "two.sided") return(qnorm(1 - p_safe / 2))
    if (alt == "less")      return(qnorm(p_safe))
    if (alt == "greater")   return(qnorm(1 - p_safe))
  }

  # Helper: enumerate candidate U range under a given alternative and pull p-values.
  # Returns the candidate tibble with attributes for the U-bounds used.
  search_candidates <- function(alt) {
    z_deep <- p_to_z(if (is.na(p_min_search)) 0 else p_min_search, alt)
    z_shallow <- p_to_z(p_max_search, alt)

    u_bound_1 <- mu + z_deep * sigma_est
    u_bound_2 <- mu + z_shallow * sigma_est

    if (is.na(u_bound_1) || is.nan(u_bound_1)) u_bound_1 <- if (alt == "greater") max_u else 0
    if (is.na(u_bound_2) || is.nan(u_bound_2)) u_bound_2 <- if (alt == "greater") mu else mu

    raw_start <- min(u_bound_1, u_bound_2)
    raw_end   <- max(u_bound_1, u_bound_2)

    u_start_est <- max(0, floor(raw_start - 2))
    u_end_est   <- min(max_u, ceiling(raw_end + 2))

    candidates <- grimu_map_pvalues(n1, n2,
                                    u_min = u_start_est, u_max = u_end_est,
                                    alternative = alt)
    candidates$alt_used <- alt
    attr(candidates, "u_range") <- c(u_start_est, u_end_est)
    candidates
  }

  # WHY: Which group is labelled "1" determines whether the right-tail or
  # left-tail p is reported under one-sided tests, and this convention
  # differs between R wilcox.test(), SPSS, and SciPy. A reporter who used
  # the opposite labelling submits the complementary one-sided p; without
  # this swap we would reject correct reports purely on labelling grounds.
  if (alternative == "two.sided") {
    primary <- search_candidates("two.sided")
    results_df <- primary
    u_range_min <- attr(primary, "u_range")[1]
    u_range_max <- attr(primary, "u_range")[2]
  } else {
    opposite <- if (alternative == "less") "greater" else "less"
    primary <- search_candidates(alternative)
    swapped <- search_candidates(opposite)
    results_df <- bind_rows(primary, swapped)
    u_range_min <- min(attr(primary, "u_range")[1], attr(swapped, "u_range")[1])
    u_range_max <- max(attr(primary, "u_range")[2], attr(swapped, "u_range")[2])
  }

  # --- Step 4: Check P-Value Consistency ---
  if (is.null(rounding)) {
    methods <- c("round")
  } else {
    methods <- match.arg(rounding, c("round", "trunc", "up"), several.ok = TRUE)
  }

  epsilon <- 10^(-digits)
  buffer <- 1e-14

  # Vectorised match check across a column of candidate p-values.
  is_match_vec <- function(val) {
    matched <- rep(FALSE, length(val))
    na_mask <- is.na(val)
    if ("round" %in% methods) {
      matched <- matched | (!na_mask &
                              val >= (p_reported - 0.5 * epsilon - buffer) &
                              val <  (p_reported + 0.5 * epsilon - buffer))
    }
    if ("trunc" %in% methods) {
      matched <- matched | (!na_mask &
                              val >= (p_reported - buffer) &
                              val <  (p_reported + epsilon - buffer))
    }
    if ("up" %in% methods) {
      matched <- matched | (!na_mask &
                              val >  (p_reported - epsilon + buffer) &
                              val <= (p_reported + buffer))
    }
    matched
  }

  check_col_vec <- function(col_val) {
    if (comparison == "equal") {
      is_match_vec(col_val)
    } else {
      !is.na(col_val) & col_val < p_reported
    }
  }

  # Identify the K-indexed asymptotic columns to distinguish K = 0 (no ties)
  # from K > 0 (ties) in the diagnostic flags.
  corr_K_cols   <- grep("^p_corr_K",   names(results_df), value = TRUE)
  uncorr_K_cols <- grep("^p_uncorr_K", names(results_df), value = TRUE)
  corr_ties_cols   <- setdiff(corr_K_cols,   "p_corr_K0")
  uncorr_ties_cols <- setdiff(uncorr_K_cols, "p_uncorr_K0")

  # Helper to OR check_col_vec across a set of columns (returns a row-vector).
  any_col_matches <- function(df, cols) {
    if (length(cols) == 0) return(rep(FALSE, nrow(df)))
    Reduce(`|`, lapply(cols, function(cn) check_col_vec(df[[cn]])))
  }

  results_checked <- results_df %>%
    mutate(
      valid_exact          = check_col_vec(p_exact),
      valid_exact_dev      = check_col_vec(p_exact_dev),
      valid_mid            = check_col_vec(p_mid),
      valid_corr_no_ties   = check_col_vec(p_corr_K0),
      valid_uncorr_no_ties = check_col_vec(p_uncorr_K0)
    ) %>%
    mutate(
      valid_corr_ties   = any_col_matches(., corr_ties_cols),
      valid_uncorr_ties = any_col_matches(., uncorr_ties_cols)
    ) %>%
    mutate(
      is_consistent = valid_exact | valid_exact_dev | valid_mid |
        valid_corr_no_ties | valid_uncorr_no_ties |
        valid_corr_ties | valid_uncorr_ties
    ) %>%
    # Filter for display (diagnostic mode): keep rows that are consistent OR
    # have any p-value in the general search window.
    filter(
      is_consistent |
        if_any(starts_with("p_"),
               ~ . >= (if (is.na(p_min_search)) 0 else p_min_search) &
                 . <= p_max_search)
    )

  # --- Step 5: Triangulate U and P ---

  # A. Is P possible at all?
  p_consistent <- any(results_checked$is_consistent)

  # B. Is the reported U consistent with the reported P?
  if (!is.na(u_reported)) {
    u_consistent <- results_checked %>%
      filter(is_consistent) %>%
      filter(abs(U - u_reported) < 1e-10) %>%
      nrow() > 0
  } else {
    u_consistent <- NA
  }

  # C. Overall Consistency
  final_consistent <- if (!is.na(u_reported)) {
    isTRUE(u_res$u_possible) && isTRUE(u_consistent)
  } else {
    p_consistent
  }

  summary_df <- tibble(
    n1 = n1,
    n2 = n2,
    u_reported = u_reported,
    p_reported = p_reported,
    alternative = alternative,
    rounding = paste(methods, collapse = "+"),

    # Overall Verdict
    consistent = final_consistent,

    # Bounds
    p_range_min = if (is.na(p_min_search)) 0 else p_min_search,
    p_range_max = p_max_search,
    u_range_min = u_range_min,
    u_range_max = u_range_max,

    # Component Checks
    u_bounds_consistent = u_res$u_bounds_consistent,
    u_granularity_consistent = u_res$u_granularity_consistent,
    u_matches_p = u_consistent,

    # NB no p_bounds_consistent: it's tautological (p min/max are derived from reported p)
    p_granularity_consistent = p_consistent,

    # Method-resolved diagnostic flags: which convention matched, if any.
    p_matches_exact_r    = any(results_checked$valid_exact),
    p_matches_exact_spss = any(results_checked$valid_exact_dev),
    p_matches_mid        = any(results_checked$valid_mid),
    p_matches_no_ties    = any(results_checked$valid_corr_no_ties | results_checked$valid_uncorr_no_ties),
    p_matches_ties       = any(results_checked$valid_corr_ties   | results_checked$valid_uncorr_ties)
  )

  return(list(summary = summary_df, details = results_checked))
}
