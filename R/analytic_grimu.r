library(tidyverse)
library(roundwork)

#' Check Consistency of Mann-Whitney U Statistic with Sample Sizes
#'
#' Verifies if a reported U statistic is mathematically possible given the group sample sizes.
#' Checks two properties:
#' 1. Bounds: U must be between 0 and n1 * n2.
#' 2. Granularity: U must be an integer (no ties) or half-integer (with ties).
#'
#' @param n1 Integer. Sample size of group 1.
#' @param n2 Integer. Sample size of group 2.
#' @param u_reported Numeric. The reported U statistic to check.
#'
#' @return A single-row tibble containing:
#'   \item{u_bounds_consistent}{Logical. TRUE if 0 <= U <= n1*n2.}
#'   \item{u_granularity_consistent}{Logical. TRUE if U is an integer or half-integer.}
#'   \item{u_possible}{Logical. TRUE if both checks pass.}
#' @export
#' check_u_consistency <- function(n1, n2, u_reported) {
  u_max <- n1 * n2
  
  # 1. Check Bounds
  # U cannot be negative or exceed n1*n2
  in_bounds <- dplyr::between(u_reported, 0, u_max)
  
  # 2. Check Granularity (Integer or Half-Integer)
  # Multiplying by 2 should result in an integer (e.g., 20.5 * 2 = 41.0)
  # Using a small float tolerance for safety
  val_doubled <- u_reported * 2
  is_valid_granularity <- abs(val_doubled - round(val_doubled)) < 1e-10
  
  res <- tibble::tibble(
    u_bounds_consistent = in_bounds,
    u_granularity_consistent = is_valid_granularity,
    u_possible = in_bounds & is_valid_granularity
  )
  
  return(res)
}

#' Generate Possible Mann-Whitney U P-values
#'
#' The core engine for GRIM-U. Generates a grid of valid U statistics (integers and half-integers)
#' within a specified range and calculates the corresponding p-values using five different methods:
#' 1. Exact Method (valid for integers only).
#' 2. Asymptotic with Continuity Correction (No Ties variance).
#' 3. Asymptotic uncorrected (No Ties variance).
#' 4. Asymptotic with Continuity Correction (Ties variance).
#' 5. Asymptotic uncorrected (Ties variance).
#'
#' @param n1 Integer. Sample size of group 1.
#' @param n2 Integer. Sample size of group 2.
#' @param u_min Numeric, optional. Lower bound for U search. Defaults to 0 or mean depending on `alternative`.
#'   Will be snapped to the nearest half-integer.
#' @param u_max Numeric, optional. Upper bound for U search. Defaults to mean or max depending on `alternative`.
#'   Will be snapped to the nearest half-integer.
#' @param alternative Character string. One of "two.sided", "less", or "greater". Defaults to "two.sided".
#'
#' @return A tibble where each row is a candidate U value, containing:
#'   \item{U}{The U statistic (integer or half-integer).}
#'   \item{is_integer}{Logical indicating if U is a whole number.}
#'   \item{p_exact}{Exact p-value (NA for non-integers).}
#'   \item{p_corr_no_ties}{Asymptotic p-value (continuity corrected, no ties).}
#'   \item{p_uncorr_no_ties}{Asymptotic p-value (uncorrected, no ties).}
#'   \item{p_corr_tied}{Asymptotic p-value (continuity corrected, ties).}
#'   \item{p_uncorr_tied}{Asymptotic p-value (uncorrected, ties).}
#' @export
grimu_map_pvalues <- function(n1, n2, u_min = NULL, u_max = NULL, alternative = "two.sided") {
  
  # Safety Check: Input Validation
  if (is.na(n1) || is.na(n2)) return(tibble(U = numeric(), is_integer = logical()))
  
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  max_u <- n1 * n2
  
  # Default bounds logic
  if (is.null(u_min) || is.null(u_max)) {
    if (alternative == "greater") {
      if (is.null(u_min)) u_min <- floor(mu) 
      if (is.null(u_max)) u_max <- max_u
    } else {
      if (is.null(u_min)) u_min <- 0
      if (is.null(u_max)) u_max <- ceiling(mu) 
    }
  }
  
  # Enforce Half-Integer Lattice
  # Snap arbitrary bounds to the nearest valid U values (integers or half-integers).
  # If u_min = 0.2 -> snaps up to 0.5
  # If u_max = 0.8 -> snaps down to 0.5
  
  # Handle NULLs before math (though logic above ensures they aren't NULL)
  if (is.null(u_min)) u_min <- 0
  if (is.null(u_max)) u_max <- max_u
  
  u_start <- ceiling(u_min * 2) / 2
  u_end   <- floor(u_max * 2) / 2
  
  # Safety Clamp
  u_start <- max(0, u_start)
  u_end   <- min(max_u, u_end)
  
  # Safety Check: If start > end, return empty tibble immediately
  if (u_start > u_end) return(tibble(U = numeric(), is_integer = logical()))
  
  vals <- roundwork::round_up(seq(u_start, u_end, by = 0.5), 1)
  
  # Constants
  sigma_no_ties <- sqrt((n1 * n2 * (N + 1)) / 12)
  correction_term <- (n1 * n2 * 6) / (12 * N * (N - 1))
  sigma_one_tie <- sqrt((n1 * n2 * (N + 1)) / 12 - correction_term)
  
  # Vectorized Calc
  results_df <- tibble(U = vals) %>%
    mutate(
      is_integer = (U %% 1 == 0),
      
      # Exact
      p_exact = if_else(is_integer, case_when(
        alternative == "less"      ~ pwilcox(U, n1, n2),
        alternative == "greater"   ~ pwilcox(U - 1, n1, n2, lower.tail = FALSE),
        alternative == "two.sided" ~ {
          p_lower <- pwilcox(U, n1, n2)
          p_upper <- pwilcox(U - 1, n1, n2, lower.tail = FALSE)
          2 * pmin(p_lower, p_upper)
        }
      ), NA_real_),
      
      # Deviations
      dev_cc = case_when(
        alternative == "two.sided" ~ pmax(0, abs(U - mu) - 0.5),
        alternative == "less"      ~ (U - mu) + 0.5,
        alternative == "greater"   ~ (U - mu) - 0.5
      ),
      dev_uncorr = case_when(
        alternative == "two.sided" ~ abs(U - mu),
        alternative == "less"      ~ (U - mu),
        alternative == "greater"   ~ (U - mu)
      ),
      
      # P-values
      # No Ties
      z_corr_no_ties = dev_cc / sigma_no_ties,
      p_corr_no_ties = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_corr_no_ties, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_corr_no_ties, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_corr_no_ties, lower.tail = FALSE)
      ),
      
      # No Ties (Uncorrected)
      z_uncorr_no_ties = dev_uncorr / sigma_no_ties,
      p_uncorr_no_ties = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_uncorr_no_ties, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_uncorr_no_ties, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_uncorr_no_ties, lower.tail = FALSE)
      ),
      
      # Tied (Corrected)
      z_corr_tied = dev_cc / sigma_one_tie,
      p_corr_tied = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_corr_tied, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_corr_tied, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_corr_tied, lower.tail = FALSE)
      ),
      
      # Tied (Uncorrected)
      z_uncorr_tied = dev_uncorr / sigma_one_tie,
      p_uncorr_tied = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_uncorr_tied, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_uncorr_tied, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_uncorr_tied, lower.tail = FALSE)
      )
    ) %>%
    # Final Clamp: Ensure p <= 1 (Standard behavior)
    mutate(across(starts_with("p_"), ~ pmin(1, .))) %>%
    select(U, is_integer, starts_with("p_"))
  
  return(results_df)
}

#' GRIM-U Forensic Consistency Check
#'
#' Checks the mathematical consistency of reported Mann-Whitney U test results.
#' Performs a "Triangulation" check across Sample Sizes (N), Test Statistic (U), and P-value (P).
#'
#' 1. N <-> U: Is the reported U physically possible? (Bounds & Granularity).
#' 2. N <-> P: Is the reported P mathematically possible? (Grid search of attainable p-values).
#' 3. U <-> P: Does the specific reported U generate the reported P?
#'
#' @param n1 Integer. Sample size of group 1.
#' @param n2 Integer. Sample size of group 2.
#' @param u_reported Numeric, optional. The reported U statistic. If NA, checks P-value consistency only.
#' @param p_reported Numeric. The reported p-value.
#' @param comparison Character. "equal" (default) or "less_than" (for reported inequalities like "p < .05").
#' @param digits Integer. The number of decimal places the p-value was reported to. Mandatory if p_reported is supplied.
#' @param rounding Character vector, optional. Allowed rounding methods: "round", "trunc", "up".
#'   If NULL, defaults to c("round", "trunc").
#' @param p_min Numeric, optional. Manually override the lower bound of the p-value search window.
#' @param p_max Numeric, optional. Manually override the upper bound of the p-value search window.
#' @param alternative Character. One of "two.sided", "less", or "greater". Defaults to "two.sided".
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
  
  # Calculate Z boundaries for the p-value window
  # Note: 
  # Small P -> Large Z magnitude (Deep Tail)
  # Large P -> Small Z magnitude (Near Mean)
  
  z_deep <- p_to_z(if(is.na(p_min_search)) 0 else p_min_search, alternative)
  z_shallow <- p_to_z(p_max_search, alternative)
  
  # Convert Z to U
  # U ~ mu + Z*sigma
  u_bound_1 <- mu + z_deep * sigma_est
  u_bound_2 <- mu + z_shallow * sigma_est
  
  # NA Safety: If bounds are NaN (e.g., extremely far out), default to full range
  if (is.na(u_bound_1) || is.nan(u_bound_1)) u_bound_1 <- if(alternative=="greater") max_u else 0
  if (is.na(u_bound_2) || is.nan(u_bound_2)) u_bound_2 <- if(alternative=="greater") mu else mu
  
  # Sort bounds
  raw_start <- min(u_bound_1, u_bound_2)
  raw_end   <- max(u_bound_1, u_bound_2)
  
  # Pad and Clamp
  u_start_est <- floor(raw_start - 2)
  u_end_est   <- ceiling(raw_end + 2)
  
  # Physical Clamping
  u_start_est <- max(0, u_start_est)
  u_end_est   <- min(max_u, u_end_est)
  
  # --- Step 4: Call Engine ---
  results_df <- grimu_map_pvalues(n1, n2, u_min = u_start_est, u_max = u_end_est, 
                                  alternative = alternative)
  
  # --- Step 5: Check P-Value Consistency ---
  if (is.null(rounding)) {
    methods <- c("round", "trunc") 
  } else {
    methods <- match.arg(rounding, c("round", "trunc", "up"), several.ok = TRUE)
  }
  
  epsilon <- 10^(-digits)
  buffer <- 1e-14 
  
  is_match <- function(val) {
    if (is.na(val)) return(FALSE)
    if ("round" %in% methods) {
      if (val >= (p_reported - 0.5 * epsilon - buffer) && 
          val <  (p_reported + 0.5 * epsilon - buffer)) return(TRUE)
    }
    if ("trunc" %in% methods) {
      if (val >= (p_reported - buffer) && 
          val <  (p_reported + epsilon - buffer)) return(TRUE)
    }
    if ("up" %in% methods) {
      if (val >  (p_reported - epsilon + buffer) && 
          val <= (p_reported + buffer)) return(TRUE)
    }
    return(FALSE)
  }
  
  check_col <- function(col_val) {
    if (comparison == "equal") {
      return(is_match(col_val))
    } else {
      return(!is.na(col_val) && col_val < p_reported)
    }
  }
  
  results_checked <- results_df %>%
    rowwise() %>%
    mutate(
      valid_exact          = check_col(p_exact),
      valid_corr_no_ties   = check_col(p_corr_no_ties),
      valid_uncorr_no_ties = check_col(p_uncorr_no_ties),
      valid_corr_tied      = check_col(p_corr_tied),
      valid_uncorr_tied    = check_col(p_uncorr_tied),
      
      is_consistent = valid_exact | 
        valid_corr_no_ties | valid_uncorr_no_ties | 
        valid_corr_tied | valid_uncorr_tied
    ) %>%
    ungroup() %>%
    # Filter for display (Diagnostic mode)
    # Show row if IT IS consistent OR if ANY p-value is in the general "search window"
    # This reveals near-misses for all 5 methods, not just tied-corrected.
    filter(
      is_consistent | 
        if_any(starts_with("p_"), 
               ~ . >= (if(is.na(p_min_search)) 0 else p_min_search) & 
                 . <= p_max_search)
    )
  
  # --- Step 6: Triangulate U and P ---
  
  # A. Is P possible at all?
  p_consistent <- any(results_checked$is_consistent)
  
  # B. Is the reported U consistent with the reported P?
  # We check if u_reported appears in the rows where 'is_consistent' is TRUE.
  if (!is.na(u_reported)) {
    # We use a small epsilon for float comparison of U values
    u_consistent <- results_checked %>%
      filter(is_consistent) %>%
      filter(abs(U - u_reported) < 1e-10) %>%
      nrow() > 0
  } else {
    u_consistent <- NA
  }
  
  # C. Overall Consistency
  # If U is reported: Must be physically possible AND match P.
  # If U is missing:  P must be possible.
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
    p_range_min = if(is.na(p_min_search)) 0 else p_min_search,
    p_range_max = p_max_search,
    u_range_min = u_start_est,
    u_range_max = u_end_est,
    
    # Component Checks
    u_bounds_consistent = u_res$u_bounds_consistent,
    u_granularity_consistent = u_res$u_granularity_consistent,
    u_matches_p = u_consistent, # Does Reported U -> Reported P?
    
    # NB No p_bounds_consistent as its tautological, p min and max are calculated from reported p
    p_granularity_consistent = p_consistent,  # Is Reported P possible at all?
    
    # Debug Flags
    p_matches_exact = any(results_checked$valid_exact),
    p_matches_no_ties = any(results_checked$valid_corr_no_ties | results_checked$valid_uncorr_no_ties),
    p_matches_ties = any(results_checked$valid_corr_tied | results_checked$valid_uncorr_tied)
  )
  
  return(list(summary = summary_df, details = results_checked))
}

#' Calculate P-Value Saturation (Granularity Density)
#'
#' Calculates the proportion of all possible rounded p-values (at a given precision) that are
#' mathematically attainable by a Mann-Whitney U test with sample sizes n1 and n2.
#' Used to estimate the False Positive Rate (FPR) of a GRIM-U check for a specific N.
#'
#' @param n1 Integer. Sample size of group 1.
#' @param n2 Integer. Sample size of group 2.
#' @param decimals Integer. The number of decimal places to test (e.g., 2 for p=0.04). Defaults to 3.
#' @param p_lower_threshold Numeric. Lower bound of the p-value range to check (default 0).
#' @param p_upper_threshold Numeric. Upper bound of the p-value range to check (default 1).
#'
#' @return A numeric value between 0 and 1 representing the saturation ratio:
#'   (Unique attainable rounded p-values) / (Total possible bins in range).
#' @export
grimu_saturation <- function(n1, n2, decimals = 3, p_lower_threshold = 0, p_upper_threshold = 1) {
  
  # 1. CALL THE ENGINE 
  # Get all U values from 0 to Mean (Sufficient because distribution is symmetric)
  # This covers every possible unique p-value the test can produce.
  mu <- (n1 * n2) / 2
  p_space <- grimu_map_pvalues(n1, n2, u_min = 0, u_max = floor(mu))
  
  # 2. Reshape and Filter
  unique_rounded_p <- p_space %>%
    pivot_longer(cols = starts_with("p_"), values_to = "p_val") %>%
    filter(!is.na(p_val)) %>%
    filter(p_val >= p_lower_threshold & p_val <= p_upper_threshold) %>%
    # Round to target precision
    mutate(p_rounded = roundwork::round_up(p_val, decimals)) %>%
    distinct(p_rounded) %>%
    nrow()
  
  # 3. Calculate Coverage
  total_slots <- length(seq(p_lower_threshold, p_upper_threshold, by = 10^-decimals))
  
  return(unique_rounded_p / total_slots)
}
