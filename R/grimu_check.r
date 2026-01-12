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

