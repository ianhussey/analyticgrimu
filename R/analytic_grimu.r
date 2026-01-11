library(tidyverse)
library(roundwork)

# core engine
grimu_map_pvalues <- function(n1, n2, u_min = 0, u_max = NULL) {
  
  # --- Constants ---
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  
  # Default to covering the lower half (sufficient for saturation) if max not specified
  if (is.null(u_max)) u_max <- floor(mu)
  
  # Ensure bounds are safe
  u_start <- max(0, u_min)
  u_end <- min(n1 * n2, u_max)
  
  # Use roundwork to ensure clean steps (avoid floating point drift)
  vals <- roundwork::round_up(seq(u_start, u_end, by = 0.5), 1)
  
  # --- Standard Errors ---
  # 1. Sigma assuming NO ties (Max variance)
  sigma_no_ties <- sqrt((n1 * n2 * (N + 1)) / 12)
  
  # 2. Sigma assuming ONE pair of ties (Minimal variance reduction)
  correction_term <- (n1 * n2 * 6) / (12 * N * (N - 1))
  sigma_one_tie <- sqrt((n1 * n2 * (N + 1)) / 12 - correction_term)
  
  # --- Vectorized Calculation ---
  results_df <- tibble(U = vals) %>%
    mutate(
      is_integer = (U %% 1 == 0),
      
      # --- A. Exact Method ---
      # STRICT: Only valid for Integers. 
      # pwilcox is mathematically undefined for fractional quantiles.
      p_exact = if_else(
        is_integer,
        2 * pwilcox(if_else(U < mu, U, n1 * n2 - U), n1, n2),
        NA_real_
      ),
      
      # --- B. Asymptotic (NO TIES Variance) ---
      # PERMISSIVE: Calculated for ALL U (Integer and Fractional).
      # Catches cases where researcher has ties (Fractional U) but used wrong Sigma.
      z_corr_no_ties   = pmax(0, abs(U - mu) - 0.5) / sigma_no_ties, # clamp continuity correction
      p_corr_no_ties   = 2 * pnorm(z_corr_no_ties, lower.tail = FALSE),
      
      z_uncorr_no_ties = abs(U - mu) / sigma_no_ties,
      p_uncorr_no_ties = 2 * pnorm(z_uncorr_no_ties, lower.tail = FALSE),
      
      # --- C. Asymptotic (TIES Variance) ---
      # PERMISSIVE: Calculated for ALL U.
      # Catches cases where researcher has Integer U but hidden ties.
      z_corr_tied      = pmax(0, abs(U - mu) - 0.5) / sigma_one_tie, # clamp continuity correction
      p_corr_tied      = 2 * pnorm(z_corr_tied, lower.tail = FALSE),
      
      z_uncorr_tied    = abs(U - mu) / sigma_one_tie,
      p_uncorr_tied    = 2 * pnorm(z_uncorr_tied, lower.tail = FALSE)
    )
  
  return(results_df)
}

# forensic tool
grimu_check <- function(n1, n2, p_reported, comparison = "equal", digits = 2, 
                        p_min = NULL, p_max = NULL) {
  
  # 1. Range Detection 
  if (is.null(p_min) || is.null(p_max)) {
    window <- 10^(-digits) * 5 
    p_max_search <- p_reported + window
    
    # Calculate raw lower bound
    raw_p_min <- p_reported - window
    
    if (comparison == "less_than") {
      p_min_search <- NA_real_ # Flag to trigger U=0 anchor
    } else if (raw_p_min <= 0) {
      # If the window touches or goes below zero, we can't use qnorm.
      # We flag this to anchor the U-search to 0 manually.
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
  
  # 2. U Bounds 
  # We need the SE to estimate bounds, so we quickly calc it or grab from engine
  # Calculating locally is faster than calling engine just for sigma
  N <- n1 + n2
  sigma_est <- sqrt((n1 * n2 * (N + 1)) / 12)
  mu <- (n1 * n2) / 2
  
  # Z bounds to U bounds
  # Upper U bound (unchanged)
  z_min <- qnorm(1 - p_max_search / 2)
  u_dev_min <- abs(z_min * sigma_est)
  u_end_est <- min(floor(mu), ceiling(mu - u_dev_min + 2)) 
  
  # Lower U bound 
  if (is.na(p_min_search)) {
    # If p_min was < 0 or "less_than", search the WHOLE tail down to 0.
    u_start_est <- 0 
  } else {
    # Standard Z-based search
    z_max <- qnorm(1 - p_min_search / 2)
    u_dev_max <- abs(z_max * sigma_est)
    u_start_est <- floor(mu - u_dev_max - 2)
  }
  
  # 3. Call Updated Engine
  results_df <- grimu_map_pvalues(n1, n2, u_min = u_start_est, u_max = u_end_est)
  
  # 4. Consistency Logic (Expanded for 5 columns)
  check_col <- function(col_val, p_rep, comp, dig) {
    # Check for NA (e.g. no-tie calc on fractional U)
    if (is.na(col_val)) return(FALSE)
    
    if (comp == "equal") {
      roundwork::round_up(col_val, dig) == p_rep
    } else {
      col_val < p_rep
    }
  }
  
  results_checked <- results_df %>%
    rowwise() %>%
    mutate(
      # Check Exact
      valid_exact = check_col(p_exact, p_reported, comparison, digits),
      
      # Check No-Ties Scenarios (Only calculated for integers)
      valid_corr_no_ties   = check_col(p_corr_no_ties, p_reported, comparison, digits),
      valid_uncorr_no_ties = check_col(p_uncorr_no_ties, p_reported, comparison, digits),
      
      # Check Tied Scenarios (Calculated for everyone)
      valid_corr_tied      = check_col(p_corr_tied, p_reported, comparison, digits),
      valid_uncorr_tied    = check_col(p_uncorr_tied, p_reported, comparison, digits),
      
      # Overall Consistency: Match ANY of the 5
      is_consistent = valid_exact | 
        valid_corr_no_ties | valid_uncorr_no_ties | 
        valid_corr_tied | valid_uncorr_tied
    ) %>%
    ungroup() %>%
    # Filter for display (Diagnostic mode)
    filter(
      is_consistent | 
        (p_corr_tied >= p_min_search & p_corr_tied <= p_max_search)
    )
  
  summary_df <- tibble(
    n1 = n1, n2 = n2, p_reported = p_reported,
    consistent = any(results_checked$is_consistent),
    # Flags for diagnosis
    matches_exact = any(results_checked$valid_exact),
    matches_no_ties = any(results_checked$valid_corr_no_ties | results_checked$valid_uncorr_no_ties),
    matches_ties = any(results_checked$valid_corr_tied | results_checked$valid_uncorr_tied)
  )
  
  return(list(summary = summary_df, details = results_checked))
}

# saturation calc
grimu_saturation <- function(n1, n2, decimals = 3, p_threshold = 0.05) {
  
  # 1. CALL THE ENGINE 
  # Get all U values from 0 to Mean (Sufficient because distribution is symmetric)
  # This covers every possible unique p-value the test can produce.
  mu <- (n1 * n2) / 2
  p_space <- grimu_map_pvalues(n1, n2, u_min = 0, u_max = floor(mu))
  
  # 2. Reshape and Filter
  unique_rounded_p <- p_space %>%
    pivot_longer(cols = starts_with("p_"), values_to = "p_val") %>%
    filter(!is.na(p_val)) %>%
    filter(p_val <= p_threshold) %>%
    # Round to target precision
    mutate(p_rounded = round(p_val, decimals)) %>%
    distinct(p_rounded) %>%
    nrow()
  
  # 3. Calculate Coverage
  total_slots <- length(seq(0, p_threshold, by = 10^-decimals))
  
  return(unique_rounded_p / total_slots)
}