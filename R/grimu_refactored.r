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
  
  vals <- roundwork::round_up(seq(u_start, u_end, by = 0.5), 1)
  
  # --- Standard Errors ---
  sigma_no_ties <- sqrt((n1 * n2 * (N + 1)) / 12)
  correction_term <- (n1 * n2 * 6) / (12 * N * (N - 1))
  sigma_one_tie <- sqrt((n1 * n2 * (N + 1)) / 12 - correction_term)
  
  # --- Vectorized Calculation ---
  results_df <- tibble(U = vals) %>%
    mutate(
      is_integer = (U %% 1 == 0),
      sigma = if_else(is_integer, sigma_no_ties, sigma_one_tie),
      
      # 1. Exact Method (Symmetric logic handles U > mean)
      p_exact = if_else(
        is_integer,
        2 * pwilcox(if_else(U < mu, U, n1 * n2 - U), n1, n2),
        NA_real_
      ),
      
      # 2. Asymptotic Corrected
      z_corrected = (abs(U - mu) - 0.5) / sigma,
      p_asymp_corrected = 2 * pnorm(z_corrected, lower.tail = FALSE),
      
      # 3. Asymptotic Uncorrected
      z_uncorrected = abs(U - mu) / sigma,
      p_asymp_uncorrected = 2 * pnorm(z_uncorrected, lower.tail = FALSE)
    )
  
  return(results_df)
}

# forensic tool
grimu_check <- function(n1, n2, p_reported, comparison = "equal", digits = 2, 
                        p_min = NULL, p_max = NULL) {
  
  # --- 1. Smart Range Detection ---
  if (is.null(p_min) || is.null(p_max)) {
    window <- 10^(-digits) * 5 
    if (comparison == "less_than") {
      p_max_search <- p_reported + window 
      p_min_search <- 1e-7 # Safe non-zero lower bound
    } else {
      p_max_search <- p_reported + window
      p_min_search <- max(1e-7, p_reported - window)
    }
  } else {
    p_min_search <- p_min
    p_max_search <- p_max
  }
  
  # --- 2. Determine U Bounds ---
  # We need the SE to estimate bounds, so we quickly calc it or grab from engine
  # (Calculating locally is faster than calling engine just for sigma)
  N <- n1 + n2
  sigma_est <- sqrt((n1 * n2 * (N + 1)) / 12)
  mu <- (n1 * n2) / 2
  
  # Z bounds to U bounds
  z_bounds <- qnorm(1 - c(p_min_search, p_max_search) / 2)
  u_deviation_max <- abs(z_bounds[1] * sigma_est)
  
  u_start_est <- floor(mu - u_deviation_max - 2)
  u_end_est   <- ceiling(mu + u_deviation_max + 2)
  
  # --- 3. CALL THE ENGINE ---
  results_df <- grimu_map_pvalues(n1, n2, u_min = u_start_est, u_max = u_end_est)
  
  # --- 4. Consistency Logic (Same as before) ---
  check_col <- function(col_val, p_rep, comp, dig) {
    if (comp == "equal") {
      # roundwork::round_up(col_val, dig) == p_rep
      round(col_val + 1e-10, dig) == p_rep 
    } else {
      col_val < p_rep
    }
  }
  
  results_checked <- results_df %>%
    rowwise() %>%
    mutate(
      valid_exact = !is.na(p_exact) && check_col(p_exact, p_reported, comparison, digits),
      valid_corrected = check_col(p_asymp_corrected, p_reported, comparison, digits),
      valid_uncorrected = check_col(p_asymp_uncorrected, p_reported, comparison, digits),
      is_consistent = valid_exact | valid_corrected | valid_uncorrected
    ) %>%
    ungroup() %>%
    filter(
      is_consistent | 
        (p_asymp_uncorrected >= p_min_search & p_asymp_uncorrected <= p_max_search)
    )
  
  summary_df <- tibble(
    n1 = n1, n2 = n2, p_reported = p_reported,
    consistent = any(results_checked$is_consistent),
    matches_exact = any(results_checked$valid_exact),
    matches_r_corrected = any(results_checked$valid_corrected),
    matches_spss_uncorrected = any(results_checked$valid_uncorrected)
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