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
