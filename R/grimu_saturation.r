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
#' @examples
#' # Proportion of 3-decimal p-values attainable for two groups of 8 and 9:
#' grimu_saturation(n1 = 8, n2 = 9, decimals = 3)
#' @export
grimu_saturation <- function(
  n1,
  n2,
  decimals = 3,
  p_lower_threshold = 0,
  p_upper_threshold = 1
) {
  # 1. CALL THE ENGINE
  # Get all U values from 0 to Mean (Sufficient because distribution is symmetric)
  # This covers every possible unique p-value the test can produce.
  #
  # Note: grimu_map_pvalues() now returns p-values under multiple methods
  # (exact R/SciPy, exact SPSS/StatXact, mid-p) and under every achievable
  # tie-correction K value (one (corr, uncorr) column pair per K). The
  # pivot_longer + distinct(p_rounded) below correctly unions over all of
  # those columns, so the saturation here reflects the full attainable
  # p-value set across every reporting convention this package models.
  mu <- (n1 * n2) / 2
  p_space <- grimu_map_pvalues(n1, n2, u_min = 0, u_max = floor(mu))

  # 2. Flatten all p_* columns into a single vector and reduce to unique
  # rounded values. We avoid pivot_longer here because at large N the
  # candidate tibble is ~thousands of U rows by ~hundreds of p_* columns,
  # and pivoting that to long form was the dominant cost.
  p_cols <- grep("^p_", names(p_space), value = TRUE)
  p_vec <- unlist(p_space[p_cols], use.names = FALSE)
  p_vec <- p_vec[
    !is.na(p_vec) &
      p_vec >= p_lower_threshold &
      p_vec <= p_upper_threshold
  ]
  unique_rounded_p <- length(unique(roundwork::round_up(p_vec, decimals)))

  # 3. Calculate Coverage
  total_slots <- length(seq(
    p_lower_threshold,
    p_upper_threshold,
    by = 10^-decimals
  ))

  return(unique_rounded_p / total_slots)
}
