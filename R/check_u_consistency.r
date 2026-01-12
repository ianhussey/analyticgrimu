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
check_u_consistency <- function(n1, n2, u_reported) {
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

