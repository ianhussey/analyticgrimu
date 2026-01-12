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


