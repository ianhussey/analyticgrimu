#' Enumerate achievable tie-correction K values
#'
#' Returns the set of K = sum(t_i^3 - t_i) values reachable by any tie-group
#' partition of a sample of size N, where each tie group has size t_i >= 2 and
#' the tie-group sizes sum to at most N. K = 0 (no ties) is always included.
#' This is the family of values that drive the tie-corrected variance of the
#' Mann-Whitney U statistic.
#'
#' Enumeration is capped at `max_partition_size` (default 15) to keep cost
#' bounded when N is large. Tie groups larger than this are uncommon in
#' typical behavioural data, and large-K corrections to sigma have
#' diminishing effect because the (N+1)/12 leading term dominates as N grows.
#' At the default cap there are at most ~115 unique K values regardless of N.
#'
#' Results are memoised per (N, max_partition_size) to avoid recomputing the
#' partition enumeration on every call to `grimu_map_pvalues()`.
#'
#' @param N Integer. Total sample size (n1 + n2).
#' @param max_partition_size Integer. Largest tie-pattern total size to
#'   enumerate. Defaults to 15.
#'
#' @return A sorted numeric vector of unique K values (including 0).
#' @keywords internal
#' @noRd
# Memoisation cache for achievable_K_values(), keyed by "N_max_partition_size".
.K_cache <- new.env(parent = emptyenv())

achievable_K_values <- function(N, max_partition_size = 15) {
  key <- paste(N, max_partition_size, sep = "_")
  cached <- .K_cache[[key]]
  if (!is.null(cached)) {
    return(cached)
  }

  Ks <- 0 # always include the no-ties case
  s_max <- min(N, max_partition_size)
  if (s_max >= 2) {
    for (s in 2:s_max) {
      p_mat <- partitions::parts(s)
      for (j in seq_len(ncol(p_mat))) {
        parts_j <- p_mat[, j]
        parts_j <- parts_j[parts_j > 0]
        if (all(parts_j >= 2)) {
          Ks <- c(Ks, sum(parts_j^3 - parts_j))
        }
      }
    }
  }
  result <- sort(unique(Ks))
  .K_cache[[key]] <- result
  result
}

#' Generate Possible Mann-Whitney U P-values
#'
#' The core engine for GRIM-U. Generates a grid of valid U statistics (integers and half-integers)
#' within a specified range and calculates the corresponding p-values under every reporting
#' convention the package models:
#' 1. Exact (R / SciPy convention): doubled smaller tail using `pwilcox()`.
#' 2. Exact (SPSS Exact / StatXact convention): deviation-from-mean two-sided p,
#'    `P(|W - mu_W| >= |w_obs - mu_W|)`. Two-sided integer U only.
#' 3. Mid-p: `P(W > w) + 0.5 * P(W = w)`. Integer U only.
#' 4. Asymptotic (corrected and uncorrected) under every achievable tie-correction K value
#'    via `achievable_K_values()`. K = 0 is the no-ties variance; K = 6 is one tied pair;
#'    larger K values cover multi-tie patterns. Columns are named `p_corr_K{value}` and
#'    `p_uncorr_K{value}`.
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
#'   \item{p_exact}{Exact p-value, R/SciPy doubled-tail convention (NA for non-integers).}
#'   \item{p_exact_dev}{Exact p-value, SPSS Exact / StatXact deviation-from-mean convention
#'     (two-sided integer U only; NA otherwise).}
#'   \item{p_mid}{Mid-p value (integer U only; NA otherwise).}
#'   \item{p_corr_K*}{Asymptotic continuity-corrected p-value under tie-correction K; one column per achievable K, named \code{p_corr_K} followed by the K value.}
#'   \item{p_uncorr_K*}{Asymptotic uncorrected p-value under tie-correction K; one column per achievable K, named \code{p_uncorr_K} followed by the K value.}
#' @examples
#' # Attainable U values and their exact p-values for two groups of 8 and 9:
#' map <- grimu_map_pvalues(n1 = 8, n2 = 9)
#' head(map[, c("U", "is_integer", "p_exact")])
#' @export
grimu_map_pvalues <- function(
  n1,
  n2,
  u_min = NULL,
  u_max = NULL,
  alternative = "two.sided"
) {
  # Safety Check: Input Validation
  if (is.na(n1) || is.na(n2)) {
    return(tibble(U = numeric(), is_integer = logical()))
  }

  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  max_u <- n1 * n2

  # Default bounds logic
  if (is.null(u_min) || is.null(u_max)) {
    if (alternative == "greater") {
      if (is.null(u_min)) {
        u_min <- floor(mu)
      }
      if (is.null(u_max)) u_max <- max_u
    } else {
      if (is.null(u_min)) {
        u_min <- 0
      }
      if (is.null(u_max)) u_max <- ceiling(mu)
    }
  }

  # Enforce Half-Integer Lattice
  # Snap arbitrary bounds to the nearest valid U values (integers or half-integers).
  if (is.null(u_min)) {
    u_min <- 0
  }
  if (is.null(u_max)) {
    u_max <- max_u
  }

  u_start <- ceiling(u_min * 2) / 2
  u_end <- floor(u_max * 2) / 2

  # Safety Clamp
  u_start <- max(0, u_start)
  u_end <- min(max_u, u_end)

  # Safety Check: If start > end, return empty tibble immediately
  if (u_start > u_end) {
    return(tibble(U = numeric(), is_integer = logical()))
  }

  vals <- roundwork::round_up(seq(u_start, u_end, by = 0.5), 1)

  # Enumerate all achievable tie-correction K values for this N
  K_vals <- achievable_K_values(N)

  # Base tibble + deviations
  results_df <- tibble(U = vals) %>%
    mutate(
      is_integer = (U %% 1 == 0),

      dev_cc = case_when(
        alternative == "two.sided" ~ pmax(0, abs(U - mu) - 0.5),
        alternative == "less" ~ (U - mu) + 0.5,
        alternative == "greater" ~ (U - mu) - 0.5
      ),
      dev_uncorr = case_when(
        alternative == "two.sided" ~ abs(U - mu),
        alternative == "less" ~ (U - mu),
        alternative == "greater" ~ (U - mu)
      ),

      # Exact (R / SciPy doubled-tail convention)
      p_exact = if_else(
        is_integer,
        case_when(
          alternative == "less" ~ pwilcox(U, n1, n2),
          alternative == "greater" ~ pwilcox(U - 1, n1, n2, lower.tail = FALSE),
          alternative == "two.sided" ~ {
            p_lower <- pwilcox(U, n1, n2)
            p_upper <- pwilcox(U - 1, n1, n2, lower.tail = FALSE)
            2 * pmin(p_lower, p_upper)
          }
        ),
        NA_real_
      ),

      # Exact (SPSS Exact / StatXact deviation-from-mean convention).
      # Defined only for two-sided integer U. At U = mu the raw sum equals
      # 1 + dwilcox(mu) due to discrete overlap at W = mu; the final
      # pmin(1, .) below clips that case to exactly 1, which matches the
      # P(|W - mu| >= 0) = 1 interpretation.
      p_exact_dev = if (alternative == "two.sided") {
        if_else(
          is_integer,
          {
            dev <- abs(U - mu)
            pwilcox(mu - dev, n1, n2) +
              pwilcox(mu + dev - 1, n1, n2, lower.tail = FALSE)
          },
          NA_real_
        )
      } else {
        NA_real_
      },

      # Mid-p: P(W > w) + 0.5 * P(W = w). Integer U only.
      # mid_lower + mid_upper = 1 exactly, so the two-sided 2*pmin never overshoots.
      p_mid = if_else(
        is_integer,
        case_when(
          alternative == "less" ~ pwilcox(U, n1, n2) - 0.5 * dwilcox(U, n1, n2),
          alternative == "greater" ~ pwilcox(
            U - 1,
            n1,
            n2,
            lower.tail = FALSE
          ) -
            0.5 * dwilcox(U, n1, n2),
          alternative == "two.sided" ~ {
            mid_lower <- pwilcox(U, n1, n2) - 0.5 * dwilcox(U, n1, n2)
            mid_upper <- pwilcox(U - 1, n1, n2, lower.tail = FALSE) -
              0.5 * dwilcox(U, n1, n2)
            2 * pmin(mid_lower, mid_upper)
          }
        ),
        NA_real_
      )
    )

  # Asymptotic p-values: one (corr, uncorr) pair per achievable K.
  # All K values are applied to both integer and half-integer U candidates;
  # this is over-coverage (some K values strictly require an even-sized tie
  # group and thus apply only to half-integer U; some apply only to integer
  # U), but a forensic check should be permissive on the consistent side.
  #
  # Vectorised via outer() to avoid an R-level loop over K_vals: at large N
  # there can be ~100 K values and ~thousands of U candidates, and the
  # column-by-column assignment was the dominant cost in the previous loop.
  sigma_vec <- sqrt(
    (n1 * n2 * (N + 1)) / 12 - (n1 * n2 * K_vals) / (12 * N * (N - 1))
  )
  # NaN protection: sigma_K = 0 (degenerate all-tied case) gives Inf z and
  # p = 0 after pnorm, which is harmless; pmin(1, .) below leaves it at 0.
  z_corr_mat <- outer(results_df$dev_cc, sigma_vec, "/")
  z_uncorr_mat <- outer(results_df$dev_uncorr, sigma_vec, "/")
  if (alternative == "two.sided") {
    p_corr_mat <- 2 * pnorm(z_corr_mat, lower.tail = FALSE)
    p_uncorr_mat <- 2 * pnorm(z_uncorr_mat, lower.tail = FALSE)
  } else if (alternative == "less") {
    p_corr_mat <- pnorm(z_corr_mat, lower.tail = TRUE)
    p_uncorr_mat <- pnorm(z_uncorr_mat, lower.tail = TRUE)
  } else {
    # "greater"
    p_corr_mat <- pnorm(z_corr_mat, lower.tail = FALSE)
    p_uncorr_mat <- pnorm(z_uncorr_mat, lower.tail = FALSE)
  }
  colnames(p_corr_mat) <- paste0("p_corr_K", K_vals)
  colnames(p_uncorr_mat) <- paste0("p_uncorr_K", K_vals)
  results_df <- bind_cols(
    results_df,
    tibble::as_tibble(p_corr_mat),
    tibble::as_tibble(p_uncorr_mat)
  )

  # Final cap at p <= 1.
  # This is required for both exact two-sided columns (p_exact and p_exact_dev),
  # which can sum to slightly over 1 at central U due to the discrete
  # distribution's overlap at W = mu. The asymptotic z columns (p_corr_K*,
  # p_uncorr_K*) never overshoot, because dev_cc = max(0, |U-mu| - 0.5) and
  # dev_uncorr = |U-mu| are both >= 0, so 2*pnorm(z, lower.tail = FALSE) <= 1
  # in the two-sided case and pnorm(z) is in [0, 1] in the one-sided cases.
  # The mid-p column never overshoots either (mid_lower + mid_upper = 1
  # exactly). The cap is applied uniformly across all p_* columns as
  # defensive coding, but only the exact columns ever activate it.
  results_df <- results_df %>%
    mutate(across(starts_with("p_"), ~ pmin(1, .))) %>%
    select(U, is_integer, starts_with("p_"))

  return(results_df)
}
