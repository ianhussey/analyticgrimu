#' analyticgrimu: Granularity Testing for Mann-Whitney U Test Results
#'
#' An implementation of the GRIM-U granularity test of Heathers and Grimes
#' (2025), a forensic meta-science tool for checking whether reported
#' Mann-Whitney U test statistics and p-values are mathematically possible
#' given the group sample sizes and the rounding of the reported values.
#'
#' @keywords internal
#' @importFrom stats pwilcox dwilcox pnorm qnorm
"_PACKAGE"

# Quiet R CMD check notes about the non-standard evaluation used inside the
# dplyr/ggplot2 pipelines (bare column names have no visible binding).
utils::globalVariables(c(
  ".",
  "U",
  "is_integer",
  "dev_cc",
  "dev_uncorr",
  "p_exact",
  "p_exact_dev",
  "p_mid",
  "p_corr_K0",
  "p_uncorr_K0",
  "valid_exact",
  "valid_exact_dev",
  "valid_mid",
  "valid_corr_no_ties",
  "valid_uncorr_no_ties",
  "valid_corr_ties",
  "valid_uncorr_ties",
  "is_consistent",
  "val",
  "y_row",
  "x_col",
  "is_possible",
  "facet_label"
))
