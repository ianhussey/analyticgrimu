#' Plot GRIMU Probability Space for Multiple Sample Size Pairs
#'
#' @description 
#' Visualizes the "granularity" of p-values for the Mann-Whitney U test across 
#' multiple sample size combinations. The function generates a heat-map style 
#' grid where each tile represents a possible 3-decimal p-value. This is useful 
#' for meta-scientific screening to identify impossible reported p-values.
#'
#' @param n1_vec A numeric vector of sample sizes for group 1.
#' @param n2_vec A numeric vector of sample sizes for group 2. Must be the 
#'   same length as `n1_vec`.
#'
#' @return A `ggplot` object displaying a faceted grid of p-value probability 
#'   spaces. Black tiles indicate mathematically possible p-values (rounded up 
#'   to 3 decimals), while gray tiles indicate impossible values.
#'
#' @details 
#' The function uses `grimu_map_pvalues` to calculate exact p-values for all 
#' possible U statistics given `n1` and `n2`. It then maps these onto a 10x100 
#' grid where the y-axis represents the first decimal place and the x-axis 
#' represents the second and third decimal places.
#' 
#' Note that this function assumes the absence of ties and uses exact 
#' probability calculations.
#'
#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @importFrom roundwork round_up
#' @importFrom tibble tibble
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' plot_grimu_probability_multi(n1_vec = c(8, 30), n2_vec = c(9, 31))
#' }
plot_grimu_probability_multi <- function(n1_vec, n2_vec) {
  
  if (length(n1_vec) != length(n2_vec)) {
    stop("n1 and n2 vectors must be of equal length.")
  }
  
  # 1. Iterate over inputs to build the master dataset
  plot_data <- map2_dfr(n1_vec, n2_vec, function(n1, n2) {
    
    # Calculate Max U for this pair
    max_u <- floor((n1 * n2) / 2)
    
    # Get all valid p-values (Exact)
    possible_p_data <- grimu_map_pvalues(n1, n2, u_min = 0, u_max = max_u)
    
    valid_p_values <- possible_p_data %>%
      filter(!is.na(p_exact)) %>%
      pull(p_exact)
    
    possible_values <- unique(roundwork::round_up(valid_p_values, 3))
    
    # Generate Grid for this specific pair
    tibble(
      val = seq(0, 0.999, by = 0.001)
    ) %>%
      mutate(
        y_row = floor(val * 10), 
        x_col = (val * 1000) %% 100,
        is_possible = roundwork::round_up(val, 3) %in% possible_values,
        # Create the Label
        facet_label = paste0("n1 = ", n1, ", n2 = ", n2)
      )
  })
  
  # 2. preserve input order
  # By default, facets sort alphabetically (so "n1=10" comes before "n1=4").
  # We force the factor levels to match the user's input order.
  ordered_labels <- paste0("n1 = ", n1_vec, ", n2 = ", n2_vec)
  plot_data$facet_label <- factor(plot_data$facet_label, levels = unique(ordered_labels))
  
  # 3. plot
  p_grid <- ggplot(plot_data, aes(x = x_col, y = y_row, fill = is_possible)) +
    geom_tile(color = "white", linewidth = 0.1) + 
    
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "gray95")) +
    
    scale_y_reverse(
      breaks = 0:9, 
      labels = paste0("0.", 0:9), 
      name = "First Decimal Place"
    ) +
    
    scale_x_continuous(
      breaks = seq(0, 99, by = 10),
      labels = function(x) paste0(".0", x),
      name = "Second & Third Decimals"
    ) +
    
    labs(
      title = "The P-Value Probability Space",
      subtitle = "Black squares represent mathematically possible p-values (exact test, no ties).\nGray squares represent impossible values.",
      fill = "Possible?"
    ) +
    
    coord_fixed(ratio = 1) + 
    facet_wrap(~facet_label, ncol = 1) +
    
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30", size = 10),
      strip.text = element_text(face = "bold", size = 11), # Style the facet labels
      legend.position = "bottom"
    )
  
  return(p_grid)
}

# Example Usage:
plot_grimu_probability_multi(n1_vec = c(8, 30, 70), n2_vec = c(9, 31, 71))