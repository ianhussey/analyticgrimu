# --- original ---
simrank <- function(n1, n2, U_target, max_iter = 100000) {
  total_ranks <- 1:(n1 + n2)
  R1_target <- U_target + n1 * (n1 + 1) / 2
  
  for (i in 1:max_iter) {
    group1_ranks <- sample(total_ranks, n1, replace = FALSE)
    if (sum(group1_ranks) == R1_target) {
      group2_ranks <- setdiff(total_ranks, group1_ranks)
      return(list(
        group1_ranks = sort(group1_ranks), # sort moved from sample() line above to here, following Nico's pull request to the original repo
        group2_ranks = group2_ranks,
        U = U_target
      ))
    }
  }
  
  # Return NA instead of stopping with an error
  return(NA)
}

simrank_optimised <- function(n1, n2, U_target, max_iter = 1000) {
  # 1. Setup targets
  total_ranks <- 1:(n1 + n2)
  min_sum <- n1 * (n1 + 1) / 2
  max_sum <- min_sum + (n1 * n2)
  
  R1_target <- U_target + min_sum
  
  # 2. Safety: Fail immediately if target is mathematically impossible
  if (R1_target < min_sum || R1_target > max_sum) return(NA)
  
  # 3. Initialization: Random start
  group1_ranks <- sample(total_ranks, n1)
  current_sum <- sum(group1_ranks)
  
  # 4. Optimization Loop (Converging Swap)
  # Instead of resampling from scratch, we swap numbers to fix the error.
  for (i in 1:max_iter) {
    diff <- R1_target - current_sum
    
    if (diff == 0) {
      # Success!
      group1_ranks <- sort(group1_ranks)
      group2_ranks <- setdiff(total_ranks, group1_ranks)
      return(list(
        group1_ranks = group1_ranks,
        group2_ranks = group2_ranks,
        U = U_target
      ))
    }
    
    # We need to swap an element u (in group1) with v (in group2)
    # such that (v - u) moves the sum closer to the target.
    group2_ranks <- setdiff(total_ranks, group1_ranks)
    
    # Shuffle to maintain randomness in which specific ranks are chosen
    g1_candidates <- sample(group1_ranks)
    g2_candidates <- sample(group2_ranks)
    
    swapped <- FALSE
    
    # Greedy search for a helpful swap
    for (u in g1_candidates) {
      for (v in g2_candidates) {
        change <- v - u
        
        # If we need to increase sum (diff > 0), we need v > u.
        # If we need to decrease sum (diff < 0), we need v < u.
        # Check if this swap reduces the absolute error
        if (abs(diff - change) < abs(diff)) {
          # Perform Swap
          group1_ranks <- c(setdiff(group1_ranks, u), v)
          current_sum <- current_sum - u + v
          swapped <- TRUE
          break
        }
      }
      if (swapped) break
    }
    
    # If we loop through all pairs and find no improvement (local optimum), 
    # force a random swap to shake it up (rarely needed for rank sums)
    if (!swapped) {
      u <- sample(group1_ranks, 1)
      v <- sample(group2_ranks, 1)
      group1_ranks <- c(setdiff(group1_ranks, u), v)
      current_sum <- current_sum - u + v
    }
  }
  
  return(NA)
}

# --- simulation helper function ---
run_simulation_step <- function(k, n1, n2) {
  
  # --- LOGIC FIX: Use ceiling() for correct tie handling ---
  # If k is 70.5, we need to generate U=71 (ceiling) and subtract 0.5.
  # The original code used round(), which failed for even integers (e.g. 70.5 -> 70).
  U0 <- if (k %% 1 != 0) ceiling(k) else k
  
  # Run simulation (requires your simrank function to be loaded)
  #j <- simrank(n1, n2, U_target = U0)
  j <- simrank_optimised(n1, n2, U_target = U0)
  
  # --- FAILURE HANDLING ---
  if (is.na(j[1])) {
    return(tibble(p_exact = NA_real_, p_approx = NA_real_))
  }
  
  # --- TIE INJECTION LOGIC ---
  if (k %% 1 == 0) {
    # Integer K: Standard Test
    m1 <- wilcox.test(j$group1_ranks, j$group2_ranks)
    m2 <- wilcox.test(j$group1_ranks, j$group2_ranks, correct = FALSE, exact = FALSE)
  } else {
    # Adjust one value to simulate a fractional rank
    w <- which(diff(j$group2_ranks) == 2)[1]
    if (!is.na(w)) {
      mv <- j$group2_ranks[w]
      j$group2_ranks[w] <- mv + 0.5
      
      w2 <- which(j$group1_ranks == mv + 1)
      if (length(w2) > 0) {
        j$group1_ranks[w2] <- mv + 0.5
      }
    }
    
    m1 <- wilcox.test(j$group1_ranks, j$group2_ranks)
    m2 <- wilcox.test(j$group1_ranks, j$group2_ranks, correct = FALSE, exact = FALSE)
  }
  
  return(tibble(p_exact = m1$p.value, 
                p_approx = m2$p.value))
}


# --- 3. Main Analysis Function ---
grimu_analyze <- function(n1, n2, p_min = 0.01, p_max = 0.10, step = 0.005) {
  
  # --- Setup Calculation ---
  # 1. Expected stats (Ties and no ties)
  uexp <- (n1 * n2) / 2
  seu  <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
  # Standard Error with tie correction factor
  seu_tie <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12 - 
                    n1 * n2 * (7) / (12 * (n1 + n2) * (n1 + n2 - 1)))
  
  # 2. Range of p-values to U-values conversion
  # Determine Z-scores for the p-value bounds
  pvals_range <- c(p_min, p_max) 
  zscores <- qnorm(1 - pvals_range / 2) # Note: z[1] is high Z (low p), z[2] is low Z (high p)
  
  # 3. Calculate U Bounds
  # Upper U (corresponds to lowest p-value / highest Z)
  u_upper <- seu * zscores[1] + uexp 
  # Lower U (corresponds to highest p-value / lowest Z, using tie SE)
  u_lower <- seu_tie * zscores[2] + uexp 
  
  # 4. Generate Sequence of U-values to Test
  vals <- seq(floor(u_lower), ceiling(u_upper), by = 0.5)
  
  result_df <- tibble(U_Values = vals) %>%
    mutate(
      sim_results = map(U_Values, run_simulation_step, n1 = n1, n2 = n2)
    ) %>%
    unnest(sim_results) %>%
    mutate(
      Diff = abs(p_approx - p_exact),
      n1 = n1,
      n2 = n2
    )
  
  return(result_df)
}

plot_grimu_steps <- function(df_results) {
  
  # 1. Reshape data to long format for easy plotting of both methods
  plot_data <- df_results %>%
    pivot_longer(
      cols = c(p_exact, p_approx),
      names_to = "Method",
      values_to = "p_value"
    ) %>%
    mutate(
      Method = recode(Method, 
                      "p_exact" = "Exact Test (R default)", 
                      "p_approx" = "Asymptotic (SPSS style)")
    )
  
  # 2. Create Plot
  ggplot(plot_data, aes(x = U_Values, y = p_value, color = Method, shape = Method)) +
    # # Add points to show the discrete calculated values
    # geom_point(size = 1, alpha = 0.8) +
    # Add step lines to emphasize the "staircase" nature of the p-values
    geom_step(direction = "vh", linetype = "solid", alpha = 0.5) +
    # Scales and Labels
    scale_y_continuous(n.breaks = 10) +
    scale_x_continuous(n.breaks = 10) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "GRIM-U: Granularity of p-values",
      subtitle = paste0("Visualizing the discrete steps in p-values for N1=", 
                        unique(df_results$n1), ", N2=", unique(df_results$n2)),
      x = "Mann-Whitney U Statistic",
      y = "Calculated p-value",
      caption = "Note: 'Half-steps' in U represent tied ranks."
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank() # Remove minor grid to see steps clearly
    )
}