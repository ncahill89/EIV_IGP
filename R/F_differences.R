## Main Analysis function

run_rsl_test <- function() {
  
  message("Running ", ref_label1, " vs ", ref_label2)
  
  # -------------------------
  
  # RSL ESTIMATES
  
  # -------------------------
  
  ref_est1 <- read_rsl(
    file.path(results_dir, ref_dir1, "RSL_Estimates.csv"),
    ref_label1
  )
  
  ref_est2 <- read_rsl(
    file.path(results_dir, ref_dir2, "RSL_Estimates.csv"),
    ref_label2
  )
  
  rsl_data <- bind_rows(ref_est1, ref_est2)
  
  p_rsl <- plot_rsl(
    rsl_data,
    SL_est, SL_lwr, SL_upr,
    "Relative Sea Level (RSL) vs Time",
    "Relative Sea Level (m)"
  )
  
  # -------------------------
  
  # RSL DIFFERENCES
  
  # -------------------------
  
  ests_ref1 <- load_samples(
    file.path(model_dir, ref_dir1, "EstsandRates.rda"),
    "pred"
  )
  
  ests_ref2 <- load_samples(
    file.path(model_dir, ref_dir2, "EstsandRates.rda"),
    "pred"
  )
  
  diff_rsl <- compute_diff(
    ests_ref1,
    ests_ref2,
    ref_est1$year,
    paste(ref_label1, "-", ref_label2)
  )
  
  write_csv(
    diff_rsl,
    paste0("RSLDiffestimates_", ref_label1,"_",ref_label2, ".csv")
  )
  
  p_diff_rsl <- plot_diff(diff_rsl, "Difference (m)")
  
  # -------------------------
  
  # RSL RATES
  
  # -------------------------
  
  ref_rates1 <- read_rsl(
    file.path(results_dir, ref_dir1, "RSL_Rates.csv"),
    ref_label1
  )
  
  ref_rates2 <- read_rsl(
    file.path(results_dir, ref_dir2, "RSL_Rates.csv"),
    ref_label2
  )
  
  rate_data <- bind_rows(ref_rates1, ref_rates2)
  
  p_rates <- plot_rsl(
    rate_data,
    rate_est, rate_lwr, rate_upr,
    "RSL Rates vs Time",
    "Rates (mm/yr)"
  )
  
  # -------------------------
  
  # RATE DIFFERENCES
  
  # -------------------------
  
  rates_ref1 <- load_samples(
    file.path(model_dir, ref_dir1, "EstsandRates.rda"),
    "dydt"
  )
  
  rates_ref2 <- load_samples(
    file.path(model_dir, ref_dir2, "EstsandRates.rda"),
    "dydt"
  )
  
  diff_rates <- compute_diff(
    rates_ref1,
    rates_ref2,
    ref_est1$year,
    paste(ref_label1, "-", ref_label2)
  )
  
  write_csv(
    diff_rates,
    paste0("RateDiffestimates_",ref_label1,"_",ref_label2, ".csv")
  )
  
  p_diff_rates <- plot_diff(diff_rates, "Difference (mm/yr)")
  
  # -------------------------
  
  # DISPLAY
  
  # -------------------------
  
  print(ggarrange(p_rsl, p_diff_rsl, nrow = 1))
  print(ggarrange(p_rates, p_diff_rates, nrow = 1))
  
  # -------------------------
  
  # SAVE FIGURES (OPTIONAL)
  
  # -------------------------
  
  if (save_plots) {
    ggsave(
      paste0("RSL_",ref_label1,"_",ref_label2, ".pdf"),
      ggarrange(p_rsl, p_diff_rsl, nrow = 1),
      width = 12, height = 5
    )
    
    
    ggsave(
      paste0("Rates_", ref_label1,"_",ref_label2, ".pdf"),
      ggarrange(p_rates, p_diff_rates, nrow = 1),
      width = 12, height = 5
    )
    
    
  }
  
  invisible(list(
    rsl_table  = diff_rsl,
    rate_table = diff_rates
  ))
}

##### HELPER FUNCTIONS



read_rsl <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(Dataset = label)
}

plot_rsl <- function(data, y, ymin, ymax, title, ylab) {
  ggplot(data, aes(x = year, y = {{ y }}, color = Dataset, fill = Dataset)) +
    geom_line(size = 1) +
    geom_ribbon(
      aes(ymin = {{ ymin }}, ymax = {{ ymax }}),
      alpha = 0.2, color = NA
    ) +
    scale_x_reverse() +
    theme_minimal(base_size = 14) +
    labs(title = title, x = "Age", y = ylab) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
}


load_samples <- function(path, element) {
  load(path)
  EstsandRates[[element]]
}

compute_diff <- function(ref, test, years, label) {
  
  diff_mat <- ref - test
  
  as_tibble(diff_mat) %>%
    setNames(as.character(years)) %>%
    pivot_longer(everything(), names_to = "year") %>%
    mutate(year = as.numeric(year)) %>%
    group_by(year) %>%
    summarise(
      diff_est = mean(value),
      diff_lwr = diff_est - 2 * sd(value),
      diff_upr = diff_est + 2 * sd(value),
      .groups = "drop"
    ) %>%
    mutate(label = label)
}

plot_diff <- function(data, ylab) {
  ggplot(data, aes(x = year, y = diff_est)) +
    geom_line(size = 1) +
    geom_ribbon(
      aes(ymin = diff_lwr, ymax = diff_upr),
      alpha = 0.2
    ) +
    geom_hline(yintercept = 0) +
    scale_x_reverse() +
    theme_minimal(base_size = 14) +
    labs(x = "Age", y = ylab) +
    facet_wrap(~ label)
}


