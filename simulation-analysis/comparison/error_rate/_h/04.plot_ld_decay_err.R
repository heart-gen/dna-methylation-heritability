#### Plot error rate across LD decay levels for boosting_hybrid vs. joint_ridge ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
})

## Styling ##
method_labels <- c(
  "boosting_hybrid" = "Boosting-hybrid",
  "joint_ridge"     = "Joint-ridge"
)
metric_colors <- c(
  "Power"       = "#7B8C99",
  "Type 1 Error" = "#B35A4E",
  "Type 2 Error" = "#D4BFAA"
)
ld_levels <- c("0.8", "0.7", "0.6", "0.5")

save_plot <- function(p, fn, w, h, dpi = 600) {
  for (ext in c(".png", ".pdf")) {
    ggsave(file = paste0(fn, ext), plot = p, width = w, height = h, dpi = dpi,
           bg = "white")
  }
}

## Main ##
out_path <- here("simulation-analysis/comparison/error_rate/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

err_file <- file.path(out_path, "ld_decay_power-analysis.tsv")
if (!file.exists(err_file)) {
  stop("Error rate file not found: ", err_file,
       "\nRun 03.ld_decay_error_rate.py first.")
}

error_rate <- read_tsv(err_file, show_col_types = FALSE) |>
  mutate(
    ld_decay = factor(ld_decay, levels = ld_levels, ordered = TRUE),
    method   = factor(method_labels[method], levels = unname(method_labels))
  ) |>
  pivot_longer(
    cols      = c(power, type1_error, type2_error),
    names_to  = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = dplyr::recode(
      metric,
      "type1_error" = "Type 1 Error",
      "type2_error" = "Type 2 Error",
      "power"       = "Power"
    ),
    metric = factor(metric, levels = names(metric_colors))
  )

p <- ggplot(error_rate,
            aes(x = ld_decay, y = value,
                color = metric, group = interaction(metric, method),
                shape = method, linetype = method)) +
  geom_point(size = 4) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = metric_colors) +
  scale_shape_manual(values = c("Boosting-hybrid" = 16, "Joint-ridge" = 17)) +
  scale_linetype_manual(values = c("Boosting-hybrid" = "solid", "Joint-ridge" = "dashed")) +
  scale_x_discrete(labels = function(x) paste0("LD = ", x)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(
    x        = "LD decay level",
    y        = "Error measurement",
    color    = "Metric",
    shape    = "Method",
    linetype = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.8),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

plot_file <- file.path(out_path, "simulated_ld_decay_error_rate")
save_plot(p, plot_file, w = 8, h = 5)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
