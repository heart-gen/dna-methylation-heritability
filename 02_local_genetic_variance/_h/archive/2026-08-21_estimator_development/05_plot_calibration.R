#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required")

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "evaluation", "calibrated-evaluation-estimates.tsv"),
    output_dir = file.path(dirname(script_path), "..", "_m", "figures")
))
data <- read_tsv(cli$input)
data$n_label <- paste0("N = ", data$n)
data$architecture <- factor(data$architecture,
                            levels = c("sparse", "oligogenic", "polygenic"))

summary_data <- aggregate(
    cbind(h2_en_calibrated, h2_calibration_lower, h2_calibration_upper) ~
        n_label + architecture + true_h2,
    data = data,
    FUN = mean
)
plot <- ggplot2::ggplot(summary_data, ggplot2::aes(
    x = true_h2, y = h2_en_calibrated, color = architecture
)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
    ggplot2::geom_linerange(ggplot2::aes(
        ymin = h2_calibration_lower, ymax = h2_calibration_upper
    ), position = ggplot2::position_dodge(width = 0.018), alpha = 0.65) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.018), size = 2) +
    ggplot2::geom_line(ggplot2::aes(group = architecture), linewidth = 0.7) +
    ggplot2::facet_wrap(~ n_label) +
    ggplot2::scale_color_manual(values = c(
        sparse = "#0072B2", oligogenic = "#E69F00", polygenic = "#009E73"
    ), drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, 0.65), ylim = c(0, 0.65)) +
    ggplot2::labs(
        x = expression("Simulated local " * h^2),
        y = expression("Calibrated elastic-net local " * h^2),
        color = "Architecture",
        title = "Held-out calibration of local SNP-explained methylation variance",
        subtitle = "Points and bars are mean estimates and mean empirical calibration limits"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
        legend.position = "bottom",
        strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
        plot.title.position = "plot"
    )

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(file.path(cli$output_dir, "calibration-truth-versus-estimate.pdf"),
                plot, width = 8.2, height = 4.6, units = "in", device = cairo_pdf)
ggplot2::ggsave(file.path(cli$output_dir, "calibration-truth-versus-estimate.png"),
                plot, width = 8.2, height = 4.6, units = "in", dpi = 300)
cat("Wrote calibration figures to", normalizePath(cli$output_dir), "\n")
