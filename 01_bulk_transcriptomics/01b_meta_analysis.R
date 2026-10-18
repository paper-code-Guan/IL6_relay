#!/usr/bin/env Rscript

# Stage 2 of the bulk DLPFC meta-analysis.
#
# Refits the random-effects model from the aggregate cohort table written by
# 01a_cohort_effect_sizes.R and renders Figure 1b. The pooled estimate is
# recomputed here rather than read from the table; the supplied pooled row is
# used only to verify that the recomputation reproduces the published values.
# Cohort estimates were derived from log10(CPM + 1)-transformed IL6 values.
#
# The horizontal scale is derived from the data and always includes the null
# value of 0, which is marked by a dashed reference line. Earlier versions fixed
# the scale to c(0.3, 1.2); that range excluded the null and silently clipped
# confidence limits once scale-correct standard errors were used.

suppressPackageStartupMessages({
  library(ggplot2)
  library(metafor)
})

parse_args <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be provided as --name value pairs.", call. = FALSE)
    }
    values[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_path))
args <- parse_args(commandArgs(trailingOnly = TRUE))

input_file <- args$input
if (is.null(input_file)) {
  input_file <- file.path(script_dir, "..", "..", "source_data",
                          "bulk_DLPFC_meta_analysis_input.csv")
}
output_dir <- args$`output-dir`
if (is.null(output_dir)) {
  output_dir <- file.path(script_dir, "output")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c(
  "cohort", "n_schizophrenia", "n_control",
  "effect_size_log2_expression_ratio", "standard_error",
  "ci_95_lower", "ci_95_upper", "weight_percent", "p_value"
)
input <- read.csv(input_file, check.names = FALSE)
missing_columns <- setdiff(required_columns, colnames(input))
if (length(missing_columns) > 0L) {
  stop(
    sprintf("Input is missing required columns: %s",
            paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}
# p_value_bound records whether a P value is exact or an upper bound, as for the
# LIBD cohort reported as P < 2e-11. Older tables without the column are treated
# as exact.
if (!"p_value_bound" %in% colnames(input)) {
  input$p_value_bound <- "exact"
}

pooled_row <- input[input$cohort == "Random-effects pooled", , drop = FALSE]
study_data <- input[input$cohort != "Random-effects pooled", , drop = FALSE]
expected_cohorts <- c("CMC", "HBCC", "LIBD")
if (!identical(study_data$cohort, expected_cohorts)) {
  stop(
    sprintf("Cohorts must be ordered as: %s",
            paste(expected_cohorts, collapse = ", ")),
    call. = FALSE
  )
}

fit <- rma(
  yi = effect_size_log2_expression_ratio,
  sei = standard_error,
  method = "REML",
  data = study_data
)

meta_row <- data.frame(
  cohort = "Meta",
  n_schizophrenia = sum(study_data$n_schizophrenia),
  n_control = sum(study_data$n_control),
  effect_size_log2_expression_ratio = as.numeric(fit$b),
  standard_error = fit$se,
  ci_95_lower = fit$ci.lb,
  ci_95_upper = fit$ci.ub,
  weight_percent = 100,
  p_value = fit$pval,
  p_value_bound = "exact",
  stringsAsFactors = FALSE
)

if (nrow(pooled_row) == 1L) {
  fields <- c(
    "effect_size_log2_expression_ratio", "standard_error",
    "ci_95_lower", "ci_95_upper", "p_value"
  )
  observed <- as.numeric(meta_row[1, fields])
  expected <- as.numeric(pooled_row[1, fields])
  if (any(abs(observed - expected) > 5e-6)) {
    stop("Recomputed pooled estimates do not match the supplied reference row.",
         call. = FALSE)
  }
}

# Cohort confidence limits are taken from the input table as released. The check
# below only guards against a corrupted or mismatched table; it does not
# overwrite the supplied values.
ci_lower_check <- study_data$effect_size_log2_expression_ratio -
  1.96 * study_data$standard_error
ci_upper_check <- study_data$effect_size_log2_expression_ratio +
  1.96 * study_data$standard_error
if (any(abs(study_data$ci_95_lower - ci_lower_check) > 1e-4) ||
    any(abs(study_data$ci_95_upper - ci_upper_check) > 1e-4)) {
  stop(
    "Supplied cohort confidence limits are inconsistent with the effect sizes and standard errors.",
    call. = FALSE
  )
}

plot_data <- rbind(study_data[, colnames(meta_row)], meta_row)
plot_data$cohort <- factor(
  plot_data$cohort,
  levels = c("Meta", "LIBD", "HBCC", "CMC")
)
# Only the pooled estimate carries a significance label. Cohort P values are reported in
# Supplementary Figure 1a-c.
plot_data$significance <- as.character(cut(
  plot_data$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", ""),
  right = FALSE
))
plot_data$significance[as.character(plot_data$cohort) != "Meta"] <- ""

stats <- data.frame(
  cohort = as.character(plot_data$cohort),
  n_schizophrenia = plot_data$n_schizophrenia,
  n_control = plot_data$n_control,
  estimate = plot_data$effect_size_log2_expression_ratio,
  standard_error = plot_data$standard_error,
  ci_95_lower = plot_data$ci_95_lower,
  ci_95_upper = plot_data$ci_95_upper,
  p_value = plot_data$p_value,
  p_value_bound = plot_data$p_value_bound,
  weight_percent = plot_data$weight_percent,
  tau_squared = c(rep(NA_real_, nrow(plot_data) - 1L), fit$tau2),
  i_squared_percent = c(rep(NA_real_, nrow(plot_data) - 1L), fit$I2),
  q_statistic = c(rep(NA_real_, nrow(plot_data) - 1L), fit$QE),
  q_p_value = c(rep(NA_real_, nrow(plot_data) - 1L), fit$QEp)
)
write.csv(
  stats,
  file.path(output_dir, "Fig1b_meta_analysis_statistics.csv"),
  row.names = FALSE
)

font_family <- Sys.getenv("PLOT_FONT_FAMILY", unset = "Helvetica")

# Horizontal scale: span every confidence limit, always include the null value
# of 0, and leave room on the right for the significance labels.
axis_lower <- min(0, plot_data$ci_95_lower)
axis_upper <- max(0, plot_data$ci_95_upper)
axis_span <- axis_upper - axis_lower
if (!is.finite(axis_span) || axis_span <= 0) {
  stop("Confidence limits do not define a usable horizontal scale.",
       call. = FALSE)
}
axis_limits <- c(axis_lower - 0.06 * axis_span, axis_upper + 0.14 * axis_span)
axis_breaks <- pretty(c(axis_lower, axis_upper), n = 6)
axis_breaks <- axis_breaks[
  axis_breaks >= axis_limits[1] & axis_breaks <= axis_limits[2]
]
label_offset <- 0.03 * diff(axis_limits)

figure <- ggplot(
  plot_data,
  aes(x = effect_size_log2_expression_ratio, y = cohort)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  geom_errorbar(
    aes(xmin = ci_95_lower, xmax = ci_95_upper),
    width = 0.06,
    linewidth = 0.8,
    colour = "#C00000",
    orientation = "y"
  ) +
  geom_point(size = 1.7, shape = 21, fill = "black", colour = "black") +
  geom_text(
    aes(x = ci_95_upper + label_offset, label = significance),
    hjust = 0,
    family = font_family,
    fontface = "bold",
    size = 5.2
  ) +
  # Heterogeneity statistics sit in a reserved band below the pooled row. Placing
  # them inside the plotting area risked collision with the significance labels,
  # whose horizontal positions follow the confidence limits and therefore move
  # with the data.
  annotate(
    "text",
    x = axis_limits[2] - 0.02 * diff(axis_limits),
    y = 0.34,
    label = sprintf("I^2 == %.2f*'%%'~~~tau^2 == %.3f", fit$I2, fit$tau2),
    hjust = 1,
    fontface = "bold",
    size = 4.1,
    parse = TRUE
  ) +
  scale_x_continuous(
    limits = axis_limits,
    breaks = axis_breaks,
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = ggplot2::expansion(add = c(0.95, 0.6))) +
  labs(x = expression(log[2]~"expression ratio (95%CI)"), y = NULL) +
  theme_classic(base_family = font_family, base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.text = element_text(colour = "black", face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    plot.margin = margin(3, 8, 3, 3)
  )

ggsave(
  file.path(output_dir, "Fig1b_meta_analysis.pdf"),
  plot = figure,
  width = 3.25,
  height = 2.4,
  units = "in",
  device = grDevices::pdf,
  useDingbats = FALSE
)

message(sprintf(
  paste(
    "Meta-analysis complete: estimate %.7f, 95%% CI %.7f to %.7f,",
    "P = %.7g, I-squared = %.2f%%, tau-squared = %.4f"
  ),
  as.numeric(fit$b), fit$ci.lb, fit$ci.ub, fit$pval, fit$I2, fit$tau2
))
