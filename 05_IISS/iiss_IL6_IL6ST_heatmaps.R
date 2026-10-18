#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(scales)
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

plot_font_family <- Sys.getenv("PLOT_FONT_FAMILY", unset = "sans")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
script_dir <- dirname(script_path)
args <- parse_args(commandArgs(trailingOnly = TRUE))
xlsx_path <- args$input
if (is.null(xlsx_path)) {
  xlsx_path <- file.path(
    script_dir,
    "..",
    "..",
    "source_data",
    "IISS",
    "IISS_IL6_IL6ST_cell_type_summary.xlsx"
  )
}
base_dir <- args$`output-dir`
if (is.null(base_dir)) {
  base_dir <- file.path(script_dir, "output")
}
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

out_il6_pdf <- file.path(base_dir, "Fig3n_IISS_Il6_heatmap.pdf")
out_il6_png <- file.path(base_dir, "Fig3n_IISS_Il6_heatmap.png")
out_il6st_pdf <- file.path(base_dir, "SFig9a_IISS_Il6st_heatmap.pdf")
out_il6st_png <- file.path(base_dir, "SFig9a_IISS_Il6st_heatmap.png")
out_combined_pdf <- file.path(base_dir, "IISS_Il6_and_Il6st_heatmaps.pdf")

sheet_df <- read_excel(
  xlsx_path,
  sheet = "Cell-type proportions",
  col_names = FALSE,
  .name_repair = "minimal"
)
sample_names <- as.character(unlist(sheet_df[2, 2:5], use.names = FALSE))
sample_keys <- gsub(" ", "", sample_names, fixed = TRUE)
sample_keys <- sub("^Syn-control", "ctrl", sample_keys)
sample_keys <- sub("^Syn-let-7bKO", "ko", sample_keys)

sample_col_map <- stats::setNames(2:5, sample_keys)
keep_ctrl <- c("ctrl1", "ctrl2")
keep_ko <- c("ko1", "ko2")

row_values <- as.character(sheet_df[[1]])

fetch_group_values <- function(marker_prefix, celltype_keyword) {
  target <- paste0(marker_prefix, " ", celltype_keyword, " (%)")
  hit <- which(row_values == target)
  if (length(hit) != 1L) {
    stop(
      "Expected exactly one row labelled '", target, "'; found ", length(hit), ".",
      call. = FALSE
    )
  }
  vals <- as.numeric(unlist(sheet_df[hit, sample_col_map[c(keep_ctrl, keep_ko)]], use.names = FALSE))
  stats::setNames(vals, c(keep_ctrl, keep_ko))
}

round_outward <- function(x, step) {
  if (x >= 0) {
    ceiling(x / step) * step
  } else {
    floor(x / step) * step
  }
}

round_inward <- function(x, step) {
  if (x >= 0) {
    floor(x / step) * step
  } else {
    ceiling(x / step) * step
  }
}

trim_num <- function(x) {
  ifelse(
    abs(x) < 1e-8,
    "0",
    sub("\\.?0+$", "", sprintf("%.2f", x))
  )
}

build_plot_data <- function(marker_prefix) {
  cell_map <- c(
    ExtNeu = "excitatory neurons",
    InbNeu = "inhibitory neurons",
    Astro = "astrocytes",
    MiGlia = "microglia"
  )

  rows <- lapply(names(cell_map), function(short_label) {
    vals <- fetch_group_values(marker_prefix, cell_map[[short_label]])
    ctrl_mean <- mean(vals[keep_ctrl])
    ko_mean <- mean(vals[keep_ko])
    data.frame(
      Celltype_label_short = short_label,
      proportion_pct_syn_control = ctrl_mean,
      proportion_pct_syn_let7b_ko = ko_mean,
      fold_change_ko_vs_ctr = ko_mean / ctrl_mean,
      stringsAsFactors = FALSE
    )
  })

  plot_df <- do.call(rbind, rows)
  plot_df$Celltype_label_short <- factor(
    plot_df$Celltype_label_short,
    levels = c("ExtNeu", "InbNeu", "Astro", "MiGlia")
  )
  plot_df
}

build_heatmap_plot <- function(heatmap_grouped_wide, caption_text) {
  heatmap_rate_df <- rbind(
    data.frame(
      Row_label = "Syn-control",
      Celltype_label_short = heatmap_grouped_wide$Celltype_label_short,
      label_text = sprintf("%.2f", heatmap_grouped_wide$proportion_pct_syn_control),
      fill_hex = "#ECECEC",
      text_colour = "black",
      stringsAsFactors = FALSE
    ),
    data.frame(
      Row_label = "Syn-let-7b KO",
      Celltype_label_short = heatmap_grouped_wide$Celltype_label_short,
      label_text = sprintf("%.2f", heatmap_grouped_wide$proportion_pct_syn_let7b_ko),
      fill_hex = "#ECECEC",
      text_colour = "black",
      stringsAsFactors = FALSE
    )
  )

  fc_values <- heatmap_grouped_wide$fold_change_ko_vs_ctr
  fc_min <- min(fc_values, na.rm = TRUE)
  fc_max <- max(fc_values, na.rm = TRUE)
  fc_lower <- min(round_inward(fc_min, 0.1), 0.5)
  fc_upper <- max(round_outward(fc_max, 0.1), 1.3)
  fc_mid_lower <- round((fc_lower + 1) / 2, 2)
  fc_mid_upper <- round((1 + fc_upper) / 2, 2)

  fc_breaks <- sort(unique(c(
    fc_lower,
    fc_mid_lower,
    1,
    fc_mid_upper,
    fc_upper
  )))
  fc_break_labels <- trim_num(fc_breaks)

  fc_palette <- c("#2166AC", "#67A9CF", "#FFFFFF", "#F4B6B0", "#B2182B")
  fc_palette_values <- scales::rescale(
    c(fc_lower, fc_mid_lower, 1, fc_mid_upper, fc_upper),
    from = c(fc_lower, fc_upper)
  )
  fc_pal_fn <- scales::gradient_n_pal(
    colours = fc_palette,
    values = fc_palette_values,
    space = "Lab"
  )

  fc_fill_hex <- function(x) {
    fc_pal_fn(scales::rescale(x, from = c(fc_lower, fc_upper)))
  }

  heatmap_fc_df <- data.frame(
    Row_label = "Fold change",
    Celltype_label_short = heatmap_grouped_wide$Celltype_label_short,
    label_text = sprintf("%.2f", fc_values),
    fill_hex = fc_fill_hex(fc_values),
    text_colour = ifelse(abs(fc_values - 1) >= 0.45 * max(1 - fc_lower, fc_upper - 1), "white", "black"),
    stringsAsFactors = FALSE
  )

  heatmap_plot_df <- rbind(heatmap_rate_df, heatmap_fc_df)
  heatmap_plot_df$Row_label <- factor(
    heatmap_plot_df$Row_label,
    levels = rev(c("Syn-control", "Syn-let-7b KO", "Fold change"))
  )

  n_cols <- length(levels(heatmap_grouped_wide$Celltype_label_short))
  heatmap_right_edge <- n_cols + 0.5
  cbar_gap <- 0.64
  cbar_width <- 0.18
  cbar_x_left <- heatmap_right_edge + cbar_gap
  cbar_x_right <- cbar_x_left + cbar_width
  cbar_y_bot <- 0.62
  cbar_height <- 1.58
  cbar_y_top <- cbar_y_bot + cbar_height
  cbar_n_steps <- 200
  cbar_colors <- fc_pal_fn(seq(0, 1, length.out = cbar_n_steps))
  cbar_raster <- as.raster(matrix(rev(cbar_colors), ncol = 1))
  cbar_tick_y_pos <- cbar_y_bot + ((fc_breaks - fc_lower) / (fc_upper - fc_lower)) * cbar_height
  caption_x <- (n_cols + 1) / 2

  plot_obj <- ggplot(
    heatmap_plot_df,
    aes(x = Celltype_label_short, y = Row_label)
  ) +
    geom_tile(
      aes(fill = fill_hex),
      colour = "white", linewidth = 0.7, height = 0.95, width = 0.96
    ) +
    geom_text(
      aes(label = label_text, colour = text_colour),
      family = plot_font_family, fontface = "bold", size = 4.4
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_discrete(position = "top", expand = c(0, 0.05)) +
    scale_y_discrete(expand = c(0, 0.05)) +
    annotate(
      "text",
      x = seq_len(n_cols), y = 3.61,
      label = levels(heatmap_grouped_wide$Celltype_label_short),
      family = plot_font_family, fontface = "bold",
      size = 4.0, colour = "black"
    ) +
    annotation_raster(
      cbar_raster,
      xmin = cbar_x_left, xmax = cbar_x_right,
      ymin = cbar_y_bot, ymax = cbar_y_top,
      interpolate = TRUE
    ) +
    annotate(
      "segment",
      x = cbar_x_right, xend = cbar_x_right + 0.05,
      y = cbar_tick_y_pos, yend = cbar_tick_y_pos,
      colour = "black", linewidth = 0.3
    ) +
    annotate(
      "text",
      x = cbar_x_right + 0.09, y = cbar_tick_y_pos,
      label = fc_break_labels,
      hjust = 0, vjust = 0.5,
      size = 3.0, family = plot_font_family, colour = "black"
    ) +
    annotate(
      "text",
      x = caption_x, y = 0.16,
      label = caption_text,
      family = plot_font_family,
      size = 4.0, hjust = 0.5, vjust = 1
    ) +
    labs(x = NULL, y = NULL, title = NULL) +
    coord_cartesian(
      xlim = c(0.5, cbar_x_right + 0.55),
      ylim = c(0.30, 3.82),
      clip = "off"
    ) +
    theme_minimal(base_size = 13, base_family = plot_font_family) +
    theme(
      aspect.ratio = 2.15 / n_cols,
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        face = "bold", size = 12, colour = "black",
        margin = margin(r = 4)
      ),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 38, 25, 10)
    )

  list(
    plot = plot_obj,
    n_cols = n_cols
  )
}

il6_df <- build_plot_data("Il6-positive")
il6st_df <- build_plot_data("Il6st-positive")

plot_il6_info <- build_heatmap_plot(il6_df, "Sample-level Il6+ mean rate (%)")
plot_il6st_info <- build_heatmap_plot(il6st_df, "Sample-level Il6st+ mean rate (%)")
plot_il6 <- plot_il6_info$plot
plot_il6st <- plot_il6st_info$plot

write.csv(
  il6_df,
  file.path(base_dir, "Fig3n_IISS_Il6_heatmap_source_data.csv"),
  row.names = FALSE
)
write.csv(
  il6st_df,
  file.path(base_dir, "SFig9a_IISS_Il6st_heatmap_source_data.csv"),
  row.names = FALSE
)

reference_cols <- 5
reference_width_cm <- 15.0
plot_width_cm <- reference_width_cm * ((plot_il6_info$n_cols + 1.37) / (reference_cols + 1.37))

ggsave(out_il6_pdf, plot_il6, width = plot_width_cm, height = 8.4, units = "cm", device = grDevices::pdf, limitsize = FALSE)
ggsave(out_il6_png, plot_il6, width = plot_width_cm, height = 8.4, units = "cm", dpi = 300, limitsize = FALSE, bg = "white")
ggsave(out_il6st_pdf, plot_il6st, width = plot_width_cm, height = 8.4, units = "cm", device = grDevices::pdf, limitsize = FALSE)
ggsave(out_il6st_png, plot_il6st, width = plot_width_cm, height = 8.4, units = "cm", dpi = 300, limitsize = FALSE, bg = "white")

grDevices::pdf(out_combined_pdf, width = plot_width_cm / 2.54, height = 8.4 / 2.54, onefile = TRUE)
print(plot_il6)
print(plot_il6st)
grDevices::dev.off()
