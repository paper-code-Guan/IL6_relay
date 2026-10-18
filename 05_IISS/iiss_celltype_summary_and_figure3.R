#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(Seurat)
  library(tibble)
  library(tidyr)
})

set.seed(123)
requested_font <- Sys.getenv("PLOT_FONT_FAMILY", unset = "")
plot_font_family <- if (nzchar(requested_font)) {
  requested_font
} else if ("ArialMT" %in% names(grDevices::pdfFonts())) {
  "ArialMT"
} else {
  "sans"
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
data_dir <- Sys.getenv(
  "IISS_SOURCE_DATA_DIR",
  unset = file.path(script_dir, "..", "..", "source_data", "IISS")
)
output_dir <- Sys.getenv(
  "IISS_OUTPUT_DIR",
  unset = file.path(script_dir, "output")
)
data_dir <- normalizePath(data_dir, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Bubble"), recursive = TRUE, showWarnings = FALSE)
dir.create(
  file.path(output_dir, "Composition"),
  recursive = TRUE,
  showWarnings = FALSE
)

sample_manifest <- tibble::tribble(
  ~file,       ~Group, ~Sample,
  "con2.txt",  "Con",  "Con2",
  "con10.txt", "Con",  "Con10",
  "ko3.txt",   "Ko",   "Ko3",
  "ko11.txt",  "Ko",   "Ko11"
)

read_iiss_sample <- function(file, Group, Sample) {
  input_path <- file.path(data_dir, file)
  if (!file.exists(input_path)) {
    stop("Missing IISS source-data file: ", input_path, call. = FALSE)
  }
  sample_data <- readr::read_tsv(
    input_path,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    as.data.frame()
  required_columns <- c("cell.name", "celltype")
  missing_columns <- setdiff(required_columns, colnames(sample_data))
  if (length(missing_columns)) {
    stop(
      file,
      " is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  sample_data$Group <- Group
  sample_data$Sample <- Sample
  sample_data$cell.name <- paste0(sample_data$cell.name, "_", Sample)
  sample_data
}

all_data <- purrr::pmap_dfr(sample_manifest, read_iiss_sample) |>
  as.data.frame()

marker_genes <- c(
  "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Slc1a2", "S100b",
  "Aif1", "Tmem119", "Cx3cr1", "Flt1", "Cldn5"
)
celltype_map <- c(
  "ExciNeuN" = "ExtNeu",
  "InhiNeuN" = "InbNeu",
  "Astrocyte" = "Astro",
  "Microglia" = "MiGlia",
  "Endo" = "Endoth"
)
celltype_levels <- unname(celltype_map)
celltype_colors <- c(
  "ExtNeu" = "#377eb8",
  "InbNeu" = "#a6cee3",
  "Astro" = "#e41a1c",
  "MiGlia" = "#4daf4a",
  "Endoth" = "#984ea3"
)

for (source_label in names(celltype_map)) {
  all_data$celltype[grepl(source_label, all_data$celltype)] <-
    celltype_map[[source_label]]
}
all_data <- all_data |>
  filter(celltype %in% celltype_levels) |>
  tidyr::drop_na()

metadata_columns <- c("cell.name", "celltype", "Group", "Sample")
gene_columns <- setdiff(colnames(all_data), metadata_columns)
if (!length(gene_columns)) {
  stop("No assigned RCP-count columns were detected.", call. = FALSE)
}
if (anyDuplicated(all_data$cell.name)) {
  stop(
    "Cell identifiers are not unique after sample suffixes were added.",
    call. = FALSE
  )
}
missing_markers <- setdiff(marker_genes, gene_columns)
if (length(missing_markers)) {
  stop(
    "Missing marker columns: ",
    paste(missing_markers, collapse = ", "),
    call. = FALSE
  )
}

all_expr <- t(as.matrix(all_data[, gene_columns, drop = FALSE]))
storage.mode(all_expr) <- "numeric"
all_meta <- all_data[, metadata_columns, drop = FALSE]
rownames(all_meta) <- all_meta$cell.name
all_meta$cell.name <- NULL
colnames(all_expr) <- rownames(all_meta)
detected_cells <- colSums(all_expr) > 0
all_expr <- all_expr[, detected_cells, drop = FALSE]
all_meta <- all_meta[detected_cells, , drop = FALSE]

seurat_object <- CreateSeuratObject(counts = all_expr, meta.data = all_meta)
seurat_object <- SCTransform(
  seurat_object,
  verbose = FALSE,
  vars.to.regress = "Sample"
)

cell_order <- unlist(
  lapply(
    celltype_levels,
    function(celltype) rownames(all_meta)[all_meta$celltype == celltype]
  ),
  use.names = FALSE
)
scaled_marker_matrix <- as.matrix(
  seurat_object@assays$SCT$counts[marker_genes, cell_order]
)
scaled_marker_matrix <- t(scale(t(scaled_marker_matrix)))

dotplot_data <- as.data.frame(scaled_marker_matrix) |>
  tibble::rownames_to_column("marker") |>
  tidyr::pivot_longer(-marker, names_to = "cell", values_to = "expr") |>
  left_join(
    all_meta |>
      tibble::rownames_to_column("cell") |>
      select(cell, celltype),
    by = "cell"
  ) |>
  group_by(celltype, marker) |>
  summarise(
    avg_expr = mean(expr, na.rm = TRUE),
    pct_expr = sum(expr > 0, na.rm = TRUE) / n() * 100,
    .groups = "drop"
  )
dotplot_data$celltype <- factor(dotplot_data$celltype, levels = celltype_levels)
dotplot_data$marker <- factor(dotplot_data$marker, levels = marker_genes)

font_size <- 50
line_width <- 2.5
dotplot <- ggplot(dotplot_data, aes(x = celltype, y = marker)) +
  geom_point(aes(color = avg_expr, size = pct_expr)) +
  scale_color_gradient2(
    low = "#7FB2D5",
    mid = "white",
    high = "#F47F72",
    midpoint = 0,
    limits = range(dotplot_data$avg_expr),
    name = "Normalized expression",
    breaks = c(-1, 0, 3)
  ) +
  scale_size(range = c(1, 20)) +
  theme_classic(base_family = plot_font_family) +
  theme(
    axis.line = element_line(colour = "black", linewidth = line_width),
    axis.text.x = element_text(
      size = font_size,
      family = plot_font_family,
      colour = "black",
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      size = font_size,
      face = "bold",
      family = plot_font_family,
      colour = "black",
      margin = margin(r = 10)
    ),
    axis.title = element_blank(),
    legend.title = element_text(
      size = font_size - 2,
      face = "bold",
      family = plot_font_family,
      colour = "black",
      margin = margin(b = 20)
    ),
    legend.text = element_text(
      size = font_size - 4,
      family = plot_font_family,
      face = "bold",
      colour = "black"
    ),
    legend.position = "right"
  ) +
  guides(
    color = guide_colorbar(title = "Expression", order = 2),
    size = guide_legend(title = "Percentage", order = 1)
  )

jpeg(
  filename = file.path(output_dir, "Bubble", "All_bubble_v7.jpg"),
  res = 300,
  width = 6050,
  height = 4000
)
print(dotplot)
dev.off()
ggsave(
  filename = file.path(output_dir, "Fig3m_IISS_marker_dotplot.pdf"),
  plot = dotplot,
  width = 6050 / 300,
  height = 4000 / 300,
  units = "in",
  device = grDevices::pdf,
  limitsize = FALSE,
  useDingbats = FALSE
)
write.csv(
  dotplot_data,
  file.path(output_dir, "Fig3m_IISS_marker_dotplot_source_data.csv"),
  row.names = FALSE
)

round_percentages_to_total <- function(proportions, digits = 1) {
  scale <- 10^digits
  raw_units <- proportions * 100 * scale
  rounded_units <- floor(raw_units)
  units_to_assign <- as.integer(round(100 * scale - sum(rounded_units)))
  if (units_to_assign > 0) {
    remainder_order <- order(raw_units - rounded_units, decreasing = TRUE)
    rounded_units[remainder_order[seq_len(units_to_assign)]] <-
      rounded_units[remainder_order[seq_len(units_to_assign)]] + 1
  }
  rounded_units / scale
}

composition_data <- all_meta |>
  tibble::rownames_to_column("cell") |>
  count(Group, celltype, name = "Count") |>
  tidyr::complete(
    Group = c("Con", "Ko"),
    celltype = celltype_levels,
    fill = list(Count = 0)
  ) |>
  group_by(Group) |>
  mutate(
    Proportion = Count / sum(Count),
    Percent = round_percentages_to_total(Proportion, digits = 1)
  ) |>
  ungroup()

for (group_id in c("Con", "Ko")) {
  group_data <- composition_data |>
    filter(Group == group_id) |>
    arrange(match(celltype, celltype_levels)) |>
    mutate(
      CellType = factor(celltype, levels = celltype_levels),
      ypos = 1 - (cumsum(Proportion) - 0.5 * Proportion),
      Label = Percent,
      LabelPos = ifelse(Percent < 4, "outer", "inner")
    )

  composition_plot <- ggplot(
    group_data,
    aes(x = 2, y = Proportion, fill = CellType)
  ) +
    geom_bar(
      stat = "identity",
      color = "black",
      linewidth = 1,
      show.legend = FALSE
    ) +
    coord_polar(theta = "y", start = 0) +
    xlim(1, 3) +
    theme_void() +
    scale_fill_manual(
      values = celltype_colors,
      breaks = celltype_levels,
      drop = FALSE
    ) +
    geom_text(
      data = subset(group_data, LabelPos == "inner"),
      aes(y = ypos, label = Label),
      color = "white",
      family = plot_font_family,
      size = 8,
      fontface = "bold"
    ) +
    geom_text(
      data = subset(group_data, LabelPos == "outer"),
      aes(x = 2.9, y = ypos, label = Label),
      color = "black",
      family = plot_font_family,
      size = 8,
      fontface = "bold"
    ) +
    geom_segment(
      data = subset(group_data, LabelPos == "outer"),
      aes(x = 2, xend = 2.7, y = ypos, yend = ypos),
      color = "black",
      linewidth = 1
    )

  output_stem <- if (group_id == "Con") {
    "Syn-control"
  } else {
    "Syn-let-7b_KO"
  }
  jpeg(
    filename = file.path(
      output_dir,
      "Composition",
      paste0(output_stem, "_composition.jpg")
    ),
    res = 300,
    width = 1500,
    height = 1500
  )
  print(composition_plot)
  dev.off()
  ggsave(
    filename = file.path(
      output_dir,
      paste0("Fig3l_", output_stem, "_composition.pdf")
    ),
    plot = composition_plot,
    width = 5,
    height = 5,
    units = "in",
    device = grDevices::pdf,
    limitsize = FALSE,
    useDingbats = FALSE
  )
  write.csv(
    group_data |>
      select(Group, CellType, Count, Proportion, Percent),
    file.path(
      output_dir,
      paste0("Fig3l_", output_stem, "_composition_source_data.csv")
    ),
    row.names = FALSE
  )
}

message(
  sprintf(
    "IISS summaries complete: %d retained cells across %d cell types.",
    ncol(seurat_object),
    length(celltype_levels)
  )
)
