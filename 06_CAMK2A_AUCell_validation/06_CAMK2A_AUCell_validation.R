# =========================================================
# Batiuk human step6: CAMK2A AUCell threshold stability validation
#
# Purpose:
#   Positive-control validation for the single-gene AUCell positivity rule
#   used in the Batiuk human IL6 analysis. CAMK2A is used as a canonical
#   excitatory-neuron marker. The script compares:
#     1) raw count positivity: CAMK2A count > 0
#     2) AUCell-based positivity: same single-gene rule used for IL6 in
#        step4/step5, i.e. rank positivity at the extended top window OR
#        raw-count detection protection.
#
# Expected input:
#   1) Annotated human Seurat object specified by SEURAT_RDS_PATH.
#   2) Formal AUCell output directory specified by BATIUK_RESULT_ROOT.
#
# Main outputs:
#   step6_final_batiuk_<timestamp>/
#     figures/
#       camk2a_count_positive_composition_donut_<timestamp>.pdf
#       camk2a_aucell_positive_composition_donut_<timestamp>.pdf
#     tables/source_data/
#       camk2a_count_positive_composition_source_<timestamp>.csv
#       camk2a_aucell_positive_composition_source_<timestamp>.csv
#       camk2a_concordance_source_<timestamp>.csv
#     tables/summary/
#       camk2a_composition_summary_<timestamp>.csv
#       camk2a_concordance_summary_<timestamp>.csv
#     logs/
#       step6_final_batiuk_run_info_<timestamp>.txt
#     rds/
#       step6_final_batiuk_objects_<timestamp>.rds
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(123)

required_aucell_version <- "1.30.0"
if (utils::compareVersion(as.character(packageVersion("AUCell")), required_aucell_version) < 0) {
  stop("AUCell >= ", required_aucell_version, " is required; installed version is ", packageVersion("AUCell"))
}

# ---------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------
human_rds <- Sys.getenv("SEURAT_RDS_PATH", unset = "")
result_root <- Sys.getenv("BATIUK_RESULT_ROOT", unset = "")
if (!nzchar(human_rds)) {
  stop("SEURAT_RDS_PATH must be set for Batiuk step6.", call. = FALSE)
}
if (!nzchar(result_root)) {
  stop("BATIUK_RESULT_ROOT must be set for Batiuk step6.", call. = FALSE)
}

tab_src_new <- file.path(result_root, "step4_formal", "tables")
tab_src_old <- file.path(result_root, "tables")
tab_src <- if (dir.exists(tab_src_new)) tab_src_new else tab_src_old

stopifnot("Human Seurat RDS not found." = file.exists(human_rds))
stopifnot("AUCell table directory not found." = dir.exists(tab_src))

timestamp_tag <- format(Sys.time(), "%Y%m%d%H%M")
out_root <- file.path(result_root, paste0("step6_final_batiuk_", timestamp_tag))
fig_dir <- file.path(out_root, "figures")
tab_dir <- file.path(out_root, "tables")
source_tab_dir <- file.path(tab_dir, "source_data")
summary_tab_dir <- file.path(tab_dir, "summary")
log_dir <- file.path(out_root, "logs")
rds_dir <- file.path(out_root, "rds")

for (dir_i in c(out_root, fig_dir, tab_dir, source_tab_dir, summary_tab_dir, log_dir, rds_dir)) {
  dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
}

# ---------------------------------------------------------
# 1. Labels, colors, and helpers
# ---------------------------------------------------------
display_levels <- c("CTR", "SCZ")
target_gene <- "CAMK2A"
expected_celltype <- "ExcitNeu"
top_percent_use <- 0.05
top_percent_use_single <- 0.10
export_cell_level_source <- TRUE
plot_font_family <- "sans"

ct_order <- c("Ast", "ExcitNeu", "InhibNeu", "OL", "Others")
ct_full_map <- c(
  "Ast" = "Astro",
  "ExcitNeu" = "Excitatory neurons",
  "InhibNeu" = "Inhibitory neurons",
  "OL" = "Oligodendrocytes",
  "Others" = "Others"
)
ct_fig_map <- c(
  "Ast" = "Astro",
  "ExcitNeu" = "ExtNeu",
  "InhibNeu" = "InbNeu",
  "OL" = "Oligo",
  "Others" = "Others"
)
cell_colors_fig <- c(
  "ExtNeu" = "#FCAC61",
  "InbNeu" = "#3E8FC6",
  "Astro" = "#776EA5",
  "Oligo" = "#9F0340",
  "Others" = "#BDBDBD"
)
fig_levels <- unname(ct_fig_map[ct_order])

canonicalize_human_ct <- function(v) {
  x <- trimws(as.character(v))
  xl <- tolower(x)
  out <- rep(NA_character_, length(x))
  out[xl %in% c("ast", "astro", "astrocyte", "astrocytes")] <- "Ast"
  out[xl %in% c("excitneu", "exn", "glutamatergic neuron", "glutamatergic neurons",
                "excitatory neuron", "excitatory neurons")] <- "ExcitNeu"
  out[xl %in% c("inhibneu", "inn", "gabaergic neuron", "gabaergic neurons",
                "inhibitory neuron", "inhibitory neurons", "interneuron", "interneurons")] <- "InhibNeu"
  out[xl %in% c("mic", "micro", "microglia", "microglial")] <- "Mic"
  out[xl %in% c("ol", "oligo", "oligodendrocyte", "oligodendrocytes")] <- "OL"
  out[xl %in% c("other", "others", "other_mix", "other cells", "other cell",
                "other_cell", "othercell")] <- "Others"
  already_standard <- x %in% c("Ast", "ExcitNeu", "InhibNeu", "Mic", "OL", "Others")
  out[is.na(out) & already_standard] <- x[is.na(out) & already_standard]
  out
}

pretty_ct <- function(x) unname(ct_full_map[as.character(x)])
pretty_ct_fig <- function(x) unname(ct_fig_map[as.character(x)])

normalize_group <- function(x) {
  g <- trimws(as.character(x))
  g[g %in% c("Ctr", "CTR", "Ctrl", "CTRL", "Control", "CONTROL", "control", "Con", "CON")] <- "CTR"
  g[g %in% c("Scz", "SCZ", "Case", "CASE", "case", "Patient", "PATIENT", "patient")] <- "SCZ"
  g
}

infer_group_from_metadata <- function(md) {
  candidate_cols <- c("Diagnosis", "Group", "group", "Group_use", "group_use", "samp.Diagnosis")
  for (col in candidate_cols) {
    if (col %in% colnames(md)) {
      grp <- normalize_group(md[[col]])
      if (!all(is.na(grp))) {
        return(grp)
      }
    }
  }

  fallback_cols <- c("orig.ident", "Patient", "sample", "Sample")
  for (col in fallback_cols) {
    if (col %in% colnames(md)) {
      x <- trimws(as.character(md[[col]]))
      grp <- ifelse(
        grepl("(^|[_-])(ctr|ctrl|control|con)($|[_-])", x, ignore.case = TRUE),
        "CTR",
        ifelse(
          grepl("(^|[_-])(scz|case|patient)($|[_-])", x, ignore.case = TRUE),
          "SCZ",
          NA_character_
        )
      )
      if (!all(is.na(grp))) {
        return(grp)
      }
    }
  }

  rep(NA_character_, nrow(md))
}

resolve_gene_ci <- function(gene_symbol, pool_rows) {
  idx <- which(tolower(pool_rows) == tolower(gene_symbol))
  if (!length(idx)) stop("Gene not found: ", gene_symbol)
  pool_rows[idx[1]]
}

resolve_features_ci <- function(symbols, pool_rows) {
  idx <- match(tolower(symbols), tolower(pool_rows))
  res <- ifelse(is.na(idx), NA_character_, pool_rows[idx])
  list(
    found = unique(stats::na.omit(res)),
    missing = unique(symbols[is.na(idx)]),
    mapping = data.frame(requested = symbols, resolved = unname(res), stringsAsFactors = FALSE)
  )
}

seurat_get_layer <- function(object, assay = "RNA", layer = "counts") {
  out <- tryCatch(
    Seurat::GetAssayData(object = object, assay = assay, layer = layer),
    error = function(e) e
  )
  if (!inherits(out, "error")) return(out)

  assay_obj <- object[[assay]]
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      "LayerData" %in% getNamespaceExports("SeuratObject")) {
    layer_out <- tryCatch(
      SeuratObject::LayerData(object = object, assay = assay, layer = layer),
      error = function(e) e
    )
    if (!inherits(layer_out, "error")) return(layer_out)

    if (inherits(assay_obj, "Assay5")) {
      layer_names <- SeuratObject::Layers(assay_obj)
      matching_layers <- grep(paste0("^", layer), layer_names, value = TRUE)
      if (!length(matching_layers)) {
        stop("Layer '", layer, "' not found in assay ", assay, call. = FALSE)
      }
      layer_mats <- lapply(matching_layers, function(layer_name) {
        SeuratObject::LayerData(object = object, assay = assay, layer = layer_name)
      })
      common_features <- Reduce(intersect, lapply(layer_mats, rownames))
      if (!length(common_features)) {
        stop("No common features across matching layers: ", paste(matching_layers, collapse = ", "))
      }
      layer_mats <- lapply(layer_mats, function(m) m[common_features, , drop = FALSE])
      return(do.call(Matrix::cbind2, layer_mats))
    }
  }

  slot_out <- tryCatch(
    Seurat::GetAssayData(object = object, assay = assay, slot = layer),
    error = function(e) e
  )
  if (!inherits(slot_out, "error")) return(slot_out)
  stop("Unable to extract assay layer/slot '", layer, "' from assay '", assay, "'.")
}

save_plot_pdf <- function(plot_obj, filename, width, height, units = "mm") {
  ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    units = units,
    device = grDevices::pdf,
    limitsize = FALSE
  )
}

compute_percent_vector <- function(x) {
  total_x <- sum(x, na.rm = TRUE)
  if (!is.finite(total_x) || total_x <= 0) {
    return(rep(0, length(x)))
  }
  100 * x / total_x
}

# ---------------------------------------------------------
# 2. Load and harmonize Seurat object metadata
# ---------------------------------------------------------
message("Loading Seurat object: ", human_rds)
obj <- readRDS(human_rds)
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
if (!"Celltype_use" %in% colnames(md)) {
  if ("Celltypes" %in% colnames(md)) {
    obj$Celltypes_raw <- as.character(md$Celltypes)
    obj$Celltype_use <- factor(canonicalize_human_ct(md$Celltypes), levels = ct_order)
  } else if ("Celltypes_raw" %in% colnames(md)) {
    obj$Celltype_use <- factor(canonicalize_human_ct(md$Celltypes_raw), levels = ct_order)
  } else {
    stop("No Celltype_use, Celltypes, or Celltypes_raw column found in Seurat metadata.")
  }
} else {
  obj$Celltype_use <- factor(canonicalize_human_ct(md$Celltype_use), levels = ct_order)
}

if (!"Celltypes_raw" %in% colnames(obj@meta.data)) {
  obj$Celltypes_raw <- as.character(obj$Celltype_use)
}
obj$Celltype_label <- factor(pretty_ct(obj$Celltype_use), levels = unname(ct_full_map[ct_order]))
obj$Celltype_fig_label <- factor(pretty_ct_fig(obj$Celltype_use), levels = fig_levels)

if ("Group" %in% colnames(obj@meta.data)) {
  obj$Group <- factor(normalize_group(obj$Group), levels = display_levels)
} else {
  inferred_group <- infer_group_from_metadata(obj@meta.data)
  if (all(is.na(inferred_group))) {
    stop("Unable to resolve Group from Seurat metadata. Checked Diagnosis/Group/group/Group_use/samp.Diagnosis and sample-like columns.")
  }
  obj$Group <- factor(inferred_group, levels = display_levels)
}
if (!"Patient" %in% colnames(obj@meta.data)) {
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    obj$Patient <- as.character(obj$orig.ident)
  } else {
    stop("Neither Patient nor orig.ident is available in Seurat metadata.")
  }
}

keep_cells <- !is.na(obj$Celltype_use)
if (!all(keep_cells)) {
  obj <- subset(obj, cells = colnames(obj)[keep_cells])
}

cts_present <- ct_order[ct_order %in% as.character(unique(obj$Celltype_use))]
if (!expected_celltype %in% cts_present) {
  warning("Expected cell type ", expected_celltype, " is not present in this object.")
}

# ---------------------------------------------------------
# 3. CAMK2A single-gene AUCell scoring using the IL6 single-gene rule
# ---------------------------------------------------------
compute_single_gene_score_df <- function(seu, ct_name, gene_symbol,
                                         top_percent_default = 0.05,
                                         top_percent_single = 0.10) {
  cells_ct <- rownames(seu@meta.data)[as.character(seu$Celltype_use) == ct_name]
  if (!length(cells_ct)) stop("No cells found for CT: ", ct_name)

  seu_ct <- subset(seu, cells = cells_ct)
  meta_ct <- seu_ct@meta.data
  counts_ct <- seurat_get_layer(seu_ct, assay = "RNA", layer = "counts")

  gene_map <- resolve_features_ci(gene_symbol, rownames(counts_ct))
  genes_found <- gene_map$found
  if (length(genes_found) != 1L) {
    stop("Target single gene not found uniquely: ", gene_symbol,
         "\nFound: ", paste(genes_found, collapse = ", "))
  }

  if (!inherits(counts_ct, "dgCMatrix")) {
    counts_ct <- as(counts_ct, "CsparseMatrix")
  }

  K_top_default <- max(1L, ceiling(top_percent_default * nrow(counts_ct)))
  K_top_single <- max(1L, ceiling(top_percent_single * nrow(counts_ct)))

  rankings <- AUCell::AUCell_buildRankings(
    counts_ct,
    plotStats = FALSE,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam()
  )

  gene_set_named <- setNames(list(genes_found), gene_symbol)
  auc_res_default <- AUCell::AUCell_calcAUC(
    geneSets = gene_set_named,
    rankings = rankings,
    aucMaxRank = K_top_default,
    verbose = FALSE
  )
  auc_res_single <- AUCell::AUCell_calcAUC(
    geneSets = gene_set_named,
    rankings = rankings,
    aucMaxRank = K_top_single,
    verbose = FALSE
  )

  cells_in_order <- colnames(seu_ct)
  score_default <- as.numeric(AUCell::getAUC(auc_res_default)[gene_symbol, cells_in_order])
  score_single <- as.numeric(AUCell::getAUC(auc_res_single)[gene_symbol, cells_in_order])
  names(score_default) <- cells_in_order
  names(score_single) <- cells_in_order

  explore_res <- tryCatch(
    AUCell::AUCell_exploreThresholds(
      auc_res_default,
      plotHist = FALSE,
      smallestPopPercent = 0.25,
      assignCells = TRUE,
      verbose = FALSE
    ),
    error = function(e) {
      message("  AUCell_exploreThresholds failed for ", ct_name, "/", gene_symbol,
              "; fallback to mean + 2*SD: ", conditionMessage(e))
      NULL
    }
  )

  thr_info <- if (!is.null(explore_res) && !is.null(explore_res[[gene_symbol]])) {
    aucThr <- explore_res[[gene_symbol]]$aucThr
    sel <- as.numeric(aucThr$selected)
    method_name <- "AUCell_default"
    candidates_local <- NULL
    if (!is.null(aucThr$thresholds)) {
      thr_df <- as.data.frame(aucThr$thresholds, stringsAsFactors = FALSE)
      thr_col <- if ("threshold" %in% colnames(thr_df)) "threshold" else colnames(thr_df)[1]
      thr_vals <- suppressWarnings(as.numeric(thr_df[[thr_col]]))
      candidates_local <- data.frame(
        method = rownames(thr_df),
        threshold = thr_vals,
        stringsAsFactors = FALSE
      )
      if (length(thr_vals) && any(is.finite(thr_vals))) {
        method_name <- rownames(thr_df)[which.min(abs(thr_vals - sel))]
      }
    }
    if (identical(method_name, "tenPercentOfMax") && !is.null(candidates_local)) {
      gk1_idx <- which(candidates_local$method == "Global_k1")
      if (length(gk1_idx)) {
        gk1_thr <- candidates_local$threshold[gk1_idx[1]]
        if (is.finite(gk1_thr)) {
          sel <- gk1_thr
          method_name <- "Global_k1"
        }
      }
    }
    list(selected = sel, method = method_name)
  } else {
    s <- score_default[is.finite(score_default)]
    sel <- if (length(s) > 1 && stats::sd(s) > 0) mean(s) + 2 * stats::sd(s) else max(s, 0) + .Machine$double.eps
    list(selected = sel, method = "mean_plus_2sd_fallback")
  }

  expr_count <- as.numeric(counts_ct[genes_found, cells_in_order, drop = TRUE])
  names(expr_count) <- cells_in_order

  count_positive <- expr_count > 0
  count_positive[!is.finite(expr_count)] <- FALSE

  aucell_threshold_positive <- score_default > thr_info$selected
  aucell_threshold_positive[!is.finite(score_default)] <- FALSE

  rank_positive <- score_single > 0
  rank_positive[!is.finite(score_single)] <- FALSE

  # This is the same protected single-gene positivity rule used for IL6 in
  # step4/step5. Threshold-only positivity is retained for audit columns.
  aucell_rule_positive <- rank_positive | count_positive

  data.frame(
    cell = cells_in_order,
    Group = factor(as.character(meta_ct[cells_in_order, "Group", drop = TRUE]), levels = display_levels),
    Patient = as.character(meta_ct[cells_in_order, "Patient", drop = TRUE]),
    Celltypes_raw = as.character(meta_ct[cells_in_order, "Celltypes_raw", drop = TRUE]),
    Celltype_use = ct_name,
    Celltype_label = pretty_ct(ct_name),
    Celltype_fig_label = pretty_ct_fig(ct_name),
    feature_set = gene_symbol,
    features_found = genes_found,
    CAMK2A_count = expr_count,
    CAMK2A_count_positive = as.logical(count_positive),
    CAMK2A_AUCell_score = score_single,
    CAMK2A_AUCell_score_default_window = score_default,
    CAMK2A_AUCell_threshold = as.numeric(thr_info$selected),
    CAMK2A_AUCell_threshold_method = thr_info$method,
    CAMK2A_AUCell_threshold_positive = as.logical(aucell_threshold_positive),
    CAMK2A_AUCell_rank_positive = as.logical(rank_positive),
    CAMK2A_AUCell_positive = as.logical(aucell_rule_positive),
    counts = as.logical(count_positive),
    rank_only_top0.05 = as.logical(score_default > 0),
    rank_only_top0.10 = as.logical(rank_positive),
    K_top_default = K_top_default,
    K_top_single = K_top_single,
    top_percent_default = top_percent_default,
    top_percent_single = top_percent_single,
    stringsAsFactors = FALSE
  )
}

message("Computing CAMK2A single-gene AUCell positivity for CTs: ", paste(cts_present, collapse = ", "))
camk2a_cell_level <- bind_rows(lapply(cts_present, function(ct) {
  compute_single_gene_score_df(
    obj,
    ct_name = ct,
    gene_symbol = target_gene,
    top_percent_default = top_percent_use,
    top_percent_single = top_percent_use_single
  )
}))

camk2a_cell_level <- camk2a_cell_level %>%
  mutate(
    Celltype_use = factor(Celltype_use, levels = ct_order),
    Celltype_fig_label = factor(Celltype_fig_label, levels = fig_levels),
    expected_celltype = expected_celltype,
    is_expected_celltype = as.character(Celltype_use) == expected_celltype
  )

if (isTRUE(export_cell_level_source)) {
  write.csv(
    camk2a_cell_level,
    file.path(source_tab_dir, sprintf("camk2a_cell_level_positivity_source_%s.csv", timestamp_tag)),
    row.names = FALSE
  )
}

# ---------------------------------------------------------
# 4. Composition source data and donut figures
# ---------------------------------------------------------
make_composition_source <- function(df, positive_col, definition_label) {
  total_positive <- sum(df[[positive_col]] %in% TRUE, na.rm = TRUE)
  out <- df %>%
    filter(.data[[positive_col]] %in% TRUE) %>%
    count(Celltype_use, Celltype_label, Celltype_fig_label, name = "positive_cells") %>%
    right_join(
      data.frame(
        Celltype_use = factor(ct_order, levels = ct_order),
        Celltype_label = unname(ct_full_map[ct_order]),
        Celltype_fig_label = factor(unname(ct_fig_map[ct_order]), levels = fig_levels),
        stringsAsFactors = FALSE
      ),
      by = c("Celltype_use", "Celltype_label", "Celltype_fig_label")
    ) %>%
    mutate(
      positive_cells = ifelse(is.na(positive_cells), 0L, as.integer(positive_cells)),
      total_positive_cells = total_positive,
      definition = definition_label,
      expected_celltype = expected_celltype,
      is_expected_celltype = as.character(Celltype_use) == expected_celltype
    ) %>%
    arrange(match(as.character(Celltype_use), ct_order))

  out$percent_of_positive_cells <- compute_percent_vector(out$positive_cells)

  out %>%
    select(
      definition,
      Celltype_use,
      Celltype_label,
      Celltype_fig_label,
      positive_cells,
      percent_of_positive_cells,
      total_positive_cells,
      expected_celltype,
      is_expected_celltype
    )
}

count_comp <- make_composition_source(camk2a_cell_level, "CAMK2A_count_positive", "count_positive")
rank05_comp <- make_composition_source(camk2a_cell_level, "rank_only_top0.05", "rank_only_top0.05")
rank10_comp <- make_composition_source(camk2a_cell_level, "rank_only_top0.10", "rank_only_top0.10")

count_comp$definition <- "counts"

write.csv(
  count_comp,
  file.path(source_tab_dir, sprintf("camk2a_counts_composition_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)
write.csv(
  rank05_comp,
  file.path(source_tab_dir, sprintf("camk2a_rank_only_top0.05_composition_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)
write.csv(
  rank10_comp,
  file.path(source_tab_dir, sprintf("camk2a_rank_only_top0.10_composition_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)

make_composition_donut <- function(comp_df, title_text, output_file) {
  plot_df <- comp_df %>%
    mutate(
      Celltype_fig_label = factor(as.character(Celltype_fig_label), levels = fig_levels),
      legend_label = paste0(
        as.character(Celltype_fig_label),
        " - ",
        format(positive_cells, big.mark = ","),
        " (",
        sprintf("%.1f%%", percent_of_positive_cells),
        ")"
      ),
      label = ifelse(positive_cells > 0, paste0(round(percent_of_positive_cells, 1), "%"), "")
    )

  legend_labels <- setNames(plot_df$legend_label, as.character(plot_df$Celltype_fig_label))

  p <- ggplot(plot_df, aes(x = 2, y = percent_of_positive_cells, fill = Celltype_fig_label)) +
    geom_col(color = "black", width = 1, linewidth = 0.3) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.7) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4.1, family = plot_font_family) +
    scale_fill_manual(
      values = cell_colors_fig,
      breaks = fig_levels,
      labels = legend_labels[fig_levels],
      drop = FALSE,
      name = NULL
    ) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 10.5, family = plot_font_family),
      legend.key.size = grid::unit(7.2, "mm"),
      legend.spacing.x = grid::unit(1, "mm"),
      legend.box.margin = margin(0, 0, 0, -6, unit = "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, size = 13.5, face = "bold", family = plot_font_family),
      plot.title.position = "plot",
      plot.margin = margin(8, 18, 5, 5, unit = "mm")
    ) +
    guides(fill = guide_legend(ncol = 1)) +
    ggtitle(title_text)

  save_plot_pdf(p, output_file, width = 178, height = 188, units = "mm")
  invisible(p)
}

fig_count <- make_composition_donut(
  count_comp,
  "Count-defined CAMK2A+ composition (%)",
  file.path(fig_dir, sprintf("camk2a_counts_composition_donut_%s.pdf", timestamp_tag))
)

fig_rank05 <- make_composition_donut(
  rank05_comp,
  "AUCell-ranked (top 5%) CAMK2A+ composition (%)",
  file.path(fig_dir, sprintf("camk2a_rank_only_top0.05_composition_donut_%s.pdf", timestamp_tag))
)

fig_rank10 <- make_composition_donut(
  rank10_comp,
  "AUCell-ranked (top 10%) CAMK2A+ composition (%)",
  file.path(fig_dir, sprintf("camk2a_rank_only_top0.10_composition_donut_%s.pdf", timestamp_tag))
)

composition_summary <- bind_rows(count_comp, rank05_comp, rank10_comp) %>%
  group_by(definition) %>%
  summarise(
    expected_celltype = expected_celltype[1],
    total_positive_cells = max(total_positive_cells, na.rm = TRUE),
    dominant_celltype = as.character(Celltype_use[which.max(positive_cells)]),
    dominant_celltype_label = as.character(Celltype_fig_label[which.max(positive_cells)]),
    dominant_celltype_positive_cells = max(positive_cells, na.rm = TRUE),
    dominant_celltype_percent = percent_of_positive_cells[which.max(positive_cells)],
    expected_celltype_positive_cells = positive_cells[as.character(Celltype_use) == expected_celltype][1],
    expected_celltype_percent = percent_of_positive_cells[as.character(Celltype_use) == expected_celltype][1],
    .groups = "drop"
  )

write.csv(
  composition_summary,
  file.path(summary_tab_dir, sprintf("camk2a_composition_summary_%s.csv", timestamp_tag)),
  row.names = FALSE
)

# ---------------------------------------------------------
# 5. Concordance table and rank-window summary figure
# ---------------------------------------------------------
calc_concordance <- function(df, scope_label, aucell_col, aucell_definition) {
  count_pos <- df$CAMK2A_count_positive %in% TRUE
  aucell_pos <- df[[aucell_col]] %in% TRUE

  both_positive <- sum(count_pos & aucell_pos, na.rm = TRUE)
  count_only <- sum(count_pos & !aucell_pos, na.rm = TRUE)
  aucell_only <- sum(!count_pos & aucell_pos, na.rm = TRUE)
  both_negative <- sum(!count_pos & !aucell_pos, na.rm = TRUE)
  count_positive <- sum(count_pos, na.rm = TRUE)
  aucell_positive <- sum(aucell_pos, na.rm = TRUE)
  union_positive <- sum(count_pos | aucell_pos, na.rm = TRUE)

  precision <- ifelse(aucell_positive > 0, both_positive / aucell_positive, NA_real_)
  recall <- ifelse(count_positive > 0, both_positive / count_positive, NA_real_)
  f1 <- ifelse(is.finite(precision + recall) && (precision + recall) > 0,
               2 * precision * recall / (precision + recall), NA_real_)
  jaccard <- ifelse(union_positive > 0, both_positive / union_positive, NA_real_)

  data.frame(
    scope = scope_label,
    aucell_definition = aucell_definition,
    n_cells = nrow(df),
    count_positive = count_positive,
    aucell_positive = aucell_positive,
    both_positive = both_positive,
    count_only = count_only,
    aucell_only = aucell_only,
    both_negative = both_negative,
    precision = precision,
    recall = recall,
    F1 = f1,
    Jaccard = jaccard,
    expected_celltype = expected_celltype,
    stringsAsFactors = FALSE
  )
}

df_all <- camk2a_cell_level
df_ext <- camk2a_cell_level %>% filter(as.character(Celltype_use) == expected_celltype)

concordance_source <- bind_rows(
  calc_concordance(df_all, "all_cells", "rank_only_top0.05", "rank_only_top0.05"),
  calc_concordance(df_ext, "ExcitNeu_only", "rank_only_top0.05", "rank_only_top0.05"),
  calc_concordance(df_all, "all_cells", "rank_only_top0.10", "rank_only_top0.10"),
  calc_concordance(df_ext, "ExcitNeu_only", "rank_only_top0.10", "rank_only_top0.10")
)

write.csv(
  concordance_source,
  file.path(source_tab_dir, sprintf("camk2a_concordance_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)

concordance_summary <- concordance_source %>%
  mutate(
    interpretation = dplyr::case_when(
      scope == "all_cells" ~ "Global concordance between count-defined positivity and the AUCell rank-only rule.",
      scope == "ExcitNeu_only" ~ "Expected-celltype concordance within excitatory neurons.",
      TRUE ~ "Concordance summary."
    )
  ) %>%
  select(scope, aucell_definition, precision, recall, F1, Jaccard, interpretation)

write.csv(
  concordance_summary,
  file.path(summary_tab_dir, sprintf("camk2a_concordance_summary_%s.csv", timestamp_tag)),
  row.names = FALSE
)

rank_window_metric_source <- concordance_source %>%
  filter(aucell_definition %in% c("rank_only_top0.05", "rank_only_top0.10")) %>%
  select(scope, aucell_definition, precision, recall, F1) %>%
  pivot_longer(cols = c(precision, recall, F1), names_to = "metric", values_to = "value") %>%
  mutate(
    scope = factor(scope, levels = c("all_cells", "ExcitNeu_only"), labels = c("All cells", "Excitatory neurons")),
    aucell_definition = factor(aucell_definition, levels = c("rank_only_top0.05", "rank_only_top0.10"),
                               labels = c("Top 5%", "Top 10%")),
    metric = factor(metric, levels = c("precision", "recall", "F1"),
                    labels = c("Precision", "Recall", "F1"))
  )

rank_window_metric_segments <- rank_window_metric_source %>%
  group_by(scope, metric) %>%
  summarise(x = min(value), xend = max(value), .groups = "drop")

count_baseline_comp <- count_comp %>%
  transmute(Celltype_use = as.character(Celltype_use), counts_percent = percent_of_positive_cells)

rank_window_shift_source <- bind_rows(rank05_comp, rank10_comp) %>%
  left_join(count_baseline_comp, by = "Celltype_use") %>%
  mutate(
    delta_vs_counts = percent_of_positive_cells - counts_percent,
    aucell_definition = factor(definition, levels = c("rank_only_top0.05", "rank_only_top0.10"),
                               labels = c("Top 5%", "Top 10%")),
    Celltype_use = factor(as.character(Celltype_use), levels = rev(c("ExcitNeu", "InhibNeu", "Ast", "OL", "Others")))
  ) %>%
  select(
    Celltype_use,
    Celltype_label,
    Celltype_fig_label,
    aucell_definition,
    positive_cells,
    percent_of_positive_cells,
    counts_percent,
    delta_vs_counts
  )

ct_topdown <- c("ExcitNeu", "InhibNeu", "Ast", "OL", "Others")
y_map <- setNames(c(5, 4, 3, 2, 1), ct_topdown)

rank_window_shift_segments <- rank_window_shift_source %>%
  mutate(Celltype_use_chr = as.character(Celltype_use)) %>%
  select(Celltype_use_chr, aucell_definition, delta_vs_counts) %>%
  tidyr::pivot_wider(names_from = aucell_definition, values_from = delta_vs_counts) %>%
  mutate(y_base = unname(y_map[Celltype_use_chr]))

rank_window_shift_points <- rank_window_shift_source %>%
  mutate(
    Celltype_use_chr = as.character(Celltype_use),
    x_plot = ifelse(as.character(aucell_definition) == "Top 5%", delta_vs_counts - 0.10, delta_vs_counts + 0.10),
    y_base = unname(y_map[Celltype_use_chr]),
    y_plot = ifelse(as.character(aucell_definition) == "Top 5%", y_base + 0.12, y_base - 0.12)
  )

write.csv(
  rank_window_metric_source,
  file.path(source_tab_dir, sprintf("camk2a_rank_window_metric_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)
write.csv(
  rank_window_shift_source,
  file.path(source_tab_dir, sprintf("camk2a_rank_window_composition_shift_source_%s.csv", timestamp_tag)),
  row.names = FALSE
)

rank_fill_cols <- c("Top 5%" = "#EBC39A", "Top 10%" = "#A9D9CF")

theme_pub <- theme_classic(base_size = 11) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 9, family = plot_font_family),
    plot.title = element_text(face = "bold", size = 12, family = plot_font_family, hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5, family = plot_font_family),
    axis.title = element_text(size = 11, family = plot_font_family),
    axis.text = element_text(size = 10, family = plot_font_family),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

rank_window_metric_plot <- ggplot(rank_window_metric_source, aes(x = value, y = metric)) +
  geom_segment(
    data = rank_window_metric_segments,
    aes(x = x, xend = xend, y = metric, yend = metric),
    inherit.aes = FALSE,
    color = "#D0D0D0",
    linewidth = 0.7
  ) +
  geom_point(
    aes(fill = aucell_definition),
    shape = 21,
    color = "black",
    stroke = 0.35,
    size = 3.2,
    position = position_dodge(width = 0.42)
  ) +
  facet_wrap(~scope, nrow = 1) +
  scale_fill_manual(values = rank_fill_cols) +
  scale_x_continuous(
    limits = c(0.85, 1.002),
    breaks = c(0.85, 0.90, 0.95, 1.00),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Concordance with count-defined positivity",
    x = "Metric value",
    y = NULL
  ) +
  theme_pub +
  theme(panel.spacing = grid::unit(12, "mm"))

rank_window_shift_plot <- ggplot() +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  geom_segment(
    data = rank_window_shift_segments,
    aes(x = `Top 5%`, xend = `Top 10%`, y = y_base, yend = y_base),
    color = "#D0D0D0",
    linewidth = 0.8
  ) +
  geom_point(
    data = rank_window_shift_points,
    aes(x = x_plot, y = y_plot, fill = aucell_definition),
    shape = 21,
    color = "black",
    stroke = 0.35,
    size = 3.25
  ) +
  scale_fill_manual(values = rank_fill_cols) +
  scale_y_continuous(
    breaks = c(5, 4, 3, 2, 1),
    labels = c("ExcitNeu", "InhibNeu", "Ast", "OL", "Others"),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(xlim = c(-8.0, 8.0), clip = "off") +
  scale_x_continuous(
    breaks = c(-8, -4, 0, 4, 8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Composition shift from count-defined baseline",
    x = "Percentage-point shift",
    y = NULL
  ) +
  theme_pub +
  theme(
    plot.title = element_text(hjust = 0.39)
  )

rank_window_summary_plot <- rank_window_metric_plot + plot_spacer() + rank_window_shift_plot +
  plot_layout(widths = c(1.18, 0.18, 1.02)) +
  plot_annotation(title = "CAMK2A comparison supports a more stable AUCell top-rank window at 10% than at 5%")

save_plot_pdf(
  rank_window_summary_plot,
  file.path(fig_dir, sprintf("camk2a_rank_window_stability_summary_%s.pdf", timestamp_tag)),
  width = 250,
  height = 110,
  units = "mm"
)

# ---------------------------------------------------------
# 6. Logs and RDS bundle
# ---------------------------------------------------------
readme_lines <- c(
  "Batiuk human step6 final outputs: CAMK2A AUCell threshold stability validation",
  paste0("Timestamp: ", timestamp_tag),
  "",
  "Input:",
  paste0("  Seurat RDS: ", normalizePath(human_rds, winslash = "/", mustWork = FALSE)),
  paste0("  Step4/step5 result root: ", normalizePath(result_root, winslash = "/", mustWork = FALSE)),
  paste0("  AUCell table directory detected: ", normalizePath(tab_src, winslash = "/", mustWork = FALSE)),
  "",
  "Definitions:",
  "  counts: raw CAMK2A count > 0.",
  "  rank_only_top0.05: CAMK2A AUC > 0 using the default top 5% rank window.",
  "  rank_only_top0.10: CAMK2A AUC > 0 using the extended top 10% rank window",
  "    used for single-gene positivity in step5.",
  "",
  "Parameters:",
  paste0("  top_percent_default = ", top_percent_use),
  paste0("  top_percent_single = ", top_percent_use_single),
  paste0("  export_cell_level_source = ", export_cell_level_source),
  paste0("  target_gene = ", target_gene),
  paste0("  expected_celltype = ", expected_celltype),
  "",
  "Directory structure:",
  "  figures/              Three composition donut PDFs and one rank-window summary PDF.",
  "  tables/source_data/   Cell-level positivity table, composition source data, concordance source data, and summary-figure source data.",
  "  tables/summary/       Short summary tables for composition and concordance.",
  "  logs/                 Run information and sessionInfo.",
  "  rds/                  Step6 object bundle."
)

writeLines(readme_lines, con = file.path(out_root, "README.txt"))

run_info <- c(
  readme_lines,
  "",
  "Cell type counts:",
  capture.output(print(table(obj$Celltype_use, useNA = "ifany"))),
  "",
  "Composition summary:",
  capture.output(print(composition_summary)),
  "",
  "Concordance summary:",
  capture.output(print(concordance_summary)),
  "",
  "Session info:",
  capture.output(sessionInfo())
)

writeLines(
  run_info,
  con = file.path(log_dir, sprintf("step6_final_batiuk_run_info_%s.txt", timestamp_tag))
)

step6_bundle <- list(
  timestamp_tag = timestamp_tag,
  human_rds = human_rds,
  result_root = result_root,
  out_root = out_root,
  target_gene = target_gene,
  expected_celltype = expected_celltype,
  top_percent_use = top_percent_use,
  top_percent_use_single = top_percent_use_single,
  camk2a_cell_level = camk2a_cell_level,
  count_composition = count_comp,
  rank05_composition = rank05_comp,
  rank10_composition = rank10_comp,
  composition_summary = composition_summary,
  concordance_source = concordance_source,
  concordance_summary = concordance_summary,
  rank_window_metric_source = rank_window_metric_source,
  rank_window_shift_source = rank_window_shift_source,
  figures = list(
    counts_composition = fig_count,
    rank05_composition = fig_rank05,
    rank10_composition = fig_rank10,
    rank_window_summary = rank_window_summary_plot
  )
)

saveRDS(
  step6_bundle,
  file.path(rds_dir, sprintf("step6_final_batiuk_objects_%s.rds", timestamp_tag))
)

message("========== STEP6 FINAL BATIUK COMPLETE ==========")
message("Output root: ", normalizePath(out_root, winslash = "/", mustWork = FALSE))
