# ==============================================================================
# Mouse single-nucleus RNA-seq annotation and marker-guided refinement pipeline
#
# Input : Harmony-integrated Seurat object and cluster marker file from Step 2
# Output: Annotated/refined Seurat object for downstream analyses
#         dim_35_res_0.8_harmony_CelltypesV3.rds
#
# Notes:
#   - Core paths, dimensions, resolution, and manual annotation decisions are kept
#     consistent with the original analysis script.
#   - Final cell-type labels are kept as short labels used by downstream scripts:
#     Ast, ExcitNeu, InhibNeu, OL, Mic, OPC, Endo, Lept, Pericyte.
#   - Three marker-driven refinement steps are applied sequentially:
#       1) Endo_mix  -> Endo / Pericyte / Mic / ExcitNeu
#       2) OL        -> OL / ExcitNeu
#       3) Mic_mix   -> Mic / ExcitNeu
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(stringr)
  library(openxlsx)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

# Optional packages used when available
has_scMayoMap <- requireNamespace("scMayoMap", quietly = TRUE)
has_harmony   <- requireNamespace("harmony", quietly = TRUE)
has_cairo     <- capabilities("cairo")

options(stringsAsFactors = FALSE)
set.seed(123)

# ==============================================================================
# 1. Configuration
# ==============================================================================

data_directory <- Sys.getenv("MOUSE_OUTPUT_DIR", unset = "")
if (!nzchar(data_directory)) {
  stop("MOUSE_OUTPUT_DIR must be set.", call. = FALSE)
}
GROUP_name            <- "ALL"
selected_resolution   <- 0.8
selected_dimensions   <- 35
species               <- "mouse"

results_directory <- file.path(
  data_directory,
  paste0("2-harmony_", GROUP_name),
  sprintf("res_%.1f", selected_resolution)
)

stopifnot("Results directory does not exist." = dir.exists(results_directory))

fig_dir   <- file.path(results_directory, "annotation_figures")
table_dir <- file.path(results_directory, "annotation_tables")
log_dir   <- file.path(results_directory, "annotation_logs")
sub_dir   <- file.path(results_directory, "subclustering_refinement")

dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(sub_dir,   recursive = TRUE, showWarnings = FALSE)

object_file <- file.path(
  results_directory,
  sprintf("dim_%d_res_%.1f.harmony.rds", selected_dimensions, selected_resolution)
)

marker_file <- file.path(
  results_directory,
  sprintf("dim_%d_res_%.1f.harmony.DE_marker.rds", selected_dimensions, selected_resolution)
)

stopifnot("Harmony Seurat object not found." = file.exists(object_file))
stopifnot("Cluster marker RDS not found."    = file.exists(marker_file))

# Key output filenames preserved for downstream scripts
annotated_object_file <- file.path(
  results_directory,
  sprintf("dim_%d_res_%.1f_harmony_Celltypes.annotated.rds", selected_dimensions, selected_resolution)
)

final_object_file <- file.path(
  results_directory,
  sprintf("dim_%d_res_%.1f_harmony_CelltypesV3.rds", selected_dimensions, selected_resolution)
)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

safe_ggsave <- function(filename, plot, width, height, units = "in") {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  device_use <- if (has_cairo) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    units    = units,
    device   = device_use,
    limitsize = FALSE
  )
}

get_umap_reduction <- function(seurat_obj, dims) {
  candidates <- c(paste0("umap_", dims), "umap")
  hit <- candidates[candidates %in% Reductions(seurat_obj)]
  if (length(hit) == 0) {
    stop("No UMAP reduction found. Available reductions: ", paste(Reductions(seurat_obj), collapse = ", "))
  }
  hit[1]
}

get_cluster_column <- function(seurat_obj, dims, res) {
  res_tag <- gsub("\\.", "p", as.character(res))
  candidates <- c(
    paste0("seurat_clusters_dim", dims, "_res", res_tag),
    "seurat_clusters"
  )
  hit <- candidates[candidates %in% colnames(seurat_obj@meta.data)]
  if (length(hit) == 0) {
    stop("No cluster column found. Checked: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

match_genes_case_insensitive <- function(genes, reference_genes) {
  idx <- match(toupper(genes), toupper(reference_genes))
  unique(reference_genes[stats::na.omit(idx)])
}

match_genes_robust <- function(genes, reference_genes) {
  reference_upper <- toupper(reference_genes)
  out <- character(0)
  missing <- character(0)

  for (gene in genes) {
    gene_upper <- toupper(gene)
    exact_idx <- which(reference_upper == gene_upper)
    if (length(exact_idx)) {
      out <- c(out, reference_genes[exact_idx[1]])
      next
    }

    boundary_idx <- grep(paste0("^", gene_upper, "($|[-_.])"), reference_upper)
    if (length(boundary_idx)) {
      out <- c(out, reference_genes[boundary_idx[1]])
    } else {
      missing <- c(missing, gene)
    }
  }

  if (length(missing)) {
    message("Unmatched genes: ", paste(missing, collapse = ", "))
  }

  unique(out)
}

pretty_celltype <- function(x) {
  dplyr::recode(
    as.character(x),
    "Ast"       = "Astrocytes",
    "ExcitNeu"  = "Excitatory neurons",
    "InhibNeu"  = "Inhibitory neurons",
    "OL"        = "Oligodendrocytes",
    "Mic"       = "Microglia",
    "OPC"       = "OPCs",
    "Endo"      = "Endothelial cells",
    "Lept"      = "Leptomeningeal cells",
    "Pericyte"  = "Pericytes",
    .default = as.character(x)
  )
}

save_celltype_table <- function(seurat_obj, label, column = "Celltypes") {
  tab <- as.data.frame(table(seurat_obj@meta.data[[column]]), stringsAsFactors = FALSE)
  colnames(tab) <- c(column, "n_cells")
  tab$celltype_label <- pretty_celltype(tab[[column]])
  utils::write.csv(
    tab,
    file.path(table_dir, paste0(label, "_celltype_counts.csv")),
    row.names = FALSE
  )
  invisible(tab)
}

plot_labeled_umap <- function(seurat_obj, group_col, filename, title = NULL) {
  reduction_use <- get_umap_reduction(seurat_obj, selected_dimensions)
  emb <- Embeddings(seurat_obj, reduction = reduction_use) %>% as.data.frame()
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
    colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  }

  df <- emb %>%
    mutate(cell_type = as.character(seurat_obj@meta.data[[group_col]]))

  label_df <- df %>%
    group_by(cell_type) %>%
    summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = pretty_celltype(cell_type))

  pal <- scales::hue_pal()(length(unique(df$cell_type)))
  names(pal) <- sort(unique(df$cell_type))

  p <- ggplot(df, aes(UMAP_1, UMAP_2, colour = cell_type)) +
    geom_point(size = 0.3, alpha = 0.85, shape = 16) +
    scale_colour_manual(values = pal) +
    ggrepel::geom_label_repel(
      data = label_df,
      aes(UMAP_1, UMAP_2, label = label),
      inherit.aes = FALSE,
      size = 3.5,
      label.size = 0.25,
      segment.color = NA,
      colour = "black"
    ) +
    labs(title = title, x = "UMAP1", y = "UMAP2") +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 16),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  safe_ggsave(filename, p, width = 8, height = 8)
  invisible(p)
}

plot_marker_dotplot <- function(seurat_obj, marker_panels, filename, width = 20, height = 12) {
  DefaultAssay(seurat_obj) <- "RNA"

  matched_panels <- lapply(marker_panels, match_genes_case_insensitive, reference_genes = rownames(seurat_obj))
  matched_panels <- matched_panels[sapply(matched_panels, length) > 0]

  if (length(matched_panels) == 0) {
    warning("No marker genes matched the Seurat object; skipping dot plot.")
    return(invisible(NULL))
  }

  p <- DotPlot(seurat_obj, features = matched_panels, scale = TRUE) +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
    theme(
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10, face = "italic"),
      legend.title = element_text(face = "bold", size = 11)
    ) +
    coord_flip() +
    labs(x = NULL, y = NULL)

  safe_ggsave(filename, p, width = width, height = height)
  invisible(p)
}

# Module-score based fallback annotation for unmapped subclusters.
auto_annotate_by_modules <- function(seurat_obj) {
  marker_sets <- list(
    ExcitNeu = c("Slc17a7", "Camk2a", "Satb2", "Rorb"),
    InhibNeu = c("Gad1", "Gad2", "Slc32a1", "Slc6a1"),
    Ast = c("Aqp4", "Gfap", "Slc1a2", "Glul"),
    Endo = c("Cldn5", "Flt1", "Pecam1", "Slc2a1"),
    Mic = c("C1qa", "Cx3cr1", "Tmem119", "Csf1r"),
    OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
    OL = c("Mog", "Mbp", "Plp1", "Cldn11"),
    Pericyte = c("Kcnj8", "Abcc9", "Pdgfrb", "Notch3"),
    Lept = c("Bnc2", "Col12a1", "Col15a1", "Lum")
  )

  matched_sets <- lapply(marker_sets, match_genes_case_insensitive, reference_genes = rownames(seurat_obj))
  matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  if (length(matched_sets) == 0) {
    stop("No module-score marker set matched the subset object.")
  }

  seurat_obj <- AddModuleScore(seurat_obj, features = matched_sets, name = "module_score_")
  score_cols <- grep("^module_score_\\d+$", colnames(seurat_obj@meta.data), value = TRUE)
  score_cols <- score_cols[order(as.integer(sub("^module_score_", "", score_cols)))]

  if (length(score_cols) != length(matched_sets)) {
    stop("Unexpected module score column count.")
  }

  score_df <- seurat_obj@meta.data[, score_cols, drop = FALSE]
  score_df$cluster <- as.character(seurat_obj$seurat_clusters)

  agg <- score_df %>%
    group_by(cluster) %>%
    summarise(across(all_of(score_cols), median, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()

  rownames(agg) <- agg$cluster
  agg$cluster <- NULL
  colnames(agg) <- names(matched_sets)

  best <- colnames(agg)[max.col(as.matrix(agg), ties.method = "first")]
  setNames(best, rownames(agg))
}

run_subcluster_refinement <- function(main_obj,
                                      target_label,
                                      selection_col,
                                      manual_map,
                                      output_tag,
                                      dims_sub = 30,
                                      resolution_sub = 0.5,
                                      k_param = 10,
                                      use_harmony = TRUE) {
  message("\n---- Refining target: ", target_label, " from column: ", selection_col, " ----")

  if (!selection_col %in% colnames(main_obj@meta.data)) {
    stop("Selection column not found: ", selection_col)
  }

  target_cells <- rownames(main_obj@meta.data)[as.character(main_obj@meta.data[[selection_col]]) == target_label]
  if (length(target_cells) == 0) {
    warning("Target label not found; skipping: ", target_label)
    return(main_obj)
  }

  message("Subset cell count: ", length(target_cells))
  subset_obj <- subset(main_obj, cells = target_cells)
  DefaultAssay(subset_obj) <- "RNA"

  subset_obj <- NormalizeData(subset_obj, verbose = FALSE)
  subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  subset_obj <- ScaleData(subset_obj, features = VariableFeatures(subset_obj), verbose = FALSE)

  max_pcs <- min(50L, ncol(subset_obj) - 1L, length(VariableFeatures(subset_obj)) - 1L)
  max_pcs <- max(10L, max_pcs)
  subset_obj <- RunPCA(subset_obj, npcs = max_pcs, verbose = FALSE)

  reduction_used <- "pca"
  if (use_harmony && has_harmony && "orig.ident" %in% colnames(subset_obj@meta.data)) {
    subset_obj <- harmony::RunHarmony(
      subset_obj,
      group.by.vars = "orig.ident",
      reduction = "pca",
      dims.use = 1:min(dims_sub, max_pcs),
      reduction.save = "harmony_sub",
      project.dim = FALSE,
      verbose = FALSE
    )
    reduction_used <- "harmony_sub"
  }

  subset_obj <- RunUMAP(
    subset_obj,
    reduction = reduction_used,
    dims = 1:min(dims_sub, max_pcs),
    reduction.name = "umap_sub",
    seed.use = 123,
    verbose = FALSE
  )

  adaptive_k <- max(10L, min(k_param, floor(ncol(subset_obj) * 0.1), ncol(subset_obj) - 1L))
  subset_obj <- FindNeighbors(
    subset_obj,
    reduction = reduction_used,
    dims = 1:min(dims_sub, max_pcs),
    k.param = adaptive_k,
    verbose = FALSE
  )
  subset_obj <- FindClusters(subset_obj, resolution = resolution_sub, algorithm = 1, verbose = FALSE)

  cluster_ids <- as.character(levels(subset_obj$seurat_clusters))
  message("Subcluster distribution:")
  print(table(subset_obj$seurat_clusters))

  auto_map <- auto_annotate_by_modules(subset_obj)
  final_map <- auto_map
  final_map[names(manual_map)] <- manual_map

  missing_manual <- setdiff(cluster_ids, names(manual_map))
  if (length(missing_manual)) {
    warning(
      "Some subclusters are not manually mapped for ", target_label,
      "; module-score fallback will be used for: ", paste(missing_manual, collapse = ", ")
    )
  }

  subset_obj$automated_celltype <- auto_map[as.character(subset_obj$seurat_clusters)]
  subset_obj$Celltypes_refined  <- final_map[as.character(subset_obj$seurat_clusters)]
  Idents(subset_obj) <- subset_obj$Celltypes_refined

  # Save subset-level outputs
  safe_ggsave(
    file.path(fig_dir, paste0("Figure_sub_", output_tag, "_umap.pdf")),
    DimPlot(subset_obj, reduction = "umap_sub", label = TRUE, label.size = 4, raster = FALSE) +
      ggtitle(paste0("Sub-UMAP: ", target_label)),
    width = 10,
    height = 8
  )

  mapping_df <- data.frame(
    target_label = target_label,
    subcluster = names(final_map),
    automated_annotation = unname(auto_map[names(final_map)]),
    refined_annotation = unname(final_map),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    mapping_df,
    file.path(table_dir, paste0("Sub_", output_tag, "_automated_vs_refined_mapping.csv")),
    row.names = FALSE
  )

  size_df <- as.data.frame(table(subset_obj$seurat_clusters), stringsAsFactors = FALSE)
  colnames(size_df) <- c("subcluster", "n_cells")
  utils::write.csv(
    size_df,
    file.path(table_dir, paste0("Sub_", output_tag, "_cluster_sizes.csv")),
    row.names = FALSE
  )

  saveRDS(
    subset_obj,
    file.path(sub_dir, paste0("subclustered_", output_tag, "_analysis.rds"))
  )

  # Integrate refined labels back into the main object.
  if (!"Celltypes_refined" %in% colnames(main_obj@meta.data)) {
    main_obj$Celltypes_refined <- as.character(main_obj$Celltypes)
  } else {
    main_obj$Celltypes_refined <- as.character(main_obj$Celltypes_refined)
  }

  main_obj$Celltypes_refined[Cells(subset_obj)] <- as.character(subset_obj$Celltypes_refined)
  main_obj$Celltypes_refined <- factor(main_obj$Celltypes_refined)
  main_obj$Celltypes <- factor(as.character(main_obj$Celltypes_refined))
  main_obj$Celltype  <- main_obj$Celltypes
  Idents(main_obj) <- main_obj$Celltypes

  save_celltype_table(main_obj, paste0("after_", output_tag), column = "Celltypes")

  main_obj
}

# ==============================================================================
# 3. Load Step 2 outputs
# ==============================================================================

message("Loading Harmony object: ", object_file)
seurat_object <- readRDS(object_file)
DefaultAssay(seurat_object) <- "RNA"

cluster_column <- get_cluster_column(seurat_object, selected_dimensions, selected_resolution)
umap_reduction <- get_umap_reduction(seurat_object, selected_dimensions)
message("Cluster column: ", cluster_column)
message("UMAP reduction: ", umap_reduction)

Idents(seurat_object) <- seurat_object@meta.data[[cluster_column]]
seurat_object$seurat_clusters <- Idents(seurat_object)

# Create group metadata used in downstream plotting.
seurat_object$sample <- seurat_object$orig.ident
seurat_object$group <- ifelse(grepl("^con", seurat_object$orig.ident, ignore.case = TRUE), "Control", "Treatment")

# ==============================================================================
# 4. Optional scMayoMap annotation
# ==============================================================================

markers_df <- readRDS(marker_file)
stopifnot("Required marker columns missing." = all(c("gene", "cluster") %in% colnames(markers_df)))

if (has_scMayoMap) {
  message("Running scMayoMap annotation...")
  data("scMayoMapDatabase", package = "scMayoMap")
  scMayo_res <- scMayoMap::scMayoMap(
    data = markers_df,
    database = scMayoMapDatabase,
    tissue = "brain"
  )

  scMayo_df <- as.data.frame(scMayo_res$res)
  utils::write.csv(scMayo_df, file.path(table_dir, "scMayoMap_results.csv"), row.names = FALSE)

  if (all(c("cluster", "celltype", "ES.norm") %in% colnames(scMayo_df))) {
    top_scMayo <- scMayo_df %>%
      group_by(cluster) %>%
      slice_max(order_by = ES.norm, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(
        cluster_chr = as.character(cluster),
        celltype_main = trimws(vapply(strsplit(as.character(celltype), ";"), function(x) tail(x, 1), character(1)))
      )

    seurat_object$scMayo_celltype_main <- top_scMayo$celltype_main[
      match(as.character(seurat_object$seurat_clusters), top_scMayo$cluster_chr)
    ]
  }

  scMayo_plot <- scMayoMap::scMayoMap.plot(scMayoMap.object = scMayo_res)
  safe_ggsave(file.path(fig_dir, "Fig_scMayoMap_cellType.pdf"), scMayo_plot, width = 20, height = 9)
} else {
  warning("Package scMayoMap is not installed; scMayoMap annotation will be skipped.")
  seurat_object$scMayo_celltype_main <- NA_character_
}

# ==============================================================================
# 5. Initial manual annotation of major clusters
# ==============================================================================

initial_cluster_map <- c(
  "0"  = "ExcitNeu",
  "1"  = "ExcitNeu",
  "2"  = "ExcitNeu",
  "3"  = "ExcitNeu",
  "4"  = "ExcitNeu",
  "5"  = "ExcitNeu",
  "6"  = "ExcitNeu",
  "7"  = "ExcitNeu",
  "8"  = "ExcitNeu",
  "9"  = "InhibNeu",
  "10" = "OL",
  "11" = "ExcitNeu",
  "12" = "ExcitNeu",
  "13" = "Ast",
  "14" = "InhibNeu",
  "15" = "ExcitNeu",
  "16" = "InhibNeu",
  "17" = "Mic_mix",
  "18" = "Endo_mix",
  "19" = "InhibNeu",
  "20" = "OPC",
  "21" = "InhibNeu",
  "22" = "InhibNeu",
  "23" = "Lept",
  "24" = "OL",
  "25" = "InhibNeu",
  "26" = "Ast"
)

cluster_ids <- as.character(seurat_object$seurat_clusters)
missing_clusters <- setdiff(unique(cluster_ids), names(initial_cluster_map))
if (length(missing_clusters)) {
  stop("The following clusters are not covered by initial_cluster_map: ", paste(missing_clusters, collapse = ", "))
}

seurat_object$original_celltypes <- initial_cluster_map[cluster_ids]
seurat_object$Celltypes_refined  <- seurat_object$original_celltypes
seurat_object$Celltypes          <- factor(seurat_object$Celltypes_refined)
seurat_object$Celltype           <- seurat_object$Celltypes
Idents(seurat_object) <- seurat_object$Celltypes

initial_map_df <- data.frame(
  cluster = names(initial_cluster_map),
  initial_celltype = unname(initial_cluster_map),
  stringsAsFactors = FALSE
)
utils::write.csv(initial_map_df, file.path(table_dir, "initial_manual_cluster_annotation_map.csv"), row.names = FALSE)
save_celltype_table(seurat_object, "initial_manual_annotation", column = "Celltypes")

safe_ggsave(
  file.path(fig_dir, "Fig_clusters_vs_initial_annotation.pdf"),
  DimPlot(seurat_object, group.by = cluster_column, reduction = umap_reduction, label = TRUE, label.size = 3, raster = FALSE) +
    ggtitle("Clusters") +
    DimPlot(seurat_object, group.by = "Celltypes", reduction = umap_reduction, label = TRUE, label.size = 3, raster = FALSE) +
    ggtitle("Initial cell-type annotation"),
  width = 16,
  height = 8
)

saveRDS(seurat_object, file = annotated_object_file)
message("Initial annotated object saved to: ", annotated_object_file)

# ==============================================================================
# 6. Marker-driven subcluster refinement
# ==============================================================================

endo_manual_map <- c(
  "0" = "ExcitNeu",
  "1" = "ExcitNeu",
  "2" = "ExcitNeu",
  "3" = "ExcitNeu",
  "4" = "Endo",
  "5" = "Pericyte",
  "6" = "Mic",
  "7" = "ExcitNeu"
)

ol_manual_map <- c(
  "0"  = "OL",
  "1"  = "OL",
  "2"  = "OL",
  "3"  = "OL",
  "4"  = "OL",
  "5"  = "OL",
  "6"  = "OL",
  "7"  = "ExcitNeu",
  "8"  = "OL",
  "9"  = "OL",
  "10" = "OL"
)

mic_manual_map <- c(
  "0"  = "Mic",
  "1"  = "Mic",
  "2"  = "Mic",
  "3"  = "Mic",
  "4"  = "ExcitNeu",
  "5"  = "Mic",
  "6"  = "Mic",
  "7"  = "Mic",
  "8"  = "Mic",
  "9"  = "Mic",
  "10" = "Mic"
)

seurat_object <- run_subcluster_refinement(
  main_obj = seurat_object,
  target_label = "Endo_mix",
  selection_col = "Celltypes_refined",
  manual_map = endo_manual_map,
  output_tag = "Endo_mix",
  dims_sub = 30,
  resolution_sub = 0.5,
  k_param = 10,
  use_harmony = TRUE
)

seurat_object <- run_subcluster_refinement(
  main_obj = seurat_object,
  target_label = "OL",
  selection_col = "Celltypes_refined",
  manual_map = ol_manual_map,
  output_tag = "OL",
  dims_sub = 30,
  resolution_sub = 0.5,
  k_param = 10,
  use_harmony = TRUE
)

seurat_object <- run_subcluster_refinement(
  main_obj = seurat_object,
  target_label = "Mic_mix",
  selection_col = "Celltypes_refined",
  manual_map = mic_manual_map,
  output_tag = "Mic_mix",
  dims_sub = 30,
  resolution_sub = 0.5,
  k_param = 10,
  use_harmony = TRUE
)

# ==============================================================================
# 7. Final object standardization and validation
# ==============================================================================

final_levels <- c(
  "ExcitNeu", "InhibNeu", "OL", "Ast", "Mic",
  "OPC", "Endo", "Lept", "Pericyte"
)

final_labels_raw <- as.character(seurat_object$Celltypes_refined)

if (any(is.na(final_labels_raw))) {
  stop("NA values detected in Celltypes_refined before final factor conversion.")
}

unexpected_labels <- setdiff(unique(final_labels_raw), final_levels)

if (length(unexpected_labels) > 0) {
  stop(
    "Unexpected final cell-type labels detected before factor conversion: ",
    paste(unexpected_labels, collapse = ", ")
  )
}

seurat_object$Celltypes_refined <- factor(final_labels_raw, levels = final_levels)
seurat_object$Celltypes <- seurat_object$Celltypes_refined
seurat_object$Celltype  <- seurat_object$Celltypes
Idents(seurat_object) <- seurat_object$Celltypes

if (any(is.na(seurat_object$Celltypes))) {
  stop("NA values detected in final Celltypes after factor conversion.")
}

final_counts <- save_celltype_table(seurat_object, "final_CelltypesV3", column = "Celltypes")
message("Final cell-type counts:")
print(final_counts)

# Export final object metadata table for traceability.
metadata_export <- seurat_object@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(
    cell,
    any_of(c("orig.ident", "group", "sample", cluster_column, "scMayo_celltype_main",
             "original_celltypes", "Celltypes_refined", "Celltypes", "Celltype"))
  )
utils::write.csv(metadata_export, file.path(table_dir, "final_CelltypesV3_metadata.csv"), row.names = FALSE)

# Save final object using the historical filename required by downstream scripts.
saveRDS(seurat_object, file = final_object_file)
message("Final refined object saved to: ", final_object_file)

# ==============================================================================
# 8. Final marker and UMAP visualizations
# ==============================================================================

astrocyte_markers       <- c("Aqp4", "Slc1a2", "Slc1a3", "Gfap", "Acsbg1", "Glul", "Atp1b2", "Ndrg2", "Htra1")
microglia_markers       <- c("C1qa", "Cx3cr1", "Csf1r", "Trem2", "Itgam", "C1qb", "P2ry12", "Aif1", "C1qc", "Tmem119", "Ptprc")
opc_markers             <- c("Pdgfra", "Olig1", "Olig2", "C1ql1", "Cspg4")
oligodendrocyte_markers <- c("Mog", "Mbp", "Mag", "Cldn11", "Mobp", "Opalin", "Aspa", "Fa2h", "Cnp")
endothelial_markers     <- c("Flt1", "Pecam1", "Vwf", "Slc2a1", "Cldn5", "Slco1a4", "Ramp2", "Epas1")
pericyte_markers        <- c("Kcnj8", "Abcc9", "Pdgfrb", "Notch3", "Higd1b", "Vtn")
leptomeningeal_markers  <- c("Bnc2", "Col12a1", "Col15a1", "Cemip", "Lum", "Cped1", "Tbx18", "Igfbp6", "Col6a2", "Dcn")
excitatory_markers      <- c("Slc17a7", "Camk2a", "Satb2", "Cux2", "Bcl11b")
inhibitory_markers      <- c("Gad2", "Gad1", "Slc32a1", "Pvalb", "Sst", "Vip", "Reln", "Lhx6")

marker_panels <- list(
  Astrocytes = astrocyte_markers,
  Excitatory_neurons = excitatory_markers,
  Inhibitory_neurons = inhibitory_markers,
  Oligodendrocytes = oligodendrocyte_markers,
  Microglia = microglia_markers,
  OPCs = opc_markers,
  Endothelial_cells = endothelial_markers,
  Leptomeningeal_cells = leptomeningeal_markers,
  Pericytes = pericyte_markers
)

plot_marker_dotplot(
  seurat_object,
  marker_panels,
  file.path(fig_dir, "Fig_DotPlot_all_markers_final_CelltypesV3.pdf"),
  width = 20,
  height = 14
)

plot_labeled_umap(
  seurat_object,
  group_col = "Celltypes",
  filename = file.path(fig_dir, "Fig_UMAP_final_CelltypesV3.pdf"),
  title = "Mouse cell-type annotation"
)

# Cell composition by experimental group.
prop_df <- as.data.frame(prop.table(table(seurat_object$Celltypes, seurat_object$group), margin = 2))
colnames(prop_df) <- c("Celltypes", "group", "proportion")
prop_df$Celltype_label <- pretty_celltype(prop_df$Celltypes)
utils::write.csv(prop_df, file.path(table_dir, "final_celltype_proportions_by_group.csv"), row.names = FALSE)

p_prop <- ggplot(prop_df, aes(x = group, y = proportion, fill = Celltype_label)) +
  geom_col(width = 0.8) +
  labs(x = NULL, y = "Fraction of cells", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(size = 9)
  )

safe_ggsave(file.path(fig_dir, "Fig_CellFraction_by_group_final_CelltypesV3.pdf"), p_prop, width = 10, height = 6)

# ==============================================================================
# 9. Session information
# ==============================================================================

writeLines(capture.output(sessionInfo()), con = file.path(log_dir, "sessionInfo_step3_annotation.txt"))

message("\n=== Step 3 annotation/refinement completed successfully ===")
message("Input object : ", object_file)
message("Final object: ", final_object_file)
message("Figures     : ", fig_dir)
message("Tables      : ", table_dir)

# ==============================================================================
# End of script
# ==============================================================================
