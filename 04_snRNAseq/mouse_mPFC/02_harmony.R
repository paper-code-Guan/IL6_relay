# ==============================================================================
# Script: Mouse scRNA-seq integration, clustering, and marker discovery
# Project: Mouse mPFC snRNA-seq analysis
# Purpose:
#   This script integrates QC-filtered mouse scRNA-seq samples using Harmony,
#   performs UMAP visualization and Louvain clustering, and identifies cluster
#   marker genes for downstream cell-type annotation and figure generation.
#
# Notes:
#   - Core input paths, output paths, dimensions, resolutions, and marker
#     discovery parameters are kept consistent with the original script.
#   - Outputs are written with explicit file paths to improve reproducibility.
#   - Comments and messages are kept in English for SCI code archiving.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# ==============================================================================
# 1. Configuration
# ==============================================================================

data_directory <- Sys.getenv("MOUSE_OUTPUT_DIR", unset = "")
if (!nzchar(data_directory)) {
  stop("MOUSE_OUTPUT_DIR must be set.", call. = FALSE)
}
GROUP_name     <- "ALL"
species_flag   <- "mouse"   # Reserved for project-level bookkeeping.

n_features             <- 2000
variation_method       <- "vst"
dimension_values       <- c(35)
clustering_resolutions <- c(0.8)
max_dimensions         <- max(dimension_values)

harmony_batch_variable <- "orig.ident"

# Marker discovery parameters. These are intentionally kept unchanged.
marker_test_use        <- "wilcox"
marker_min_pct         <- 0.1
marker_logfc_threshold <- 0.25
marker_only_pos        <- TRUE
marker_min_cells_group <- 3
marker_adj_p_cutoff    <- 0.05
top_marker_n           <- 50

input_rds <- file.path(data_directory, "all_samples_merged.rds")

results_directory <- file.path(
  data_directory,
  paste0("2-harmony_", GROUP_name),
  paste0("res_", paste(clustering_resolutions, collapse = "_"))
)

dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Utility functions
# ==============================================================================

message_header <- function(text) {
  cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
  cat(text, "\n", sep = "")
  cat(paste(rep("=", 72), collapse = ""), "\n", sep = "")
}

sanitize_resolution <- function(x) {
  gsub("\\.", "p", as.character(x))
}

save_session_info <- function(out_file) {
  writeLines(capture.output(sessionInfo()), con = out_file)
}

run_harmony_integration <- function(seurat_object, dimensions_to_use, batch_variable) {
  if (!batch_variable %in% colnames(seurat_object@meta.data)) {
    stop("Batch variable not found in metadata: ", batch_variable, call. = FALSE)
  }

  cat("Running Harmony integration using batch variable: ", batch_variable, "\n", sep = "")

  if ("RunHarmony" %in% getNamespaceExports("harmony")) {
    return(
      harmony::RunHarmony(
        object         = seurat_object,
        reduction      = "pca",
        group.by.vars  = batch_variable,
        reduction.save = "harmony",
        dims.use       = dimensions_to_use,
        verbose        = FALSE
      )
    )
  }

  stop(
    "RunHarmony was not found in the installed harmony package.",
    call. = FALSE
  )
}

identify_cluster_markers <- function(seurat_object, cluster_id) {
  tryCatch({
    cluster_markers <- FindMarkers(
      object          = seurat_object,
      ident.1         = cluster_id,
      test.use        = marker_test_use,
      min.pct         = marker_min_pct,
      logfc.threshold = marker_logfc_threshold,
      only.pos        = marker_only_pos,
      min.cells.group = marker_min_cells_group
    )

    if (nrow(cluster_markers) == 0) {
      return(data.frame())
    }

    cluster_markers$gene    <- rownames(cluster_markers)
    cluster_markers$cluster <- as.character(cluster_id)

    cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), , drop = FALSE]
  }, error = function(e) {
    message(sprintf("FindMarkers failed for cluster %s: %s", cluster_id, e$message))
    data.frame()
  })
}

write_marker_outputs <- function(markers_collection, current_dimensions, current_resolution, out_dir) {
  resolution_label <- as.character(current_resolution)

  all_markers_combined <- dplyr::bind_rows(markers_collection)

  rds_filename <- file.path(
    out_dir,
    sprintf("dim_%d_res_%s.harmony.DE_marker.rds", current_dimensions, resolution_label)
  )
  saveRDS(all_markers_combined, rds_filename)

  significant_markers <- lapply(markers_collection, function(marker_df) {
    if (nrow(marker_df) > 0 && "p_val_adj" %in% colnames(marker_df)) {
      marker_df[marker_df$p_val_adj < marker_adj_p_cutoff, , drop = FALSE]
    } else {
      marker_df
    }
  })

  excel_filename <- file.path(
    out_dir,
    sprintf("dim_%d_res_%s.harmony.markers.xlsx", current_dimensions, resolution_label)
  )
  openxlsx::write.xlsx(significant_markers, file = excel_filename, asTable = TRUE)

  significant_counts <- sapply(significant_markers, nrow)
  summary_dataframe <- data.frame(
    cluster       = names(significant_counts),
    n_significant = as.integer(significant_counts),
    stringsAsFactors = FALSE
  )

  summary_filename <- file.path(
    out_dir,
    sprintf("dim_%d_res_%s.harmony.markers_summary.csv", current_dimensions, resolution_label)
  )
  utils::write.csv(summary_dataframe, file = summary_filename, row.names = FALSE)

  cat("Significant markers per cluster:\n")
  print(significant_counts)

  list(
    all_markers_combined = all_markers_combined,
    significant_markers  = significant_markers,
    summary_dataframe    = summary_dataframe,
    rds_file             = rds_filename,
    excel_file           = excel_filename,
    summary_file         = summary_filename
  )
}

write_top_marker_outputs <- function(all_markers_combined, current_dimensions, current_resolution, out_dir) {
  resolution_label <- as.character(current_resolution)

  required_columns <- c("cluster", "avg_log2FC", "p_val_adj", "gene")
  missing_columns <- setdiff(required_columns, colnames(all_markers_combined))

  if (length(missing_columns) > 0 || nrow(all_markers_combined) == 0) {
    warning(
      "Top marker extraction skipped because marker table is empty or missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
    return(invisible(NULL))
  }

  filtered_markers <- all_markers_combined %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(p_val_adj < marker_adj_p_cutoff) %>%
    filter(!is.na(avg_log2FC))

  top_markers <- filtered_markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = top_marker_n) %>%
    ungroup()

  excel_sheets <- split(top_markers, paste0("Cluster_", top_markers$cluster))

  all_cluster_ids <- sort(unique(as.character(all_markers_combined$cluster)))
  missing_clusters <- setdiff(all_cluster_ids, unique(top_markers$cluster))

  if (length(missing_clusters) > 0) {
    for (cluster_id in missing_clusters) {
      excel_sheets[[paste0("Cluster_", cluster_id)]] <- data.frame(
        info = sprintf(
          "No genes passed filters or fewer than %d markers were available.",
          top_marker_n
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(excel_sheets) > 0) {
    sheet_order <- order(suppressWarnings(as.numeric(gsub("Cluster_", "", names(excel_sheets)))))
    excel_sheets <- excel_sheets[sheet_order]
  }

  top_markers_excel <- file.path(
    out_dir,
    sprintf("dim_%d_res_%s.harmony.top50.xlsx", current_dimensions, resolution_label)
  )
  openxlsx::write.xlsx(excel_sheets, top_markers_excel, asTable = TRUE)

  top_markers_csv <- file.path(
    out_dir,
    sprintf("dim_%d_res_%s.harmony.top50.csv", current_dimensions, resolution_label)
  )
  utils::write.csv(top_markers, top_markers_csv, row.names = FALSE)

  cat("Top markers saved to:\n")
  cat("  ", top_markers_excel, "\n", sep = "")
  cat("  ", top_markers_csv, "\n", sep = "")

  invisible(list(excel = top_markers_excel, csv = top_markers_csv))
}

# ==============================================================================
# 3. Load QC-filtered merged object
# ==============================================================================

message_header("Loading QC-filtered merged single-cell object")

if (!file.exists(input_rds)) {
  stop("Input RDS file not found: ", input_rds, call. = FALSE)
}

combined_object <- readRDS(input_rds)
if (!inherits(combined_object, "Seurat")) {
  stop("Loaded object is not a Seurat object.", call. = FALSE)
}

DefaultAssay(combined_object) <- "RNA"

cat("Input file: ", input_rds, "\n", sep = "")
cat("Initial object dimensions: ", paste(dim(combined_object), collapse = " genes x cells = "), "\n", sep = "")
cat("Output directory: ", results_directory, "\n", sep = "")

# Save a copy of analysis parameters for reproducibility.
analysis_parameters <- data.frame(
  parameter = c(
    "GROUP_name",
    "species_flag",
    "n_features",
    "variation_method",
    "dimension_values",
    "clustering_resolutions",
    "max_dimensions",
    "harmony_batch_variable",
    "marker_test_use",
    "marker_min_pct",
    "marker_logfc_threshold",
    "marker_only_pos",
    "marker_min_cells_group",
    "marker_adj_p_cutoff",
    "top_marker_n"
  ),
  value = c(
    GROUP_name,
    species_flag,
    n_features,
    variation_method,
    paste(dimension_values, collapse = ","),
    paste(clustering_resolutions, collapse = ","),
    max_dimensions,
    harmony_batch_variable,
    marker_test_use,
    marker_min_pct,
    marker_logfc_threshold,
    marker_only_pos,
    marker_min_cells_group,
    marker_adj_p_cutoff,
    top_marker_n
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  analysis_parameters,
  file.path(results_directory, "analysis_parameters.csv"),
  row.names = FALSE
)

# ==============================================================================
# 4. Standard preprocessing and Harmony integration
# ==============================================================================

message_header("Preprocessing: normalization, HVG selection, scaling, and PCA")

combined_object <- NormalizeData(combined_object, verbose = FALSE)
combined_object <- FindVariableFeatures(
  combined_object,
  selection.method = variation_method,
  nfeatures = n_features,
  verbose = FALSE
)
combined_object <- ScaleData(
  combined_object,
  features = VariableFeatures(combined_object),
  verbose = FALSE
)
combined_object <- RunPCA(
  combined_object,
  reduction.name = "pca",
  npcs = max_dimensions,
  verbose = FALSE
)

message_header("Harmony integration")
combined_object <- run_harmony_integration(
  seurat_object     = combined_object,
  dimensions_to_use = 1:max_dimensions,
  batch_variable    = harmony_batch_variable
)

# ==============================================================================
# 5. Dimension reduction, clustering, visualization, and marker discovery
# ==============================================================================

all_run_summaries <- list()

for (current_dimensions in dimension_values) {
  message_header(sprintf("Processing Harmony dimensions: %d", current_dimensions))
  set.seed(1000 + current_dimensions)

  umap_reduction_name <- paste0("umap_", current_dimensions)
  graph_name          <- paste0("RNA_snn_", current_dimensions)

  combined_object <- RunUMAP(
    combined_object,
    reduction = "harmony",
    dims = 1:current_dimensions,
    reduction.name = umap_reduction_name,
    seed.use = 123,
    verbose = FALSE
  )

  combined_object <- FindNeighbors(
    combined_object,
    reduction = "harmony",
    dims = 1:current_dimensions,
    graph.name = graph_name,
    verbose = FALSE
  )

  for (current_resolution in clustering_resolutions) {
    message_header(sprintf(
      "Clustering and marker discovery: dims=%d, resolution=%s",
      current_dimensions,
      as.character(current_resolution)
    ))

    set.seed(2000 + current_dimensions * 100 + as.integer(current_resolution * 100))

    combined_object <- FindClusters(
      combined_object,
      graph.name = graph_name,
      resolution = current_resolution,
      algorithm = 1,
      verbose = FALSE
    )

    cluster_column_name <- paste0(
      "seurat_clusters_dim", current_dimensions,
      "_res", sanitize_resolution(current_resolution)
    )
    combined_object@meta.data[[cluster_column_name]] <- combined_object$seurat_clusters
    Idents(combined_object) <- cluster_column_name

    cluster_identifiers <- levels(Idents(combined_object))
    cluster_cell_counts <- table(Idents(combined_object))

    cat("Identified clusters: ", length(cluster_identifiers), "\n", sep = "")
    cat("Cluster cell counts:\n")
    print(sort(cluster_cell_counts, decreasing = TRUE))

    cluster_counts_df <- data.frame(
      cluster = names(cluster_cell_counts),
      n_cells = as.integer(cluster_cell_counts),
      stringsAsFactors = FALSE
    )
    utils::write.csv(
      cluster_counts_df,
      file.path(
        results_directory,
        sprintf("dim_%d_res_%s.harmony.cluster_cell_counts.csv", current_dimensions, as.character(current_resolution))
      ),
      row.names = FALSE
    )

    # UMAP quality-control plots.
    sample_umap_plot <- DimPlot(
      combined_object,
      group.by = "orig.ident",
      reduction = umap_reduction_name,
      raster = FALSE
    ) +
      ggtitle("By sample") +
      theme_classic()

    cluster_umap_plot <- DimPlot(
      combined_object,
      group.by = cluster_column_name,
      label = TRUE,
      label.size = 3,
      reduction = umap_reduction_name,
      raster = FALSE
    ) +
      ggtitle("By cluster") +
      theme_classic()

    combined_plot <- gridExtra::arrangeGrob(sample_umap_plot, cluster_umap_plot, ncol = 2)
    plot_filename <- file.path(
      results_directory,
      sprintf("dim_%d_res_%s.harmony.umap.pdf", current_dimensions, as.character(current_resolution))
    )
    ggsave(
      filename = plot_filename,
      plot = combined_plot,
      width = 14,
      height = 7
    )

    # Save the processed Seurat object for this parameter set.
    object_filename <- file.path(
      results_directory,
      sprintf("dim_%d_res_%s.harmony.rds", current_dimensions, as.character(current_resolution))
    )
    saveRDS(combined_object, file = object_filename)

    cat("Saved processed object: ", object_filename, "\n", sep = "")

    # Marker discovery.
    cat("Identifying positive cluster marker genes using FindMarkers...\n")

    valid_clusters <- names(cluster_cell_counts[cluster_cell_counts >= marker_min_cells_group])
    small_clusters <- names(cluster_cell_counts[cluster_cell_counts < marker_min_cells_group])

    if (length(small_clusters) > 0) {
      warning(
        "Skipping marker discovery for clusters with < ", marker_min_cells_group,
        " cells: ", paste(small_clusters, collapse = ", ")
      )
    }

    markers_collection <- setNames(
      vector("list", length(cluster_identifiers)),
      paste0("Cluster_", cluster_identifiers)
    )

    for (cluster_id in cluster_identifiers) {
      sheet_name <- paste0("Cluster_", cluster_id)
      if (cluster_id %in% valid_clusters) {
        markers_collection[[sheet_name]] <- identify_cluster_markers(combined_object, cluster_id)
      } else {
        markers_collection[[sheet_name]] <- data.frame()
      }
    }

    marker_outputs <- write_marker_outputs(
      markers_collection  = markers_collection,
      current_dimensions  = current_dimensions,
      current_resolution  = current_resolution,
      out_dir             = results_directory
    )

    write_top_marker_outputs(
      all_markers_combined = marker_outputs$all_markers_combined,
      current_dimensions   = current_dimensions,
      current_resolution   = current_resolution,
      out_dir              = results_directory
    )

    all_run_summaries[[paste0("dim", current_dimensions, "_res", current_resolution)]] <- data.frame(
      dimensions     = current_dimensions,
      resolution     = current_resolution,
      n_clusters     = length(cluster_identifiers),
      n_cells        = ncol(combined_object),
      n_genes        = nrow(combined_object),
      object_file    = object_filename,
      umap_file      = plot_filename,
      marker_rds     = marker_outputs$rds_file,
      marker_excel   = marker_outputs$excel_file,
      marker_summary = marker_outputs$summary_file,
      stringsAsFactors = FALSE
    )
  }
}

# ==============================================================================
# 6. Final summary and reproducibility records
# ==============================================================================

message_header("Final summary")

run_summary_df <- dplyr::bind_rows(all_run_summaries)
utils::write.csv(
  run_summary_df,
  file.path(results_directory, "integration_run_summary.csv"),
  row.names = FALSE
)

cat("All output files were saved in:\n", results_directory, "\n", sep = "")
cat("Run summary:\n")
print(run_summary_df)

last_dimensions <- tail(dimension_values, 1)
last_resolution <- tail(clustering_resolutions, 1)
last_cluster_column <- paste0(
  "seurat_clusters_dim", last_dimensions,
  "_res", sanitize_resolution(last_resolution)
)

if (last_cluster_column %in% colnames(combined_object@meta.data)) {
  cat("Final cluster distribution for ", last_cluster_column, ":\n", sep = "")
  print(table(combined_object@meta.data[[last_cluster_column]]))
}

save_session_info(file.path(results_directory, "sessionInfo_step2_harmony.txt"))

cat("\nIntegration and marker discovery completed successfully.\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
