# ==============================================================================
# Human dlPFC snRNA-seq: Cell-type annotation after Harmony integration (SCI V4)
#
# Purpose:
#   1) Load the Harmony-integrated object and cluster marker table from Step 2.
#   2) Run scMayoMap-assisted annotation when available.
#   3) Apply the finalized global manual cluster-to-cell-type annotation first,
#      including cluster 18 = Other_mix.
#   4) Save the pre-refinement annotated2 object for traceability.
#   5) Immediately extract cells labeled Other_mix and perform the in-line local
#      subcluster refinement: subcluster 8 = Ast, all other subclusters = Other.
#   6) Export annotation tables, QC figures, and the final annotated3 Seurat object.
#
# Output object for downstream human AUCell and figure scripts:
#   dim_35_res_0.8_harmony_Celltypes.annotated3.rds
#
# Environment:
#   R 4.x, Seurat 4.2.1, harmony, dplyr, ggplot2, readr, stringr, tibble, scales
#   Optional: scMayoMap, tidyr
# ===============================================================================


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# ------------------------------------------------------------------------------
# 0) Configuration: keep paths and analysis parameters consistent with Step 2
# ------------------------------------------------------------------------------
reproduce_dir <- Sys.getenv("BATIUK_OUTPUT_DIR", unset = "")
if (!nzchar(reproduce_dir)) {
  stop("BATIUK_OUTPUT_DIR must be set.", call. = FALSE)
}
GROUP_name    <- "ALL"
res_select    <- 0.8
res_label     <- sprintf("%.1f", res_select)
dim_select    <- 35

res_dir <- file.path(
  reproduce_dir,
  paste0("2-harmony_", GROUP_name),
  paste0("res_", res_label)
)
stopifnot("Harmony result directory not found." = dir.exists(res_dir))

obj_file <- file.path(
  res_dir,
  sprintf("dim_%d_res_%s.harmony.rds", dim_select, res_label)
)
de_file <- file.path(
  res_dir,
  sprintf("dim_%d_res_%s.harmony.DE_marker.rds", dim_select, res_label)
)
stopifnot("Harmony Seurat object not found." = file.exists(obj_file))
stopifnot("Cluster marker RDS not found." = file.exists(de_file))

sample_manifest_csv <- Sys.getenv("BATIUK_SAMPLE_MANIFEST", unset = "")
if (!nzchar(sample_manifest_csv)) {
  stop("BATIUK_SAMPLE_MANIFEST must be set.", call. = FALSE)
}

annotation_table_dir  <- file.path(res_dir, "annotation_tables")
annotation_figure_dir <- file.path(res_dir, "annotation_figures")
annotation_log_dir    <- file.path(res_dir, "annotation_logs")
subcluster_dir        <- file.path(annotation_table_dir, "Other_mix_refinement")
for (dd in c(annotation_table_dir, annotation_figure_dir, annotation_log_dir, subcluster_dir)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

# Output names preserved for downstream compatibility.
annotated2_file <- file.path(
  res_dir,
  sprintf("dim_%d_res_%s_harmony_Celltypes.annotated2.rds", dim_select, res_label)
)
annotated3_file <- file.path(
  res_dir,
  sprintf("dim_%d_res_%s_harmony_Celltypes.annotated3.rds", dim_select, res_label)
)

# Allowed final cell-type categories defined by the annotation strategy.
# Note: final cell counts are intentionally NOT hard-coded in SCI V4, because they
# are outputs of the current run rather than prior assumptions. The script checks
# label validity and total cell conservation, then exports the observed counts.
final_levels <- c("Ast", "ExcitNeu", "InhibNeu", "OL", "Other")
intermediate_levels <- c(final_levels, "Other_mix")

# ------------------------------------------------------------------------------
# 1) Helper functions
# ------------------------------------------------------------------------------
safe_ggsave <- function(filename, plot, width, height, units = "in") {
  tryCatch(
    {
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = units,
        device = grDevices::cairo_pdf,
        limitsize = FALSE
      )
    },
    error = function(e) {
      message("cairo_pdf failed; falling back to pdf device: ", e$message)
      ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = width,
        height = height,
        units = units,
        device = "pdf",
        limitsize = FALSE
      )
    }
  )
}

get_umap_name <- function(object, dim_select) {
  candidates <- c(paste0("umap_", dim_select), "umap")
  hit <- candidates[candidates %in% Seurat::Reductions(object)]
  if (length(hit) == 0L) {
    stop("No UMAP reduction found. Available reductions: ",
         paste(Seurat::Reductions(object), collapse = ", "))
  }
  hit[1]
}

match_genes_case_insensitive <- function(genes, available_genes) {
  idx <- match(toupper(genes), toupper(available_genes))
  unique(stats::na.omit(available_genes[idx]))
}

resolve_cluster_column <- function(object, dim_select, res_label) {
  snapshot_col <- paste0(
    "seurat_clusters_dim",
    dim_select,
    "_res",
    gsub("\\.", "p", as.character(res_label))
  )
  if (snapshot_col %in% colnames(object@meta.data)) return(snapshot_col)
  if ("seurat_clusters" %in% colnames(object@meta.data)) return("seurat_clusters")
  stop("No clustering column found: neither ", snapshot_col, " nor seurat_clusters exists.")
}

add_sample_metadata_if_needed <- function(object, sample_manifest_csv) {
  if (!file.exists(sample_manifest_csv)) {
    warning("The Batiuk sample manifest was not found; returning the object without additional sample metadata: ",
            sample_manifest_csv)
    return(object)
  }
  if (!"sample" %in% colnames(object@meta.data)) {
    if ("orig.ident" %in% colnames(object@meta.data)) {
      object$sample <- object$orig.ident
    } else {
      warning("No sample/orig.ident column was found; the Batiuk sample manifest cannot be merged.")
      return(object)
    }
  }

  sample_info <- readr::read_csv(sample_manifest_csv, show_col_types = FALSE)
  if ("sample_id" %in% colnames(sample_info) &&
      !"Identifier" %in% colnames(sample_info)) {
    sample_info <- sample_info %>%
      dplyr::rename(Identifier = sample_id)
  }
  if ("diagnosis" %in% colnames(sample_info) &&
      !"Diagnosis" %in% colnames(sample_info)) {
    sample_info <- sample_info %>%
      dplyr::rename(Diagnosis = diagnosis)
  }
  if (!"Identifier" %in% colnames(sample_info)) {
    warning("The Batiuk sample manifest lacks a 'sample_id' or 'Identifier' column; sample metadata was not added.")
    return(object)
  }

  sample_vec <- as.character(object$sample)
  missing_ids <- setdiff(unique(sample_vec), sample_info$Identifier)
  if (length(missing_ids) > 0L) {
    warning("The following samples are absent from the Batiuk sample manifest: ",
            paste(missing_ids, collapse = ", "))
  }

  sample_info_indexed <- sample_info %>% tibble::column_to_rownames("Identifier")
  new_meta <- sample_info_indexed[sample_vec, , drop = FALSE]
  rownames(new_meta) <- colnames(object)

  # Avoid silently overwriting existing metadata columns from Step 1.
  for (nm in colnames(new_meta)) {
    if (!nm %in% colnames(object@meta.data) || all(is.na(object@meta.data[[nm]]))) {
      object@meta.data[[nm]] <- new_meta[[nm]]
    }
  }
  object
}

standardize_diagnosis_group <- function(object) {
  diagnosis_col <- NULL
  if ("Diagnosis" %in% colnames(object@meta.data)) {
    diagnosis_col <- "Diagnosis"
  } else if ("samp.Diagnosis" %in% colnames(object@meta.data)) {
    diagnosis_col <- "samp.Diagnosis"
  }

  if (is.null(diagnosis_col)) {
    warning("No Diagnosis/samp.Diagnosis column found; Group was not created.")
    return(object)
  }

  diagnosis_raw <- trimws(as.character(object@meta.data[[diagnosis_col]]))
  diagnosis_lower <- tolower(diagnosis_raw)

  group_std <- rep(NA_character_, length(diagnosis_raw))
  group_std[grepl("ctr|ctrl|control|healthy|hc", diagnosis_lower)] <- "CTR"
  group_std[grepl("scz|schizo|case", diagnosis_lower)] <- "SCZ"

  if (any(is.na(group_std))) {
    stop("Unrecognized diagnosis labels in column '", diagnosis_col, "': ",
         paste(unique(diagnosis_raw[is.na(group_std)]), collapse = ", "))
  }

  object$Group <- factor(group_std, levels = c("CTR", "SCZ"))
  object$GROUP <- object$Group
  object$group <- object$Group
  object
}

plot_basic_umap <- function(object, group_col, reduction, filename, width = 8, height = 7) {
  if (!group_col %in% colnames(object@meta.data)) return(invisible(NULL))
  p <- Seurat::DimPlot(
    object,
    group.by = group_col,
    reduction = reduction,
    label = TRUE,
    label.size = 3,
    raster = FALSE
  ) + ggplot2::ggtitle(group_col)
  safe_ggsave(filename, p, width = width, height = height)
  invisible(p)
}

validate_final_labels <- function(labels, allowed_levels, label_name = "Celltypes") {
  labels_chr <- as.character(labels)
  if (any(is.na(labels_chr))) {
    stop("NA values detected in ", label_name, " before final factor conversion.")
  }
  unexpected <- setdiff(unique(labels_chr), allowed_levels)
  if (length(unexpected) > 0L) {
    stop("Unexpected labels detected in ", label_name, ": ",
         paste(unexpected, collapse = ", "))
  }
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# 2) Load Harmony object and marker table
# ------------------------------------------------------------------------------
message("Loading Harmony object: ", obj_file)
scobj <- readRDS(obj_file)
stopifnot("Loaded object is not a Seurat object." = inherits(scobj, "Seurat"))
DefaultAssay(scobj) <- "RNA"

umap_name <- get_umap_name(scobj, dim_select)
cluster_col <- resolve_cluster_column(scobj, dim_select, res_label)
message("Using UMAP reduction: ", umap_name)
message("Using cluster column: ", cluster_col)

scobj$seurat_clusters <- factor(as.character(scobj@meta.data[[cluster_col]]))
Idents(scobj) <- scobj$seurat_clusters

message("Loading cluster markers: ", de_file)
de_markers <- readRDS(de_file)
stopifnot("Marker table must contain gene and cluster columns." =
            all(c("gene", "cluster") %in% colnames(de_markers)))
de_markers$cluster <- as.character(de_markers$cluster)

# Quick overview plots for annotation QC.
p_by_sample <- Seurat::DimPlot(
  scobj,
  group.by = if ("sample" %in% colnames(scobj@meta.data)) "sample" else "orig.ident",
  reduction = umap_name,
  raster = FALSE
) + ggplot2::ggtitle("By sample")
safe_ggsave(
  file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_by_sample_overview.pdf", dim_select, res_label)),
  p_by_sample,
  width = 8,
  height = 7
)

p_by_cluster <- Seurat::DimPlot(
  scobj,
  group.by = "seurat_clusters",
  reduction = umap_name,
  label = TRUE,
  label.size = 3,
  raster = FALSE
) + ggplot2::ggtitle("By cluster")
safe_ggsave(
  file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_by_cluster_overview.pdf", dim_select, res_label)),
  p_by_cluster,
  width = 8,
  height = 7
)

# ------------------------------------------------------------------------------
# 3) Optional scMayoMap-assisted annotation
# ------------------------------------------------------------------------------
scobj$scMayo_celltype_main <- NA_character_

if (requireNamespace("scMayoMap", quietly = TRUE)) {
  message("Running scMayoMap annotation...")
  tryCatch(
    {
      data("scMayoMapDatabase", package = "scMayoMap")
      scmm <- scMayoMap::scMayoMap(
        data = de_markers,
        database = scMayoMapDatabase,
        tissue = "brain"
      )

      mayo_plot <- scMayoMap::scMayoMap.plot(scMayoMap.object = scmm)
      safe_ggsave(
        file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_scMayoMap_plot.pdf", dim_select, res_label)),
        mayo_plot,
        width = 20,
        height = 9
      )

      mayo_res <- as.data.frame(scmm$res)
      if (all(c("cluster", "celltype", "ES.norm") %in% colnames(mayo_res)) &&
          requireNamespace("tidyr", quietly = TRUE)) {
        mayo_top <- mayo_res %>%
          mutate(across(everything(), as.character)) %>%
          tidyr::separate_rows(celltype, ES.norm, ES.raw, ES.se, sep = ";") %>%
          mutate(
            across(c(celltype, ES.norm, ES.raw, ES.se), ~ stringr::str_trim(.x)),
            across(c(ES.norm, ES.raw, ES.se), ~ suppressWarnings(as.numeric(.x))),
            cluster = as.character(cluster)
          ) %>%
          group_by(cluster) %>%
          slice_max(order_by = ES.norm, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          transmute(cluster, scMayo_celltype_main = celltype, ES.norm, ES.raw, ES.se)

        scobj$scMayo_celltype_main <- mayo_top$scMayo_celltype_main[
          match(as.character(scobj$seurat_clusters), mayo_top$cluster)
        ]

        readr::write_csv(
          mayo_top,
          file.path(annotation_table_dir, sprintf("dim_%d_res_%s_scMayoMap_top_annotation.csv", dim_select, res_label))
        )
      } else {
        message("scMayoMap results were produced, but top-label parsing was skipped because required columns or tidyr were unavailable.")
      }
    },
    error = function(e) {
      warning("scMayoMap annotation failed but the pipeline will continue with manual annotation: ", e$message)
    }
  )
} else {
  warning("Package 'scMayoMap' is not installed. Manual annotation will still be applied.")
}

# ------------------------------------------------------------------------------
# 4) Manual cluster-to-cell-type annotation before Other_mix refinement
# ------------------------------------------------------------------------------
# Finalized human annotation for ALL group, dim35, resolution 0.8.
# Cluster 18 is intentionally labeled as Other_mix here and refined in Section 6.
manual_annotation <- c(
  "0"  = "ExcitNeu",
  "1"  = "InhibNeu",
  "2"  = "ExcitNeu",
  "3"  = "ExcitNeu",
  "4"  = "InhibNeu",
  "5"  = "ExcitNeu",
  "6"  = "ExcitNeu",
  "7"  = "ExcitNeu",
  "8"  = "InhibNeu",
  "9"  = "InhibNeu",
  "10" = "InhibNeu",
  "11" = "InhibNeu",
  "12" = "ExcitNeu",
  "13" = "InhibNeu",
  "14" = "OL",
  "15" = "ExcitNeu",
  "16" = "ExcitNeu",
  "17" = "InhibNeu",
  "18" = "Other_mix",
  "19" = "InhibNeu",
  "20" = "ExcitNeu",
  "21" = "InhibNeu",
  "22" = "ExcitNeu",
  "23" = "Ast",
  "24" = "InhibNeu",
  "25" = "ExcitNeu",
  "26" = "ExcitNeu",
  "27" = "InhibNeu",
  "28" = "OL"
)

cluster_levels <- levels(scobj$seurat_clusters)
unknown_in_map <- setdiff(names(manual_annotation), cluster_levels)
missing_in_map <- setdiff(cluster_levels, names(manual_annotation))
if (length(unknown_in_map) > 0L) {
  stop("Manual map contains cluster IDs not present in object: ",
       paste(unknown_in_map, collapse = ", "))
}
if (length(missing_in_map) > 0L) {
  stop("Manual map does not cover all clusters: ",
       paste(missing_in_map, collapse = ", "))
}

manual_lookup <- manual_annotation[as.character(scobj$seurat_clusters)]
validate_final_labels(manual_lookup, intermediate_levels, label_name = "manual Celltypes")

scobj$Celltypes <- factor(manual_lookup, levels = intermediate_levels)
scobj$Celltypes_refined <- as.character(scobj$Celltypes)
scobj$Celltype <- scobj$Celltypes
scobj$original_celltypes <- as.character(scobj$Celltypes)
Idents(scobj) <- scobj$Celltypes

manual_map_df <- data.frame(
  cluster = names(manual_annotation),
  manual_celltype = unname(manual_annotation),
  stringsAsFactors = FALSE
)
readr::write_csv(
  manual_map_df,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_manual_cluster_annotation_map_before_refinement.csv", dim_select, res_label))
)

cluster_annotation_table_before <- data.frame(
  cell = colnames(scobj),
  cluster = as.character(scobj$seurat_clusters),
  scMayo_celltype_main = as.character(scobj$scMayo_celltype_main),
  Celltypes_before_refinement = as.character(scobj$Celltypes),
  stringsAsFactors = FALSE
)
readr::write_csv(
  cluster_annotation_table_before,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_celllevel_annotation_before_refinement.csv", dim_select, res_label))
)

manual_vs_scmayo <- as.data.frame(table(
  Cluster = as.character(scobj$seurat_clusters),
  scMayo = as.character(scobj$scMayo_celltype_main),
  ManualBeforeRefinement = as.character(scobj$Celltypes),
  useNA = "ifany"
))
readr::write_csv(
  manual_vs_scmayo,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_scMayo_vs_manual_before_refinement.csv", dim_select, res_label))
)

# ------------------------------------------------------------------------------
# 5) Add or verify sample metadata
# ------------------------------------------------------------------------------
scobj <- add_sample_metadata_if_needed(scobj, sample_manifest_csv)
scobj <- standardize_diagnosis_group(scobj)

# Save the pre-refinement annotated object for traceability after metadata checks.
saveRDS(scobj, file = annotated2_file)
message("Pre-refinement annotated object saved to: ", annotated2_file)

# ------------------------------------------------------------------------------
# 6) Restore Other_mix subcluster refinement from the original workflow
# ------------------------------------------------------------------------------
# This section is deliberately written in-line, directly after the global manual
# annotation and annotated2 saving step, to make the execution order unambiguous:
#   global clusters 0-28 are annotated first; cluster 18 is labeled as Other_mix;
#   then only cells labeled Other_mix are extracted for local subclustering.
target_label      <- "Other_mix"
meta_col_for_pick <- "Celltypes"
dims_sub          <- 30L
res_sub           <- 0.5
k_param_sub       <- 10L
use_harmony_sub   <- TRUE

if (!meta_col_for_pick %in% colnames(scobj@meta.data)) {
  stop("Column '", meta_col_for_pick, "' not found in object metadata.")
}

if (!target_label %in% unique(as.character(scobj@meta.data[[meta_col_for_pick]]))) {
  stop("Target label '", target_label, "' not found in '", meta_col_for_pick, "'.")
}

cells_use <- rownames(scobj@meta.data)[
  as.character(scobj@meta.data[[meta_col_for_pick]]) == target_label
]
n_cells_subset <- length(cells_use)

if (n_cells_subset < 20L) {
  stop("Too few cells for reliable subclustering: ", target_label,
       " (n = ", n_cells_subset, ")")
}

message("Subclustering ", target_label, ": ", n_cells_subset, " cells.")
sub_other <- subset(scobj, cells = cells_use)
DefaultAssay(sub_other) <- "RNA"

sub_other <- NormalizeData(sub_other, verbose = FALSE)
sub_other <- FindVariableFeatures(sub_other, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
sub_other <- ScaleData(sub_other, features = VariableFeatures(sub_other), verbose = FALSE)

n_var_features <- length(VariableFeatures(sub_other))
max_pcs <- min(50L, n_cells_subset - 1L, n_var_features - 1L)
n_pcs_use <- min(as.integer(dims_sub), max_pcs)

if (n_pcs_use < 5L) {
  stop("Too few usable PCs for subclustering: ", target_label,
       " (usable PCs = ", n_pcs_use, ")")
}

sub_other <- RunPCA(sub_other, npcs = max_pcs, verbose = FALSE)

batch_var <- NULL
if ("sample" %in% colnames(sub_other@meta.data) && length(unique(sub_other$sample)) > 1L) {
  batch_var <- "sample"
} else if ("orig.ident" %in% colnames(sub_other@meta.data) && length(unique(sub_other$orig.ident)) > 1L) {
  batch_var <- "orig.ident"
}

if (use_harmony_sub && !is.null(batch_var) && requireNamespace("harmony", quietly = TRUE)) {
  sub_other <- harmony::RunHarmony(
    sub_other,
    group.by.vars = batch_var,
    reduction = "pca",
    dims.use = 1:n_pcs_use,
    reduction.save = "harmony_sub",
    project.dim = FALSE,
    verbose = FALSE
  )
  reduction_use <- "harmony_sub"
} else {
  if (use_harmony_sub && is.null(batch_var)) {
    message("Harmony skipped for Other_mix subclustering: no multi-level sample/orig.ident variable.")
  }
  if (use_harmony_sub && !requireNamespace("harmony", quietly = TRUE)) {
    message("Harmony skipped for Other_mix subclustering: package 'harmony' is not installed.")
  }
  reduction_use <- "pca"
}

adaptive_k <- min(
  as.integer(k_param_sub),
  max(2L, floor(n_cells_subset * 0.1)),
  n_cells_subset - 1L
)
if (adaptive_k < 2L) stop("k.param is too small for subclustering: ", target_label)

sub_other <- RunUMAP(
  sub_other,
  reduction = reduction_use,
  dims = 1:n_pcs_use,
  reduction.name = "umap_sub",
  seed.use = 123,
  verbose = FALSE
)
sub_other <- FindNeighbors(
  sub_other,
  reduction = reduction_use,
  dims = 1:n_pcs_use,
  k.param = adaptive_k,
  verbose = FALSE
)
sub_other <- FindClusters(
  sub_other,
  resolution = res_sub,
  algorithm = 1,
  verbose = FALSE
)

sub_other$seurat_clusters <- factor(as.character(sub_other$seurat_clusters))
sub_levels <- levels(sub_other$seurat_clusters)

# Finalized local manual refinement from the original human annotation workflow.
manual_refinement <- c(
  "0"  = "Other",
  "1"  = "Other",
  "2"  = "Other",
  "3"  = "Other",
  "4"  = "Other",
  "5"  = "Other",
  "6"  = "Other",
  "7"  = "Other",
  "8"  = "Ast",
  "9"  = "Other",
  "10" = "Other",
  "11" = "Other"
)

extra_in_map <- setdiff(names(manual_refinement), sub_levels)
missing_in_map <- setdiff(sub_levels, names(manual_refinement))
if (length(extra_in_map) > 0L) {
  stop("manual_refinement contains subclusters not present in ", target_label, ": ",
       paste(extra_in_map, collapse = ", "))
}
if (length(missing_in_map) > 0L) {
  stop("manual_refinement does not cover all ", target_label, " subclusters: ",
       paste(missing_in_map, collapse = ", "))
}

sub_other$Celltypes_refined <- unname(manual_refinement[as.character(sub_other$seurat_clusters)])
Idents(sub_other) <- sub_other$Celltypes_refined

p_sub_cluster <- DimPlot(
  sub_other,
  reduction = "umap_sub",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 4,
  raster = FALSE
) + ggtitle(paste0("Sub-UMAP clusters: ", target_label))
safe_ggsave(
  file.path(annotation_figure_dir, sprintf("Sub_%s_umap_clusters.pdf", target_label)),
  p_sub_cluster,
  width = 8,
  height = 7
)

p_sub_final <- DimPlot(
  sub_other,
  reduction = "umap_sub",
  group.by = "Celltypes_refined",
  label = TRUE,
  label.size = 4,
  raster = FALSE
) + ggtitle(paste0("Final refined labels: ", target_label))
safe_ggsave(
  file.path(annotation_figure_dir, sprintf("Sub_%s_umap_final_refined.pdf", target_label)),
  p_sub_final,
  width = 8,
  height = 7
)

subcluster_sizes <- as.data.frame(table(sub_other$seurat_clusters))
colnames(subcluster_sizes) <- c("subcluster", "n_cells")
subcluster_sizes$manual_refined_celltype <- manual_refinement[as.character(subcluster_sizes$subcluster)]
readr::write_csv(
  subcluster_sizes,
  file.path(subcluster_dir, sprintf("Sub_%s_cluster_sizes_and_refinement.csv", target_label))
)

celllevel_refinement <- data.frame(
  cell = colnames(sub_other),
  parent_label = target_label,
  subcluster = as.character(sub_other$seurat_clusters),
  refined_celltype = as.character(sub_other$Celltypes_refined),
  stringsAsFactors = FALSE
)
readr::write_csv(
  celllevel_refinement,
  file.path(subcluster_dir, sprintf("Sub_%s_celllevel_refinement.csv", target_label))
)

# Write the local refinement result back to the global object.
scobj$Celltypes_refined <- as.character(scobj$Celltypes_refined)
scobj@meta.data[colnames(sub_other), "Celltypes_refined"] <- as.character(sub_other$Celltypes_refined)

# Global overlay for checking the refined subset in the original UMAP.
umap_global <- get_umap_name(scobj, dim_select)
scobj$Other_mix_refinement_in_context <- as.character(scobj$Celltypes)
refined_tags <- paste0(
  target_label,
  "_",
  as.character(sub_other$seurat_clusters),
  "_",
  as.character(sub_other$Celltypes_refined)
)
names(refined_tags) <- colnames(sub_other)
scobj@meta.data[names(refined_tags), "Other_mix_refinement_in_context"] <- refined_tags

p_global_context <- DimPlot(
  scobj,
  reduction = umap_global,
  group.by = "Other_mix_refinement_in_context",
  label = TRUE,
  label.size = 3,
  raster = FALSE
) + ggtitle("Other_mix subclusters in global UMAP context") +
  theme(legend.text = element_text(size = 7))
safe_ggsave(
  file.path(annotation_figure_dir, sprintf("Sub_%s_global_context.pdf", target_label)),
  p_global_context,
  width = 12,
  height = 9
)

saveRDS(
  sub_other,
  file = file.path(subcluster_dir, sprintf("Sub_%s_refined_object.rds", target_label))
)

# Finalize cell-type labels after refinement.
if (any(as.character(scobj$Celltypes_refined) == "Other_mix")) {
  stop("Other_mix labels remain after refinement; final annotation is incomplete.")
}
validate_final_labels(scobj$Celltypes_refined, final_levels, label_name = "Celltypes_refined")
scobj$Celltypes_refined <- factor(as.character(scobj$Celltypes_refined), levels = final_levels)
scobj$Celltypes <- scobj$Celltypes_refined
scobj$Celltype  <- scobj$Celltypes
Idents(scobj) <- scobj$Celltypes

if (any(is.na(scobj$Celltypes))) {
  stop("NA values detected in final Celltypes after factor conversion.")
}

# ------------------------------------------------------------------------------
# 7) Export final annotation summaries and basic figures
# ------------------------------------------------------------------------------
celltype_counts <- as.data.frame(table(scobj$Celltypes))
colnames(celltype_counts) <- c("Celltypes", "n_cells")
celltype_counts$proportion <- celltype_counts$n_cells / sum(celltype_counts$n_cells)

# Keep the output order fixed even when a category has zero cells.
celltype_counts$Celltypes <- factor(as.character(celltype_counts$Celltypes), levels = final_levels)
celltype_counts <- celltype_counts %>% arrange(Celltypes)

readr::write_csv(
  celltype_counts,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_final_celltype_counts.csv", dim_select, res_label))
)

message("Final observed cell-type counts:")
print(celltype_counts)

if (sum(celltype_counts$n_cells) != ncol(scobj)) {
  stop(
    "Final cell-type counts do not sum to the total number of cells in the object. ",
    "Please inspect Celltypes and Celltypes_refined."
  )
}

# Final cell-level annotation table.
cluster_annotation_table_final <- data.frame(
  cell = colnames(scobj),
  sample = if ("sample" %in% colnames(scobj@meta.data)) as.character(scobj$sample) else NA_character_,
  cluster = as.character(scobj$seurat_clusters),
  Celltypes_before_refinement = as.character(scobj$original_celltypes),
  Celltypes_final = as.character(scobj$Celltypes),
  stringsAsFactors = FALSE
)
readr::write_csv(
  cluster_annotation_table_final,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_celllevel_annotation_final.csv", dim_select, res_label))
)

if ("Group" %in% colnames(scobj@meta.data)) {
  group_by_celltype <- as.data.frame(table(Celltypes = scobj$Celltypes, Group = scobj$Group))
  colnames(group_by_celltype) <- c("Celltypes", "Group", "n_cells")
  readr::write_csv(
    group_by_celltype,
    file.path(annotation_table_dir, sprintf("dim_%d_res_%s_final_celltype_by_group.csv", dim_select, res_label))
  )
}

if ("sample" %in% colnames(scobj@meta.data)) {
  sample_by_celltype <- as.data.frame(table(sample = scobj$sample, Celltypes = scobj$Celltypes))
  colnames(sample_by_celltype) <- c("sample", "Celltypes", "n_cells")
  readr::write_csv(
    sample_by_celltype,
    file.path(annotation_table_dir, sprintf("dim_%d_res_%s_final_celltype_by_sample.csv", dim_select, res_label))
  )
}

plot_basic_umap(
  scobj,
  group_col = "Celltypes",
  reduction = umap_name,
  filename = file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_final_Celltypes_umap.pdf", dim_select, res_label)),
  width = 8,
  height = 8
)

if ("Group" %in% colnames(scobj@meta.data)) {
  plot_basic_umap(
    scobj,
    group_col = "Group",
    reduction = umap_name,
    filename = file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_Group_umap.pdf", dim_select, res_label)),
    width = 8,
    height = 8
  )
}

# Compact marker DotPlot for annotation QC.
marker_panels <- list(
  Ast      = c("GFAP", "AQP4", "SLC1A2", "GLUL", "ATP1B2"),
  ExcitNeu = c("SLC17A6", "SLC17A7", "CAMK2A", "SATB2", "TBR1", "CUX2"),
  InhibNeu = c("GAD1", "GAD2", "SLC32A1", "VIP", "RELN"),
  OL       = c("MBP", "PLP1", "MOG", "MAG", "MOBP")
)
marker_genes <- unique(unlist(lapply(marker_panels, match_genes_case_insensitive, available_genes = rownames(scobj))))
if (length(marker_genes) > 0L) {
  DefaultAssay(scobj) <- "RNA"
  p_dot <- Seurat::DotPlot(
    scobj,
    features = marker_genes,
    group.by = "Celltypes",
    scale = TRUE
  ) +
    ggplot2::scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold.italic"),
      axis.title = ggplot2::element_blank()
    )
  safe_ggsave(
    file.path(annotation_figure_dir, sprintf("dim_%d_res_%s_final_marker_dotplot.pdf", dim_select, res_label)),
    p_dot,
    width = 12,
    height = 8
  )
}

# ------------------------------------------------------------------------------
# 8) Save final annotated object and reproducibility logs
# ------------------------------------------------------------------------------
saveRDS(scobj, file = annotated3_file)
message("Final annotated object saved to: ", annotated3_file)

metadata_export <- scobj@meta.data %>%
  tibble::rownames_to_column("cell")
readr::write_csv(
  metadata_export,
  file.path(annotation_table_dir, sprintf("dim_%d_res_%s_final_metadata.csv", dim_select, res_label))
)

writeLines(
  c(
    paste0("Input object: ", normalizePath(obj_file, winslash = "/", mustWork = FALSE)),
    paste0("Input marker file: ", normalizePath(de_file, winslash = "/", mustWork = FALSE)),
    paste0("Output annotated2 before Other_mix refinement: ", normalizePath(annotated2_file, winslash = "/", mustWork = FALSE)),
    paste0("Output annotated3 after Other_mix refinement: ", normalizePath(annotated3_file, winslash = "/", mustWork = FALSE)),
    paste0("UMAP reduction: ", umap_name),
    paste0("Cluster column: ", cluster_col),
    "Other_mix refinement rule: subcluster 8 = Ast; subclusters 0-7 and 9-11 = Other",
    "SCI V4 note: final counts are observed outputs and are not hard-coded as expected values."
  ),
  con = file.path(annotation_log_dir, sprintf("dim_%d_res_%s_annotation_exports_SCIV4.txt", dim_select, res_label))
)

sink(file.path(annotation_log_dir, sprintf("sessionInfo_step3_annotation_human_SCIV4_dim%d_res%s.txt", dim_select, res_label)))
print(sessionInfo())
sink()

message("\n=== Human cell-type annotation SCI V4 finished successfully. ===")
message("Pre-refinement annotated object: ", annotated2_file)
message("Final annotated object: ", annotated3_file)
message("Annotation tables: ", annotation_table_dir)
message("Annotation figures: ", annotation_figure_dir)
