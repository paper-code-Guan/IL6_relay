# =============================================================================
# Script: Mouse single-nucleus RNA-seq quality control and sample merging
# Project: Mouse mPFC snRNA-seq analysis
# Description:
#   This script reproduces per-sample QC for mouse snRNA-seq data using Seurat
#   and DoubletFinder. Each sample is processed independently, including QC
#   metric calculation, cell/gene filtering, cell-cycle scoring, automatic pK
#   selection by BCmetric, DoubletFinder-based doublet removal, and final merging
#   into control, KO, and all-sample Seurat objects.
#
# Output filenames:
#       1-QC/<sample>/<sample>_qc_final.rds
#       con_group_merged.rds
#       ko_group_merged.rds
#       all_samples_merged.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(DoubletFinder)
  library(stringr)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# =============================================================================
# 0. Input/output paths and global parameters
# =============================================================================

samples_dir <- Sys.getenv("MOUSE_INPUT_DIR", unset = "")
reproduce_dir <- Sys.getenv("MOUSE_OUTPUT_DIR", unset = "")
if (!nzchar(samples_dir)) {
  stop("MOUSE_INPUT_DIR must be set.", call. = FALSE)
}
if (!nzchar(reproduce_dir)) {
  stop("MOUSE_OUTPUT_DIR must be set.", call. = FALSE)
}

qc_root <- file.path(reproduce_dir, "1-QC")

if (!dir.exists(reproduce_dir)) dir.create(reproduce_dir, recursive = TRUE)
if (!dir.exists(qc_root))      dir.create(qc_root, recursive = TRUE)

cat("[INFO] Output directory:", reproduce_dir, "\n")

# Per-sample QC thresholds. These values are retained from the validated workflow.
qc_params <- list(
  "con1" = list(nFeature_min = 300, nFeature_max = 5500, percent_mito = 10,  percent_ribo = 2,   percent_hb = 0.5),
  "con2" = list(nFeature_min = 300, nFeature_max = 6000, percent_mito = 2.5, percent_ribo = 0.6, percent_hb = 0.05),
  "con3" = list(nFeature_min = 300, nFeature_max = 5500, percent_mito = 10,  percent_ribo = 2,   percent_hb = 0.5),
  "ko1"  = list(nFeature_min = 300, nFeature_max = 5500, percent_mito = 10,  percent_ribo = 2,   percent_hb = 0.5),
  "ko2"  = list(nFeature_min = 300, nFeature_max = 5500, percent_mito = 10,  percent_ribo = 2,   percent_hb = 0.5),
  "ko3"  = list(nFeature_min = 300, nFeature_max = 5500, percent_mito = 10,  percent_ribo = 2,   percent_hb = 0.5)
)

# QC gene patterns are kept unchanged from the validated workflow.
mito_pattern <- "^mt-"
ribo_pattern <- "^Rp[sl]"
hb_pattern   <- "^Hb[^(p)]"

if (!dir.exists(samples_dir)) stop("Input samples directory does not exist: ", samples_dir)

samples <- list.dirs(samples_dir, full.names = FALSE, recursive = FALSE)
samples <- samples[nzchar(samples)]
if (length(samples) == 0L) stop("No samples found in: ", samples_dir)

samples_simple <- gsub("^sample\\d+_", "", samples)
if (anyDuplicated(samples_simple)) {
  stop("Duplicated simplified sample names detected: ",
       paste(unique(samples_simple[duplicated(samples_simple)]), collapse = ", "))
}

cat("[INFO] Samples detected (n = ", length(samples_simple), "):\n", sep = "")
print(samples_simple)

# =============================================================================
# 1. Helper functions
# =============================================================================

get_assay_counts <- function(seurat_obj, assay = "RNA") {
  Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
}

save_pdf <- function(plot, filename, width, height, units = "cm") {
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    units    = units,
    device   = grDevices::cairo_pdf,
    limitsize = FALSE
  )
}

generate_qc_plots <- function(seurat_obj, prefix, outdir) {
  cat("[PLOT] QC plots -", prefix, "\n")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  p1 <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA"),
    pt.size  = 0,
    ncol     = 2
  ) + NoLegend()
  save_pdf(plot = p1,
           filename = file.path(outdir, paste0("Vlnplot_features_", prefix, ".pdf")),
           width = 16, height = 12, units = "cm")

  p2 <- VlnPlot(
    seurat_obj,
    features = c("percent_mito", "percent_ribo", "percent_hb"),
    pt.size  = 0,
    ncol     = 3
  ) + NoLegend()
  save_pdf(plot = p2,
           filename = file.path(outdir, paste0("Vlnplot_metrics_", prefix, ".pdf")),
           width = 24, height = 12, units = "cm")

  p3 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    pt.size  = 0.5
  )
  save_pdf(plot = p3,
           filename = file.path(outdir, paste0("Scatterplot_", prefix, ".pdf")),
           width = 12, height = 9, units = "cm")
}

get_doubletfinder_function <- function(name_v3, name_current) {
  ns <- asNamespace("DoubletFinder")
  if (exists(name_v3, envir = ns, inherits = FALSE)) {
    return(get(name_v3, envir = ns))
  }
  if (exists(name_current, envir = ns, inherits = FALSE)) {
    return(get(name_current, envir = ns))
  }
  stop("Neither ", name_v3, " nor ", name_current, " was found in DoubletFinder.")
}

process_for_doublet_detection <- function(seurat_obj, dims_use = 1:30, cluster_resolution = 0.6) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_use, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution, verbose = FALSE)
  seurat_obj
}

detect_doublets <- function(seurat_obj, dims_use = 1:30, pN = 0.25) {
  cat("[DoubletFinder] Running automatic pK selection by BCmetric\n")

  param_sweep_fun  <- get_doubletfinder_function("paramSweep_v3", "paramSweep")
  summarize_fun    <- get_doubletfinder_function("summarizeSweep", "summarizeSweep")
  find_pk_fun      <- get_doubletfinder_function("find.pK", "find.pK")
  doublet_find_fun <- get_doubletfinder_function("doubletFinder_v3", "doubletFinder")

  sweep_list  <- param_sweep_fun(seurat_obj, PCs = dims_use, sct = FALSE)
  sweep_stats <- summarize_fun(sweep_list, GT = FALSE)
  bcmvn       <- find_pk_fun(sweep_stats)

  if (!all(c("pK", "BCmetric") %in% colnames(bcmvn))) {
    stop("find.pK output does not contain required columns: pK and BCmetric.")
  }

  pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  if (!is.finite(pK_bcmvn)) stop("Automatic pK selection failed.")

  # Expected doublet rate retained from the validated workflow.
  doublet_rate <- ncol(seurat_obj) * 8 * 1e-6
  annotations  <- seurat_obj@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  nExp_poi     <- round(doublet_rate * ncol(seurat_obj))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

  cat("[DoubletFinder] pK:", pK_bcmvn,
      "| doublet_rate:", doublet_rate,
      "| nExp_poi:", nExp_poi,
      "| nExp_poi_adj:", nExp_poi_adj, "\n")

  seurat_obj <- doublet_find_fun(
    seurat_obj,
    PCs       = dims_use,
    pN        = pN,
    pK        = pK_bcmvn,
    nExp      = nExp_poi_adj,
    reuse.pANN = FALSE,
    sct       = FALSE
  )

  df_cols <- grep("^DF\\.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(df_cols) == 0L) stop("DoubletFinder classification column was not found.")

  df_class_col <- tail(df_cols, 1)
  seurat_obj$DoubletFinder <- seurat_obj@meta.data[[df_class_col]]

  doublet_cells <- colnames(seurat_obj)[seurat_obj$DoubletFinder == "Doublet"]
  cat("[DoubletFinder] Doublets detected:", length(doublet_cells), "\n")

  if (length(doublet_cells) > 0L) {
    seurat_obj <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), doublet_cells))
    cat("[DoubletFinder] Dimensions after doublet removal:", paste(dim(seurat_obj), collapse = " x "), "\n")
  }

  attr(seurat_obj, "doubletfinder_qc") <- list(
    pK = pK_bcmvn,
    doublet_rate = doublet_rate,
    homotypic_prop = homotypic_prop,
    nExp_poi = nExp_poi,
    nExp_poi_adj = nExp_poi_adj,
    classification_column = df_class_col,
    n_doublets = length(doublet_cells)
  )

  seurat_obj
}

# =============================================================================
# 2. Per-sample QC
# =============================================================================

perform_sample_qc <- function(sample_index) {
  sample_name <- samples_simple[sample_index]
  sample_dir_raw <- samples[sample_index]

  cat("\n============================================================\n")
  cat("[SAMPLE] Processing:", sample_name, "\n")
  cat("============================================================\n")

  params <- qc_params[[sample_name]]
  if (is.null(params)) {
    stop("No QC parameters were defined for sample: ", sample_name)
  }

  qc_dir <- file.path(qc_root, sample_name)
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

  data_path <- file.path(samples_dir, sample_dir_raw, "filtered_cell_gene_matrix")
  if (!dir.exists(data_path)) stop("Input data path does not exist: ", data_path)

  cat("[LOAD]", data_path, "\n")
  counts <- Read10X(data.dir = data_path)

  sce <- CreateSeuratObject(
    counts       = counts,
    project      = sample_name,
    min.cells    = 0,
    min.features = 0
  )
  rm(counts)
  gc()

  initial_genes <- nrow(sce)
  initial_cells <- ncol(sce)
  cat("[QC] Initial dimensions:", paste(dim(sce), collapse = " x "), "\n")

  sce <- PercentageFeatureSet(sce, pattern = mito_pattern, col.name = "percent_mito")
  sce <- PercentageFeatureSet(sce, pattern = ribo_pattern, col.name = "percent_ribo")
  sce <- PercentageFeatureSet(sce, pattern = hb_pattern,   col.name = "percent_hb")

  generate_qc_plots(sce, "before_filtering", qc_dir)

  cat("[QC] Thresholds:\n")
  print(params)

  keep_cells_nfeature <- WhichCells(
    sce,
    expression = nFeature_RNA > params$nFeature_min & nFeature_RNA < params$nFeature_max
  )

  counts_mat <- get_assay_counts(sce)
  keep_features <- rownames(sce)[Matrix::rowSums(counts_mat > 0) > 3]

  if (!is.null(params$umi_max)) {
    gene_total_counts <- Matrix::rowSums(counts_mat)
    keep_features <- rownames(sce)[gene_total_counts > 3 & gene_total_counts < params$umi_max]
  }

  sce <- subset(sce, features = keep_features, cells = keep_cells_nfeature)
  cat("[QC] After nFeature/gene filtering:", paste(dim(sce), collapse = " x "), "\n")

  keep_mito <- WhichCells(sce, expression = percent_mito < params$percent_mito)
  keep_ribo <- WhichCells(sce, expression = percent_ribo < params$percent_ribo)
  keep_hb   <- WhichCells(sce, expression = percent_hb   < params$percent_hb)

  final_cells <- Reduce(intersect, list(keep_mito, keep_ribo, keep_hb))
  sce <- subset(sce, cells = final_cells)
  cat("[QC] After mito/ribo/HB filtering:", paste(dim(sce), collapse = " x "), "\n")

  generate_qc_plots(sce, "after_filtering", qc_dir)

  cat("[QC] Cell-cycle scoring\n")
  sce <- NormalizeData(sce, verbose = FALSE)
  s_genes   <- stringr::str_to_title(tolower(cc.genes.updated.2019$s.genes))
  g2m_genes <- stringr::str_to_title(tolower(cc.genes.updated.2019$g2m.genes))
  s_genes   <- intersect(s_genes, rownames(sce))
  g2m_genes <- intersect(g2m_genes, rownames(sce))

  sce <- CellCycleScoring(
    object       = sce,
    s.features   = s_genes,
    g2m.features = g2m_genes
  )

  cat("[QC] Doublet detection\n")
  sce <- process_for_doublet_detection(sce, dims_use = 1:30, cluster_resolution = 0.6)
  sce <- detect_doublets(sce, dims_use = 1:30, pN = 0.25)

  final_genes <- nrow(sce)
  final_cells_n <- ncol(sce)
  df_qc <- attr(sce, "doubletfinder_qc")

  output_file <- file.path(qc_dir, paste0(sample_name, "_qc_final.rds"))
  saveRDS(sce, output_file)
  cat("[SAVE]", output_file, "\n")

  summary_row <- data.frame(
    sample = sample_name,
    raw_sample_dir = sample_dir_raw,
    input_path = data_path,
    initial_genes = initial_genes,
    initial_cells = initial_cells,
    final_genes = final_genes,
    final_cells = final_cells_n,
    pK = if (!is.null(df_qc)) df_qc$pK else NA_real_,
    doublet_rate = if (!is.null(df_qc)) df_qc$doublet_rate else NA_real_,
    homotypic_prop = if (!is.null(df_qc)) df_qc$homotypic_prop else NA_real_,
    nExp_poi = if (!is.null(df_qc)) df_qc$nExp_poi else NA_integer_,
    nExp_poi_adj = if (!is.null(df_qc)) df_qc$nExp_poi_adj else NA_integer_,
    n_doublets = if (!is.null(df_qc)) df_qc$n_doublets else NA_integer_,
    output_file = output_file,
    stringsAsFactors = FALSE
  )

  summary_file <- file.path(qc_dir, paste0(sample_name, "_qc_summary.csv"))
  write.csv(summary_row, summary_file, row.names = FALSE)
  cat("[SAVE]", summary_file, "\n")

  invisible(summary_row)
}

# =============================================================================
# 3. Merge QC-filtered samples
# =============================================================================

merge_samples <- function(sample_names, output_name) {
  cat("\n============================================================\n")
  cat("[MERGE]", output_name, ": ", paste(sample_names, collapse = ", "), "\n", sep = "")
  cat("============================================================\n")

  if (length(sample_names) == 0L) stop("No samples provided for merging: ", output_name)

  sce_list <- lapply(sample_names, function(sample) {
    file_path <- file.path(qc_root, sample, paste0(sample, "_qc_final.rds"))
    if (!file.exists(file_path)) stop("QC object not found: ", file_path)
    readRDS(file_path)
  })
  names(sce_list) <- sample_names

  merged_sce <- if (length(sce_list) == 1L) {
    sce_list[[1]]
  } else {
    merge(x = sce_list[[1]], y = sce_list[-1], add.cell.ids = sample_names)
  }

  output_file <- file.path(reproduce_dir, paste0(output_name, "_merged.rds"))
  saveRDS(merged_sce, output_file)
  cat("[SAVE]", output_file, "| dimensions:", paste(dim(merged_sce), collapse = " x "), "\n")

  invisible(list(object = merged_sce, output_file = output_file))
}

# =============================================================================
# 4. Main execution
# =============================================================================

cat("\n[MAIN] Starting mouse QC reproduction\n")
cat("[MAIN] Input samples directory:", samples_dir, "\n")
cat("[MAIN] Reproduction directory:", reproduce_dir, "\n")

sample_summary_list <- vector("list", length(samples_simple))
names(sample_summary_list) <- samples_simple

for (i in seq_along(samples_simple)) {
  sample_summary_list[[i]] <- perform_sample_qc(i)
  gc()
}

sample_summary <- dplyr::bind_rows(sample_summary_list)
summary_out <- file.path(qc_root, "mouse_qc_summary_all_samples.csv")
write.csv(sample_summary, summary_out, row.names = FALSE)
cat("[SAVE]", summary_out, "\n")

con_samples <- grep("^con", samples_simple, value = TRUE)
ko_samples  <- grep("^ko",  samples_simple, value = TRUE)

con_merged <- merge_samples(con_samples, "con_group")
ko_merged  <- merge_samples(ko_samples,  "ko_group")
all_merged <- merge_samples(samples_simple, "all_samples")

cat("\n==================== QC REPRODUCTION SUMMARY ====================\n")
cat("Reproduction directory:", reproduce_dir, "\n")
cat("Total samples processed:", length(samples_simple), "\n")
cat("Control samples:", length(con_samples), "-", paste(con_samples, collapse = ", "), "\n")
cat("KO samples:", length(ko_samples), "-", paste(ko_samples, collapse = ", "), "\n")
cat("Merged datasets:\n")
cat("  con_group_merged.rds:", paste(dim(con_merged$object), collapse = " x "), "\n")
cat("  ko_group_merged.rds :", paste(dim(ko_merged$object), collapse = " x "), "\n")
cat("  all_samples_merged.rds:", paste(dim(all_merged$object), collapse = " x "), "\n")
cat("Sample-level QC summary:", summary_out, "\n")
cat("=================================================================\n")

session_file <- file.path(qc_root, "sessionInfo_step1_qc_mouse.txt")
sink(session_file)
print(sessionInfo())
sink()
cat("[SAVE]", session_file, "\n")
cat("[DONE] Mouse QC reproduction completed successfully.\n")
