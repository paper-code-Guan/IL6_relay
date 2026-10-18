# =============================================================================
# Script: Step 1 - Human dlPFC snRNA-seq quality control and sample merging
# Project: Human schizophrenia dlPFC snRNA-seq reanalysis
# Purpose:
#   1) Build per-sample Seurat objects from the released raw count matrices.
#   2) Compute standard snRNA-seq QC metrics.
#   3) Apply per-sample QC filters and DoubletFinder-based doublet removal.
#   4) Merge QC-filtered samples and attach sample-level metadata.
#   5) Export post-QC visualization files and summary tables.
#
# Notes:
#   - Paths and key thresholds are kept consistent with the validated workflow.
#   - The script is designed for Linux/RStudio Server or terminal Rscript usage.
#   - Analyses start from released count matrices, not raw FASTQ files.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(DoubletFinder)
  library(readr)
  library(tibble)
  library(stringr)
})

set.seed(123)
options(stringsAsFactors = FALSE)

# =============================================================================
# 0. Input/output paths and global parameters
# =============================================================================

required_env_path <- function(name, must_exist = TRUE) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop(name, " must be set.", call. = FALSE)
  }
  normalizePath(value, winslash = "/", mustWork = must_exist)
}

# Named list of dgCMatrix count matrices, one matrix per released sample.
counts_list_rds <- required_env_path("BATIUK_COUNTS_RDS")

# Optional named annotation list keyed by cell barcode.
anno_rds <- Sys.getenv("BATIUK_ANNOTATION_RDS", unset = "")
if (nzchar(anno_rds)) {
  anno_rds <- normalizePath(anno_rds, winslash = "/", mustWork = TRUE)
}

sample_manifest_path <- required_env_path("BATIUK_SAMPLE_MANIFEST")
project_name  <- "human_brain"
reproduce_dir <- required_env_path("BATIUK_OUTPUT_DIR", must_exist = FALSE)
qc_root       <- file.path(reproduce_dir, "1-QC")
qc_plot_dir   <- file.path(reproduce_dir, "qc_plots_after_filter")
log_dir       <- file.path(reproduce_dir, "logs")

for (d in c(reproduce_dir, qc_root, qc_plot_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Per-sample QC thresholds. A sample-specific override can be added later by
# adding an entry named by sample ID. The default thresholds are retained from
# the validated human workflow.
qc_params <- list(
  "default" = list(
    nFeature_min = 300,
    nFeature_max = 10000,
    percent_mito = 25,
    percent_ribo = 2,
    percent_hb   = 0.5
  )
)

# Curated human hemoglobin genes. This avoids matching unrelated genes such as HBP1.
RBC_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

# Important genes retained even when lowly expressed in a sample. This prevents
# downstream IL6/IL-6 pathway analyses from failing because of sample-level gene filtering.
important_genes <- c("IL6", "SOCS1", "SOCS3", "IRF9", "EBI3")

# DoubletFinder settings retained from the validated workflow.
df_pN          <- 0.25
df_dims_pk     <- 1:20
df_dims_full   <- 1:30
df_nfeatures   <- 2000
df_resolution  <- 0.5
max_pk_cells   <- 15000L

cat("[INFO] Reproduce directory: ", reproduce_dir, "\n", sep = "")
cat("[INFO] QC root: ", qc_root, "\n", sep = "")

# =============================================================================
# 1. Helper functions
# =============================================================================

write_session_info <- function(path) {
  sink(path)
  print(sessionInfo())
  sink()
}

get_doubletfinder_function <- function(fun_v3, fun_plain) {
  if (exists(fun_v3, mode = "function")) {
    return(get(fun_v3, mode = "function"))
  }
  if (exists(fun_plain, mode = "function")) {
    return(get(fun_plain, mode = "function"))
  }
  stop("Cannot find DoubletFinder function: ", fun_v3, " or ", fun_plain)
}

summarize_qc_metrics <- function(obj, stage, sample_id) {
  md <- obj@meta.data
  data.frame(
    sample = sample_id,
    stage = stage,
    n_genes = nrow(obj),
    n_cells = ncol(obj),
    median_nCount_RNA = median(md$nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(md$nFeature_RNA, na.rm = TRUE),
    median_percent_mito = if ("percent_mito" %in% colnames(md)) median(md$percent_mito, na.rm = TRUE) else NA_real_,
    median_percent_ribo = if ("percent_ribo" %in% colnames(md)) median(md$percent_ribo, na.rm = TRUE) else NA_real_,
    median_percent_hb = if ("percent_hb" %in% colnames(md)) median(md$percent_hb, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

generate_qc_plots <- function(seurat_obj, prefix, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  p1 <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA"),
    pt.size = 0,
    ncol = 2
  ) + NoLegend()
  ggsave(
    file.path(outdir, paste0("Vlnplot_features_", prefix, ".pdf")),
    plot = p1,
    width = 16,
    height = 12,
    units = "cm",
    device = cairo_pdf
  )

  metric_features <- intersect(c("percent_mito", "percent_ribo", "percent_hb"), colnames(seurat_obj[[]]))
  if (length(metric_features) > 0) {
    p2 <- VlnPlot(
      seurat_obj,
      features = metric_features,
      pt.size = 0,
      ncol = length(metric_features)
    ) + NoLegend()
    ggsave(
      file.path(outdir, paste0("Vlnplot_metrics_", prefix, ".pdf")),
      plot = p2,
      width = 24,
      height = 12,
      units = "cm",
      device = cairo_pdf
    )
  }

  p3 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    pt.size = 0.5
  )
  ggsave(
    file.path(outdir, paste0("Scatterplot_", prefix, ".pdf")),
    plot = p3,
    width = 12,
    height = 9,
    units = "cm",
    device = cairo_pdf
  )
}

add_annotations_to_object <- function(sobj, annos) {
  if (!is.list(annos) || length(annos) == 0) {
    return(sobj)
  }
  for (nm in names(annos)) {
    vec <- annos[[nm]]
    if (!is.null(names(vec))) {
      sobj[[nm]] <- vec[colnames(sobj)]
    }
  }
  sobj
}

run_doubletfinder_memory <- function(obj,
                                     sample_id,
                                     pN = df_pN,
                                     dims_pk = df_dims_pk,
                                     dims_df = df_dims_full,
                                     max_pk_cells_use = max_pk_cells) {
  cat("[DF] Preparing DoubletFinder for sample: ", sample_id, "\n", sep = "")

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = df_nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = max(dims_df), verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_df, verbose = FALSE)
  obj <- FindClusters(obj, resolution = df_resolution, verbose = FALSE)

  if (ncol(obj) > max_pk_cells_use) {
    cat("[DF] Downsampling to ", max_pk_cells_use, " cells for pK sweep.\n", sep = "")
    cells_ds <- colnames(obj)[sample.int(ncol(obj), max_pk_cells_use)]
    obj_ds <- subset(obj, cells = cells_ds)
  } else {
    obj_ds <- obj
  }

  param_sweep_fun <- get_doubletfinder_function("paramSweep_v3", "paramSweep")
  doublet_fun     <- get_doubletfinder_function("doubletFinder_v3", "doubletFinder")

  cat("[DF] Running pK sweep on ", ncol(obj_ds), " cells.\n", sep = "")
  sweep_list  <- param_sweep_fun(obj_ds, PCs = dims_pk, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
  bcmvn       <- find.pK(sweep_stats)

  if (is.null(bcmvn) || nrow(bcmvn) == 0) {
    pK_bcmvn <- 0.05
    cat("[DF] pK table is empty; fallback pK = ", pK_bcmvn, "\n", sep = "")
  } else {
    pK_bcmvn <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    if (is.na(pK_bcmvn)) {
      pK_bcmvn <- 0.05
    }
    cat("[DF] Selected pK by BCmetric: ", pK_bcmvn, "\n", sep = "")
  }

  # Retained from the validated workflow.
  doublet_rate    <- ncol(obj) * 8e-6
  homotypic_prop  <- modelHomotypic(obj$seurat_clusters)
  nExp_poi        <- round(doublet_rate * ncol(obj))
  nExp_poi_adj    <- round(nExp_poi * (1 - homotypic_prop))

  cat(
    "[DF] DoubletRate = ", round(doublet_rate, 6),
    "; nExp = ", nExp_poi,
    "; nExp adjusted = ", nExp_poi_adj, "\n",
    sep = ""
  )

  obj <- doublet_fun(
    obj,
    PCs = dims_df,
    pN = pN,
    pK = pK_bcmvn,
    nExp = nExp_poi_adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  df_cols <- grep("^DF\\.classifications", colnames(obj@meta.data), value = TRUE)
  if (length(df_cols) < 1) {
    stop("DoubletFinder classification column was not generated for sample: ", sample_id)
  }
  obj$DoubletFinder <- obj@meta.data[[tail(df_cols, 1)]]

  n_before <- ncol(obj)
  obj <- obj[, obj$DoubletFinder != "Doublet"]
  n_removed <- n_before - ncol(obj)

  cat("[DF] Removed ", n_removed, " doublets. Kept ", ncol(obj), " cells.\n", sep = "")

  list(
    object = obj,
    summary = data.frame(
      sample = sample_id,
      pN = pN,
      selected_pK = pK_bcmvn,
      doublet_rate = doublet_rate,
      nExp_poi = nExp_poi,
      nExp_poi_adj = nExp_poi_adj,
      doublets_removed = n_removed,
      cells_before_doublet_removal = n_before,
      cells_after_doublet_removal = ncol(obj),
      stringsAsFactors = FALSE
    )
  )
}

qc_one_sample <- function(sample_id, counts, annos = NULL) {
  cat("\n================ SAMPLE: ", sample_id, " ================\n", sep = "")

  outdir <- file.path(qc_root, sample_id)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  sobj <- CreateSeuratObject(
    counts = counts,
    project = project_name,
    min.cells = 0,
    min.features = 0
  )
  sobj$orig.ident <- sample_id
  sobj$sample <- sample_id

  cat("[QC] Initial dimensions: ", paste(dim(sobj), collapse = " x "), "\n", sep = "")

  sobj <- add_annotations_to_object(sobj, annos)

  sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent_mito")
  sobj <- PercentageFeatureSet(sobj, pattern = "^RP[SL]", col.name = "percent_ribo")
  sobj <- PercentageFeatureSet(
    sobj,
    features = intersect(RBC_GENES, rownames(sobj)),
    col.name = "percent_hb"
  )

  generate_qc_plots(sobj, "before_filtering", outdir)

  params <- qc_params[[sample_id]]
  if (is.null(params)) {
    params <- qc_params[["default"]]
  }

  qc_summary <- list()
  qc_summary[["initial"]] <- summarize_qc_metrics(sobj, "initial", sample_id)

  keep_cells <- WhichCells(
    sobj,
    expression = nFeature_RNA > params$nFeature_min &
      nFeature_RNA < params$nFeature_max
  )

  expressed_in_cells <- Matrix::rowSums(sobj@assays$RNA@counts > 0)
  keep_genes <- unique(c(
    names(expressed_in_cells)[expressed_in_cells >= 1],
    intersect(important_genes, rownames(sobj))
  ))

  sobj <- subset(sobj, features = keep_genes, cells = keep_cells)
  if (!"IL6" %in% rownames(sobj)) {
    stop("IL6 is missing after sample-level feature filtering: ", sample_id)
  }
  cat("[QC] After nFeature filtering: ", paste(dim(sobj), collapse = " x "), "\n", sep = "")
  qc_summary[["after_nFeature"]] <- summarize_qc_metrics(sobj, "after_nFeature_filtering", sample_id)

  keep_mito <- WhichCells(sobj, expression = percent_mito < params$percent_mito)
  keep_ribo <- WhichCells(sobj, expression = percent_ribo < params$percent_ribo)
  keep_hb   <- WhichCells(sobj, expression = percent_hb < params$percent_hb)

  sobj <- subset(sobj, cells = Reduce(intersect, list(keep_mito, keep_ribo, keep_hb)))
  cat("[QC] After metric filtering: ", paste(dim(sobj), collapse = " x "), "\n", sep = "")
  qc_summary[["after_metrics"]] <- summarize_qc_metrics(sobj, "after_metric_filtering", sample_id)

  generate_qc_plots(sobj, "after_filtering", outdir)

  cat("[QC] Cell-cycle scoring\n")
  sobj <- NormalizeData(sobj, verbose = FALSE)
  sobj <- CellCycleScoring(
    sobj,
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes
  )

  df_res <- run_doubletfinder_memory(
    sobj,
    sample_id = sample_id,
    pN = df_pN,
    dims_pk = df_dims_pk,
    dims_df = df_dims_full
  )
  sobj <- df_res$object

  qc_summary[["after_doublet_removal"]] <- summarize_qc_metrics(sobj, "after_doublet_removal", sample_id)

  sample_qc_summary <- dplyr::bind_rows(qc_summary)
  write.csv(
    sample_qc_summary,
    file.path(outdir, paste0(sample_id, "_qc_summary_stages.csv")),
    row.names = FALSE
  )
  write.csv(
    df_res$summary,
    file.path(outdir, paste0(sample_id, "_doubletfinder_summary.csv")),
    row.names = FALSE
  )

  outfile <- file.path(outdir, paste0(sample_id, "_qc_final.rds"))
  saveRDS(sobj, outfile)
  cat("[SAVE] ", outfile, "\n", sep = "")

  list(
    sample = sample_id,
    final_rds = outfile,
    qc_summary = sample_qc_summary,
    doubletfinder_summary = df_res$summary
  )
}

merge_samples <- function(sample_ids, tag) {
  objs <- lapply(sample_ids, function(sid) {
    path <- file.path(qc_root, sid, paste0(sid, "_qc_final.rds"))
    if (!file.exists(path)) {
      stop("Missing per-sample QC RDS: ", path)
    }
    readRDS(path)
  })
  names(objs) <- sample_ids

  if (length(objs) == 1) {
    merged <- objs[[1]]
  } else {
    merged <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = sample_ids)
  }

  out <- file.path(reproduce_dir, paste0(tag, "_merged.rds"))
  saveRDS(merged, out)
  cat("[SAVE] ", out, " | dims: ", paste(dim(merged), collapse = " x "), "\n", sep = "")
  merged
}

attach_sample_metadata <- function(scobj, sample_manifest_path) {
  if (!file.exists(sample_manifest_path)) {
    warning("The Batiuk sample manifest was not found. The clean merged object will be saved without sample metadata.")
    return(scobj)
  }

  samdat <- readr::read_csv(sample_manifest_path, show_col_types = FALSE)
  if ("sample_id" %in% colnames(samdat) && !"Identifier" %in% colnames(samdat)) {
    samdat <- samdat |>
      dplyr::rename(Identifier = sample_id)
  }
  if ("diagnosis" %in% colnames(samdat) && !"Diagnosis" %in% colnames(samdat)) {
    samdat <- samdat |>
      dplyr::rename(Diagnosis = diagnosis)
  }
  if (!"Identifier" %in% colnames(samdat)) {
    stop("The Batiuk sample manifest must contain a 'sample_id' or 'Identifier' column.")
  }

  sample_vec <- if ("sample" %in% colnames(scobj@meta.data)) {
    as.character(scobj$sample)
  } else if ("orig.ident" %in% colnames(scobj@meta.data)) {
    as.character(scobj$orig.ident)
  } else {
    stop("Neither 'sample' nor 'orig.ident' exists in meta.data.")
  }

  not_found <- setdiff(unique(sample_vec), samdat$Identifier)
  if (length(not_found) > 0) {
    warning("Samples not found in the Batiuk sample manifest: ", paste(not_found, collapse = ", "))
  }

  samdat_reindexed <- samdat %>% column_to_rownames("Identifier")
  new_meta <- samdat_reindexed[sample_vec, , drop = FALSE]
  rownames(new_meta) <- colnames(scobj)

  duplicate_cols <- intersect(colnames(scobj@meta.data), colnames(new_meta))
  if (length(duplicate_cols) > 0) {
    colnames(new_meta)[match(duplicate_cols, colnames(new_meta))] <- paste0("samp.", duplicate_cols)
  }

  scobj <- AddMetaData(scobj, metadata = new_meta)

  meta <- scobj@meta.data
  df_cols <- grep("^DF\\.classifications_", colnames(meta), value = TRUE)
  if (length(df_cols) > 0) {
    meta$DoubletFinder_class <- apply(meta[, df_cols, drop = FALSE], 1, function(x) {
      x <- as.character(x)
      idx <- which(!is.na(x) & nzchar(x))
      if (length(idx) > 0) x[idx[1]] else NA_character_
    })
  }

  drop_cols <- c(
    grep("^pANN_", colnames(meta), value = TRUE),
    df_cols
  )
  if ("percent_hb" %in% colnames(meta) && all(meta$percent_hb == 0, na.rm = TRUE)) {
    drop_cols <- c(drop_cols, "percent_hb")
  }

  scobj@meta.data <- meta[, setdiff(colnames(meta), unique(drop_cols)), drop = FALSE]

  if ("Diagnosis" %in% colnames(scobj@meta.data)) {
    diagnosis <- as.character(scobj$Diagnosis)
  } else if ("samp.Diagnosis" %in% colnames(scobj@meta.data)) {
    diagnosis <- as.character(scobj$`samp.Diagnosis`)
  } else {
    diagnosis <- rep(NA_character_, ncol(scobj))
  }

  if (!all(is.na(diagnosis))) {
    ctrl_pattern <- "(?i)^(ctr|control|hc|healthy|ctrl)$|\\bcontrol\\b|^ctr"
    scobj$GROUP <- ifelse(stringr::str_detect(diagnosis, ctrl_pattern), "CTR", "SCZ")
    scobj$group <- scobj$GROUP
  }

  scobj
}

make_post_qc_visualizations <- function(scobj, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  assay_use <- DefaultAssay(scobj)
  counts <- GetAssayData(scobj, assay = assay_use, slot = "counts")

  add_pct_if_missing <- function(obj, mat, col_name, pattern_vec) {
    if (!col_name %in% colnames(obj@meta.data)) {
      idx <- unique(unlist(lapply(pattern_vec, function(p) grep(p, rownames(mat), ignore.case = TRUE))))
      if (length(idx) > 0) {
        total_counts <- Matrix::colSums(mat)
        sub_counts   <- Matrix::colSums(mat[idx, , drop = FALSE])
        obj[[col_name]] <- as.numeric(sub_counts / pmax(total_counts, 1) * 100)
        message("[QC] Computed ", col_name, " using ", length(idx), " genes.")
      } else {
        message("[QC] Skipped ", col_name, ": no matching genes found.")
      }
    }
    obj
  }

  scobj <- add_pct_if_missing(scobj, counts, "percent_mito", c("^MT-"))
  scobj <- add_pct_if_missing(scobj, counts, "percent_ribo", c("^RPL", "^RPS"))
  scobj <- add_pct_if_missing(scobj, counts, "percent_hb", c("^HBA", "^HBB", "^HBD", "^HBG", "^HBM", "^HBQ", "^HBZ"))

  if ("percent.mt" %in% colnames(scobj@meta.data) && !"percent_mito" %in% colnames(scobj@meta.data)) {
    scobj$percent_mito <- scobj$percent.mt
  }

  features_basic <- intersect(c("nFeature_RNA", "nCount_RNA"), colnames(scobj@meta.data))
  if (length(features_basic) > 0) {
    p_violin_basic <- VlnPlot(
      scobj,
      group.by = "orig.ident",
      features = features_basic,
      pt.size = 0,
      raster = FALSE,
      ncol = length(features_basic)
    ) + NoLegend()
    ggsave(
      file.path(outdir, "Vln_features_nCount_nFeature_bySample.pdf"),
      p_violin_basic,
      width = 24,
      height = 10,
      units = "cm",
      device = cairo_pdf
    )
  }

  features_metrics <- intersect(c("percent_mito", "percent_ribo", "percent_hb"), colnames(scobj@meta.data))
  if (length(features_metrics) > 0) {
    p_violin_metrics <- VlnPlot(
      scobj,
      group.by = "orig.ident",
      features = features_metrics,
      pt.size = 0,
      raster = FALSE,
      ncol = length(features_metrics)
    ) + NoLegend()
    ggsave(
      file.path(outdir, "Vln_metrics_mito_ribo_hb_bySample.pdf"),
      p_violin_metrics,
      width = 28,
      height = 12,
      units = "cm",
      device = cairo_pdf
    )
  }

  scatter_vars <- intersect(c("nCount_RNA", "nFeature_RNA", "percent_mito", "orig.ident"), colnames(scobj@meta.data))
  if (all(c("nCount_RNA", "nFeature_RNA", "orig.ident") %in% scatter_vars)) {
    df_scatter <- FetchData(scobj, vars = c("nCount_RNA", "nFeature_RNA", "percent_mito", "orig.ident")) %>%
      dplyr::rename(sample = orig.ident)

    p_scatter <- ggplot(df_scatter, aes(nCount_RNA, nFeature_RNA)) +
      geom_bin2d(bins = 60) +
      facet_wrap(~ sample, scales = "free") +
      xlab("nCount_RNA") +
      ylab("nFeature_RNA") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold"))

    ggsave(
      file.path(outdir, "Scatter_nCount_vs_nFeature_bySample_bin2d.pdf"),
      p_scatter,
      width = 30,
      height = 18,
      units = "cm",
      device = cairo_pdf
    )
  }

  meta <- scobj@meta.data
  has_df <- "DoubletFinder_class" %in% colnames(meta)
  summary_table <- meta %>%
    mutate(
      is_doublet = if (has_df) ifelse(tolower(DoubletFinder_class) == "doublet", 1L, 0L) else NA_integer_
    ) %>%
    group_by(orig.ident) %>%
    summarise(
      n_cells = dplyr::n(),
      median_nCount = median(nCount_RNA, na.rm = TRUE),
      median_nFeature = median(nFeature_RNA, na.rm = TRUE),
      median_mito = if ("percent_mito" %in% colnames(.)) median(percent_mito, na.rm = TRUE) else NA_real_,
      median_ribo = if ("percent_ribo" %in% colnames(.)) median(percent_ribo, na.rm = TRUE) else NA_real_,
      median_hb = if ("percent_hb" %in% colnames(.)) median(percent_hb, na.rm = TRUE) else NA_real_,
      doublet_rate = if (has_df) mean(is_doublet, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(desc(n_cells))

  write.csv(summary_table, file.path(outdir, "QC_summary_by_sample.csv"), row.names = FALSE)

  p_bar <- ggplot(summary_table, aes(x = reorder(orig.ident, -n_cells), y = n_cells)) +
    geom_col() +
    theme_bw() +
    xlab("Sample") +
    ylab("Number of nuclei") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

  ggsave(
    file.path(outdir, "Bar_nCells_bySample.pdf"),
    p_bar,
    width = 24,
    height = 10,
    units = "cm",
    device = cairo_pdf
  )

  invisible(summary_table)
}

# =============================================================================
# 2. Load input files
# =============================================================================

if (!file.exists(counts_list_rds)) {
  stop("Count matrix RDS not found: ", counts_list_rds)
}
counts_list <- readRDS(counts_list_rds)
if (!is.list(counts_list) || length(counts_list) < 1) {
  stop("counts_list_rds must contain a non-empty named list of dgCMatrix objects.")
}
if (is.null(names(counts_list)) || any(names(counts_list) == "")) {
  names(counts_list) <- paste0("sample", seq_along(counts_list))
}
sample_names <- names(counts_list)
cat("[INFO] Samples detected: ", paste(sample_names, collapse = ", "), "\n", sep = "")

anno_list <- list()
if (nzchar(anno_rds) && file.exists(anno_rds)) {
  anno_list <- readRDS(anno_rds)
  if (!is.list(anno_list)) {
    warning("Annotation RDS is not a list. It will be ignored.")
    anno_list <- list()
  }
  cat("[INFO] Loaded annotation fields: ", paste(names(anno_list), collapse = ", "), "\n", sep = "")
} else {
  cat("[INFO] No external annotation list was loaded.\n")
}

# =============================================================================
# 3. Per-sample QC
# =============================================================================

cat("\n[MAIN] Starting human per-sample QC.\n")
qc_results <- vector("list", length(sample_names))
names(qc_results) <- sample_names

for (sid in sample_names) {
  counts <- counts_list[[sid]]
  if (!inherits(counts, "dgCMatrix")) {
    stop("Counts for sample ", sid, " are not a dgCMatrix object.")
  }
  qc_results[[sid]] <- qc_one_sample(sample_id = sid, counts = counts, annos = anno_list)
  gc()
}

qc_summary_all <- dplyr::bind_rows(lapply(qc_results, `[[`, "qc_summary"))
df_summary_all <- dplyr::bind_rows(lapply(qc_results, `[[`, "doubletfinder_summary"))

write.csv(qc_summary_all, file.path(qc_root, "human_qc_summary_all_samples.csv"), row.names = FALSE)
write.csv(df_summary_all, file.path(qc_root, "human_doubletfinder_summary_all_samples.csv"), row.names = FALSE)

cat("[MAIN] Human per-sample QC finished.\n")

# =============================================================================
# 4. Merge QC-filtered samples
# =============================================================================

all_merged <- merge_samples(sample_names, "human_all_samples")
if (!"IL6" %in% rownames(all_merged)) {
  stop("IL6 was lost after merging.")
}

# =============================================================================
# 5. Attach sample-level metadata and save clean merged object
# =============================================================================

scobj <- attach_sample_metadata(all_merged, sample_manifest_path)

out_clean <- file.path(reproduce_dir, "human_all_samples_merged.clean.rds")
saveRDS(scobj, out_clean)
cat("[SAVE] ", out_clean, " | dims: ", paste(dim(scobj), collapse = " x "), "\n", sep = "")

# =============================================================================
# 6. Post-QC visualization bundle and summary tables
# =============================================================================

post_qc_summary <- make_post_qc_visualizations(scobj, qc_plot_dir)

# =============================================================================
# 7. Session information
# =============================================================================

write_session_info(file.path(log_dir, "sessionInfo_step1_qc_human.txt"))

cat("\n============================================================\n")
cat("Human Step 1 QC finished successfully.\n")
cat("Reproduce dir: ", reproduce_dir, "\n", sep = "")
cat("QC root      : ", qc_root, "\n", sep = "")
cat("Clean RDS    : ", out_clean, "\n", sep = "")
cat("QC plots dir : ", qc_plot_dir, "\n", sep = "")
cat("============================================================\n")

cat("END!\n")
