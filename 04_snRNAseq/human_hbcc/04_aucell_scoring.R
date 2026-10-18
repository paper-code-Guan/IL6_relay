# =========================================================
# HBCC formal analysis (AUCell-based)
# IL6 single-gene score, IL-6-associated signaling module score,
# GFAP score (Ast-only), and AIF1 score (Mic-only).
#
# Workflow per cell type:
#   1) subset Seurat object by curated cell type
#   2) extract RNA raw count matrix
#   3) build within-cell-type rankings via AUCell::AUCell_buildRankings
#   4) compute AUC scores via AUCell::AUCell_calcAUC at two top windows:
#        - default window (aucMaxRank = top_percent_use * G)        : multi-gene module
#        - extended window (aucMaxRank = top_percent_use_single * G) : single-gene cytokine / marker
#   5) determine positivity:
#        - single-gene (IL6, GFAP, AIF1) : (AUC > 0 at extended window) OR (raw count > 0)
#        - multi-gene module             : AUC > AUCell-selected threshold (default window)
#
# Important interpretation:
#   - Single-gene analyses (IL6, GFAP, AIF1) use raw-count detection protection
#     (rank-and-detection) under sparse-signature transcript dropout (Hicks 2018).
#   - The AUCell threshold for the multi-gene module is selected by the package's
#     default logic (smallestPopPercent = 0.25; candidates include minimumDens,
#     L_k2, R_k3, Global_k1). When AUCell selects the heuristic `tenPercentOfMax`,
#     it is automatically replaced by the data-driven `Global_k1` (mean + 2*SD)
#     per Aibar 2017 vignette guidance for sparse-signature distributions.
#
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(SummarizedExperiment)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
set.seed(123)

# Lock BiocParallel global backend to SerialParam so that any AUCell call
# that does not explicitly pass BPPARAM (AUCell_calcAUC,
# AUCell_exploreThresholds, internal helpers) cannot fork worker processes
# off the registered MulticoreParam default. On HPC nodes the default
# backend is typically MulticoreParam(workers = ~46, cluster type FORK),
# which would copy the per-chunk dense ranking matrix into every worker
# and inflate the resident memory by an order of magnitude. SerialParam
# keeps everything in a single process; chunked AUCell already bounds
# peak memory to one chunk's ranking matrix at a time.
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
message("Active BiocParallel backend: ", class(BiocParallel::bpparam())[1],
        " (workers = ", BiocParallel::bpworkers(BiocParallel::bpparam()), ")")

# =========================================================
# 0. Package checks
# =========================================================
required_aucell_version <- "1.30.0"
if (utils::compareVersion(as.character(packageVersion("AUCell")), required_aucell_version) < 0) {
  stop("AUCell >= ", required_aucell_version, " is required; installed version is ", packageVersion("AUCell"))
}

message("AUCell version: ", as.character(packageVersion("AUCell")))
message("Default argument checks: PASS\n")

# Cell-chunking budget for AUCell ranking on large cell types. Each cell's
# gene rank is computed from that cell's count vector alone (intra-cell
# operation), so chunking by cells does not change rank values, AUC values,
# or the multi-gene-module AUC denominator (which depends only on the
# gene_set_list passed in, identical across chunks). Chunking lowers peak
# memory by bounding the ranking matrix size per AUCell::AUCell_buildRankings
# call. Default 30000 cells per chunk -> peak per chunk approximately
# (30000 cells * G genes * 4 bytes int32) + AUCell working memory; with
# G ~ 30000 this is ~3.6 GB ranking + ~3-4 GB working memory ~= 7-8 GB,
# comfortable on 32 GB systems. Override with `AUCELL_CHUNK_NCELLS` env var
# to scale up (e.g. 80000 on a >=64 GB workstation) or down for tighter
# memory. Cell types with fewer cells than the chunk size are processed in
# a single AUCell call (no chunking overhead).
chunk_ncells_default <- 30000L
chunk_ncells_min     <- 1000L  # defensive floor: AUCell edge-case stability on tiny chunks
.chunk_ncells_env <- suppressWarnings(as.integer(Sys.getenv("AUCELL_CHUNK_NCELLS", unset = "")))
chunk_ncells_use <- if (length(.chunk_ncells_env) == 1L && !is.na(.chunk_ncells_env) && .chunk_ncells_env > 0L) {
  max(.chunk_ncells_env, chunk_ncells_min)
} else {
  chunk_ncells_default
}
message("AUCell chunk_ncells (cell-chunking budget): ", chunk_ncells_use,
        if (length(.chunk_ncells_env) == 1L && !is.na(.chunk_ncells_env) && .chunk_ncells_env > 0L && .chunk_ncells_env < chunk_ncells_min) {
          sprintf(" (env value %d clamped up to floor %d)", .chunk_ncells_env, chunk_ncells_min)
        } else {
          ""
        })

# =========================================================
# 1. Paths and global parameters
# =========================================================
hbcc_rds_env <- Sys.getenv("SEURAT_RDS_PATH", unset = "")
hbcc_base_env <- Sys.getenv("HBCC_BASE_PATH", unset = "")
hbcc_result_root_env <- Sys.getenv("HBCC_RESULT_ROOT", unset = "")

if (!nzchar(hbcc_rds_env)) {
  stop("SEURAT_RDS_PATH must be set for HBCC step4.", call. = FALSE)
}
if (!nzchar(hbcc_base_env) && !nzchar(hbcc_result_root_env)) {
  stop("Set HBCC_RESULT_ROOT or HBCC_BASE_PATH for HBCC step 4.", call. = FALSE)
}

hbcc_rds <- hbcc_rds_env
hbcc_base <- hbcc_base_env

timestamp_tag2 <- format(Sys.time(), "%Y%m%d%H%M")

out_root <- if (nzchar(hbcc_result_root_env)) {
  normalizePath(hbcc_result_root_env, winslash = "/", mustWork = FALSE)
} else {
  file.path(hbcc_base, paste0("results/hbcc_IL6_", timestamp_tag2))
}

step4_dir <- file.path(out_root, "step4_formal")
tab_dir <- file.path(step4_dir, "tables")
score_tab_dir <- file.path(tab_dir, "scores")
summary_tab_dir <- file.path(tab_dir, "summaries")
threshold_tab_dir <- file.path(tab_dir, "thresholds")
mapping_tab_dir <- file.path(tab_dir, "gene_mapping")
stats_cell_tab_dir <- file.path(tab_dir, "stats_celllevel")
qc_tab_dir <- file.path(tab_dir, "qc")
run_tab_dir <- file.path(tab_dir, "run_level")
rds_dir <- file.path(step4_dir, "rds")
log_dir <- file.path(step4_dir, "logs")

for (dir_i in c(out_root, step4_dir, tab_dir, score_tab_dir, summary_tab_dir,
                threshold_tab_dir, mapping_tab_dir, stats_cell_tab_dir,
                qc_tab_dir, run_tab_dir, rds_dir, log_dir)) {
  dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)
}

writeLines(c(
  "HBCC IL6 run output",
  paste0("Run timestamp: ", timestamp_tag2),
  "",
  "Directory structure:",
  "  step4_formal/            Formal AUCell scoring outputs.",
  "  step4_formal/tables/     Classified scoring, summary, threshold, mapping, statistics and QC tables.",
  "  step4_formal/rds/        Per-cell-type and run-level RDS objects.",
  "  step4_formal/logs/       Session information.",
  "  step5_final_<timestamp>/ Final figures and post-processing tables generated by step5.",
  "",
  "File naming uses content descriptions rather than manuscript figure numbering."
), con = file.path(out_root, "README.txt"))

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
        stop("No shared features across ", layer, " layers in assay ", assay, call. = FALSE)
      }
      layer_mats <- lapply(layer_mats, function(mat) mat[common_features, , drop = FALSE])
      merged_layer <- do.call(cbind, layer_mats)
      if (anyDuplicated(colnames(merged_layer))) {
        stop("Duplicated cell names detected after combining ", layer, " layers in assay ", assay, call. = FALSE)
      }
      return(merged_layer)
    }
  }

  if (inherits(assay_obj, "Assay")) {
    legacy_out <- tryCatch(
      Seurat::GetAssayData(object = object, assay = assay, slot = layer),
      error = function(e) e
    )
    if (!inherits(legacy_out, "error")) return(legacy_out)
  }

  stop(conditionMessage(out), call. = FALSE)
}

# =========================================================
# 2. Load Seurat object and harmonize metadata
# =========================================================
stopifnot("HBCC Seurat RDS not found." = file.exists(hbcc_rds))
obj <- readRDS(hbcc_rds)
stopifnot("Loaded object is not a Seurat object." = inherits(obj, "Seurat"))

DefaultAssay(obj) <- "RNA"

# ---------------------------------------------------------------------------
# HBCC-style metadata column compatibility (5-CT-per-RDS workflow).
#
# When the input RDS uses HBCC published-dataset metadata column names
# (Celltype_use / ct4 / Celltype_label, Group / group, donor_id) instead of
# the canonical step4-expected names (Celltypes, Diagnosis, Patient), copy
# the values across so the rest of step4 sees the canonical names. This is
# safe to apply unconditionally because:
#   - The 'if (!"X" %in% colnames(md))' guard short-circuits on legacy
#     RDSes that already have the canonical names.
#   - We never overwrite an existing canonical column.
# ---------------------------------------------------------------------------
{
  .md_pre <- obj@meta.data
  if (!"Celltypes" %in% colnames(.md_pre)) {
    if ("Celltype_use" %in% colnames(.md_pre)) {
      obj$Celltypes <- as.character(.md_pre$Celltype_use)
    } else if ("ct4" %in% colnames(.md_pre)) {
      obj$Celltypes <- as.character(.md_pre$ct4)
    } else if ("Celltype_label" %in% colnames(.md_pre)) {
      obj$Celltypes <- as.character(.md_pre$Celltype_label)
    }
  }
  if (!"Diagnosis" %in% colnames(.md_pre)) {
    if ("Diagnosis" %in% colnames(.md_pre)) {
      # already canonical (no-op)
    } else if ("Group" %in% colnames(.md_pre)) {
      obj$Diagnosis <- as.character(.md_pre$Group)
    } else if ("group" %in% colnames(.md_pre)) {
      obj$Diagnosis <- as.character(.md_pre$group)
    }
  }
  # HBCC: donor_id is the canonical donor ID. Always overwrite Patient
  # from donor_id when available, because HBCC's per-CT RDS frequently
  # carries an uninformative literal Patient column (e.g. "Donor")
  # inherited from upstream pipelines, while the true donor identifier
  # lives in donor_id.
  if ("donor_id" %in% colnames(.md_pre)) {
    obj$Patient <- as.character(.md_pre$donor_id)
  } else if ("Patient" %in% colnames(.md_pre)) {
    obj$Patient <- as.character(.md_pre$Patient)
  } else {
    obj$Patient <- NA_character_
  }
  # Defensive: if Patient is degenerate (single literal value, e.g. all
  # cells == "Donor"), reconstruct from cell-barcode prefix.
  # HBCC convention: '<DonorID>-<lane>-<barcode>-<rep>'.
  if (length(unique(stats::na.omit(obj$Patient))) <= 1L) {
    bc_donor <- sub("-.*", "", colnames(obj))
    if (length(unique(bc_donor)) > 1L) {
      message(sprintf(
        "[Patient] Degenerate Patient column detected (n_unique=%d); reconstructed from cell-barcode prefix (n_unique=%d).",
        length(unique(stats::na.omit(obj$Patient))), length(unique(bc_donor))
      ))
      obj$Patient <- bc_donor
    }
  }
  rm(.md_pre)
}

md <- obj@meta.data

stopifnot("'Diagnosis' column is missing." = "Diagnosis" %in% colnames(md))
stopifnot("'Celltypes' column is missing." = "Celltypes" %in% colnames(md))

display_levels <- c("CTR", "SCZ")
group_colors   <- c("CTR" = "#7FB2D5", "SCZ" = "#F47F72")

grp <- trimws(as.character(md$Diagnosis))
grp[grp %in% c("Ctr", "CTR", "Ctrl", "CTRL", "Control", "CONTROL", "control",
                "Normal", "NORMAL", "normal", "Healthy", "HEALTHY", "healthy")] <- "CTR"
grp[grp %in% c("Scz", "SCZ", "Sz", "SZ", "sz", "Case", "CASE", "case",
                "Schizophrenia", "schizophrenia", "SCHIZOPHRENIA")] <- "SCZ"
# Anything not yet canonicalized is treated as missing for the strict check.
grp[!grp %in% c("CTR", "SCZ")] <- NA_character_
stopifnot("Unresolved Diagnosis values detected." = !any(is.na(grp)))
obj$Group <- factor(grp, levels = display_levels)

# =========================================================
# 3. Cell-type harmonization
# =========================================================
canonicalize_hbcc_ct <- function(v) {
  x  <- trimws(as.character(v))
  xl <- tolower(x)
  out <- rep(NA_character_, length(x))

  out[xl %in% c("ast", "astro", "astrocyte", "astrocytes")] <- "Ast"
  out[xl %in% c("excitneu", "exn", "glutamatergic neuron", "glutamatergic neurons",
                "excitatory neuron", "excitatory neurons")] <- "ExcitNeu"
  out[xl %in% c("inhibneu", "inn", "gabaergic neuron", "gabaergic neurons",
                "inhibitory neuron", "inhibitory neurons", "interneuron", "interneurons")] <- "InhibNeu"
  out[xl %in% c("mic", "micro", "microglia", "microglial")] <- "Mic"
  out[xl %in% c("ol", "oligo", "oligodendrocyte", "oligodendrocytes")] <- "OL"
  out[xl %in% c("other", "others", "other_mix", "other cells", "other cell", "other_cell", "othercell")] <- "Others"

  already_standard <- x %in% c("Ast", "ExcitNeu", "InhibNeu", "Mic", "OL", "Others")
  out[is.na(out) & already_standard] <- x[is.na(out) & already_standard]

  out
}

ct_label_map <- c(
  "Ast"      = "Astrocytes",
  "ExcitNeu" = "Excitatory neurons",
  "InhibNeu" = "Inhibitory neurons",
  "Mic"      = "Microglia",
  "OL"       = "Oligodendrocytes",
  "Others"   = "Others"
)

ct_order <- c("Ast", "ExcitNeu", "InhibNeu", "Mic", "OL", "Others")
# Module-score eligible candidate cell types. Oligo (OL) is included so the
# MSigDB Hallmark IL6/JAK/STAT3 multi-gene module is computed for OL whenever
# the dataset contains it. The actual main_module_cts (computed below) keeps
# only those candidates that are actually present in the data.
main4_candidate_cts <- c("Ast", "ExcitNeu", "InhibNeu", "Mic", "OL")

# Locked per-cell-type nucleus counts of the five HBCC input objects, before the
# donor-level age and comorbidity filter that step 5 applies. Quality-control
# snapshot only; no downstream computation depends on it.
expected_counts_locked <- data.frame(
  Celltype_use = c("Ast", "ExcitNeu", "InhibNeu", "Mic", "OL"),
  expected_n = c(111006, 220466, 179831, 50212, 282705),
  stringsAsFactors = FALSE
)
expected_count_tolerance <- 2L

obj$Celltypes_raw <- as.character(md$Celltypes)
obj$Celltype_use <- canonicalize_hbcc_ct(md$Celltypes)
obj$Celltype_use <- factor(obj$Celltype_use, levels = ct_order)

raw_mapping_check <- data.frame(
  Celltypes_raw = obj$Celltypes_raw,
  Celltype_use  = as.character(obj$Celltype_use),
  stringsAsFactors = FALSE
) %>%
  group_by(Celltypes_raw, Celltype_use) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  arrange(is.na(Celltype_use), desc(n), Celltypes_raw)

write.csv(
  raw_mapping_check,
  file.path(qc_tab_dir, sprintf("hbcc_raw_celltypes_to_celltype_use_mapping_%s.csv", timestamp_tag2)),
  row.names = FALSE
)

message("========== Raw Celltypes -> Celltype_use mapping ==========")
print(raw_mapping_check)
message("")

keep_cells <- !is.na(obj$Celltype_use)
obj <- obj[, keep_cells]
obj$Celltype_use <- droplevels(obj$Celltype_use)
cts_present <- ct_order[ct_order %in% as.character(unique(obj$Celltype_use))]
all_ct_tag <- sprintf("ALL%dCT", length(cts_present))
main_module_cts <- main4_candidate_cts[main4_candidate_cts %in% as.character(unique(obj$Celltype_use))]

celltype_count_table <- data.frame(
  Celltype_use = names(table(obj$Celltype_use)),
  n_cells = as.integer(table(obj$Celltype_use)),
  stringsAsFactors = FALSE
) %>%
  right_join(data.frame(Celltype_use = ct_order, stringsAsFactors = FALSE), by = "Celltype_use") %>%
  left_join(expected_counts_locked, by = "Celltype_use") %>%
  mutate(
    n_cells = ifelse(is.na(n_cells), 0L, n_cells),
    Celltype_label = unname(ct_label_map[Celltype_use]),
    abs_diff = ifelse(is.na(expected_n), NA_real_, abs(n_cells - expected_n)),
    rel_diff = ifelse(is.na(expected_n), NA_real_, ifelse(expected_n > 0, abs_diff / expected_n, NA_real_)),
    count_tolerance = expected_count_tolerance,
    count_matches_expected = ifelse(is.na(expected_n), NA, n_cells == expected_n),
    count_within_tolerance = ifelse(is.na(expected_n), TRUE, abs_diff <= count_tolerance)
  ) %>%
  arrange(match(Celltype_use, ct_order))

message("========== HBCC celltype count table ==========")
print(celltype_count_table)
message("")

write.csv(
  celltype_count_table,
  file.path(qc_tab_dir, sprintf("hbcc_celltype_counts_expected_check_%s.csv", timestamp_tag2)),
  row.names = FALSE
)

if (!all(celltype_count_table$count_within_tolerance, na.rm = TRUE)) {
  warning("Some cell-type counts differ from the locked analysis snapshot beyond tolerance. Please inspect the exported count table.")
} else if (!all(celltype_count_table$count_matches_expected, na.rm = TRUE)) {
  message("Cell-type counts show only minor deviations within the locked-snapshot tolerance.")
}

# =========================================================
# 4. Helper functions
# =========================================================
resolve_features_ci <- function(symbols, pool_rows) {
  idx <- match(tolower(symbols), tolower(pool_rows))
  resolved <- ifelse(is.na(idx), NA_character_, pool_rows[idx])
  list(
    found = unique(na.omit(resolved)),
    missing = unique(symbols[is.na(idx)]),
    mapping = data.frame(
      requested = symbols,
      resolved = unname(resolved),
      stringsAsFactors = FALSE
    )
  )
}

safe_file_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

extract_stats_table <- function(x, ct_name, feature_name, level_name) {
  if (is.null(x)) return(NULL)
  out <- x
  out$Celltype_use   <- ct_name
  out$Celltype_label <- unname(ct_label_map[ct_name])
  out$feature_set    <- feature_name
  out$stats_level    <- level_name
  out
}

# Compute cell-level stats: per-group mean / median, log2 fold change,
# Wilcoxon rank-sum p, and per-group AUCell-positive cell counts.
compute_celllevel_stats <- function(score_vec, group_vec, positive_mask,
                                    feature_label, group_levels) {
  groups <- group_levels
  group1 <- groups[1]
  group2 <- groups[2]

  s1 <- score_vec[group_vec == group1]
  s2 <- score_vec[group_vec == group2]
  s1 <- s1[is.finite(s1)]
  s2 <- s2[is.finite(s2)]

  m1 <- if (length(s1)) mean(s1) else NA_real_
  m2 <- if (length(s2)) mean(s2) else NA_real_
  med1 <- if (length(s1)) stats::median(s1) else NA_real_
  med2 <- if (length(s2)) stats::median(s2) else NA_real_

  p_value <- if (length(s1) > 0 && length(s2) > 0) {
    tryCatch(suppressWarnings(stats::wilcox.test(s1, s2)$p.value),
             error = function(e) NA_real_)
  } else NA_real_

  eps <- 1e-5
  log2fc <- log2((m1 + eps) / (m2 + eps))

  pm <- as.logical(positive_mask)
  n_pos1 <- sum(pm & group_vec == group1, na.rm = TRUE)
  n_pos2 <- sum(pm & group_vec == group2, na.rm = TRUE)

  setNames(
    data.frame(
      feature_label, log2fc, m1, m2, med1, med2, p_value, n_pos1, n_pos2,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    c(
      "feature",
      paste0("log2FC_", group1, "vs", group2, "_cell"),
      paste0("mean_", group1, "_cell"),
      paste0("mean_", group2, "_cell"),
      paste0("median_", group1, "_cell"),
      paste0("median_", group2, "_cell"),
      "p_cell",
      paste0("n_aucell_pos_", group1),
      paste0("n_aucell_pos_", group2)
    )
  )
}

# Per-group summary of total cells and positive cells across multiple
# positivity criteria. Reports four mask-derived counts so the audit
# table records each criterion independently:
#   - aucell_mask   : AUC at default top window > AUCell-selected threshold
#   - rank_mask     : AUC at the active top window > 0
#                     (extended window for single-gene; equals aucell for multi-gene)
#   - detection_mask: raw transcript count > 0 (single-gene only; FALSE elsewhere)
#   - positive_mask : the binary positivity actually used downstream
summarize_aucell_positives <- function(aucell_mask, rank_mask, detection_mask, positive_mask,
                                       group_vec, group_levels,
                                       threshold_selected, threshold_method,
                                       positive_rule) {
  am <- as.logical(aucell_mask)
  rm <- as.logical(rank_mask)
  dm <- as.logical(detection_mask)
  pm <- as.logical(positive_mask)
  rows <- lapply(group_levels, function(g) {
    in_group <- group_vec == g
    total <- sum(in_group)
    pos_combined <- sum(pm & in_group, na.rm = TRUE)
    pos_aucell   <- sum(am & in_group, na.rm = TRUE)
    pos_rank     <- sum(rm & in_group, na.rm = TRUE)
    pos_detect   <- sum(dm & in_group, na.rm = TRUE)
    data.frame(
      Group = g,
      total_cells = total,
      score_pos = pos_combined,
      score_prop = if (total > 0) pos_combined / total else NA_real_,
      score_aucell_pos = pos_aucell,
      score_aucell_prop = if (total > 0) pos_aucell / total else NA_real_,
      score_rank_pos = pos_rank,
      score_rank_prop = if (total > 0) pos_rank / total else NA_real_,
      score_detection_pos = pos_detect,
      score_detection_prop = if (total > 0) pos_detect / total else NA_real_,
      score_combined_pos = pos_combined,
      score_combined_prop = if (total > 0) pos_combined / total else NA_real_,
      threshold_selected = threshold_selected,
      threshold_method = threshold_method,
      positive_rule = positive_rule,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Robust extraction of the AUCell-selected threshold and method name from
# AUCell::AUCell_exploreThresholds output. The function tolerates several
# possible internal structures across AUCell minor versions.
extract_aucell_threshold_info <- function(explore_res, gs_name, auc_vector) {
  fallback_global_k1 <- function() {
    s <- auc_vector[is.finite(auc_vector)]
    if (length(s) > 1 && stats::sd(s) > 0) {
      mean(s) + 2 * stats::sd(s)
    } else {
      max(s, 0, na.rm = TRUE) + .Machine$double.eps
    }
  }

  if (is.null(explore_res) || is.null(explore_res[[gs_name]])) {
    thr <- fallback_global_k1()
    return(list(
      selected = thr,
      method = "Global_k1_fallback",
      candidates = data.frame(
        method = "Global_k1_fallback",
        threshold = thr,
        nCells = sum(auc_vector > thr, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    ))
  }

  res <- explore_res[[gs_name]]
  aucThr <- res$aucThr

  selected_value <- NA_real_
  if (!is.null(aucThr$selected)) {
    sv <- suppressWarnings(as.numeric(aucThr$selected))
    if (length(sv) == 1L && is.finite(sv)) selected_value <- sv
  }
  if (!is.finite(selected_value) && !is.null(aucThr$threshold)) {
    sv <- suppressWarnings(as.numeric(aucThr$threshold))
    if (length(sv) == 1L && is.finite(sv)) selected_value <- sv
  }
  if (!is.finite(selected_value)) {
    selected_value <- fallback_global_k1()
  }

  candidates <- NULL
  method_name <- NA_character_
  thr_df <- aucThr$thresholds
  if (!is.null(thr_df)) {
    thr_df_df <- as.data.frame(thr_df, stringsAsFactors = FALSE)
    method_rn <- rownames(thr_df_df)
    thr_col <- if ("threshold" %in% colnames(thr_df_df)) "threshold" else colnames(thr_df_df)[1]
    thr_vals <- suppressWarnings(as.numeric(thr_df_df[[thr_col]]))
    nc_col <- if ("nCells_AUC_above_thr" %in% colnames(thr_df_df)) {
      "nCells_AUC_above_thr"
    } else if (ncol(thr_df_df) >= 2) {
      colnames(thr_df_df)[2]
    } else {
      NA_character_
    }
    n_cells_vec <- if (!is.na(nc_col)) suppressWarnings(as.numeric(thr_df_df[[nc_col]])) else rep(NA_real_, length(thr_vals))
    candidates <- data.frame(
      method = method_rn,
      threshold = thr_vals,
      nCells = n_cells_vec,
      stringsAsFactors = FALSE
    )
    if (length(thr_vals) && any(is.finite(thr_vals))) {
      idx <- which.min(abs(thr_vals - selected_value))
      method_name <- method_rn[idx]
    }
  }
  if (is.na(method_name)) method_name <- "AUCell_default"
  if (is.null(candidates)) {
    candidates <- data.frame(
      method = method_name,
      threshold = selected_value,
      nCells = sum(auc_vector > selected_value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  # When AUCell selects the non-data-driven `tenPercentOfMax` heuristic
  # (10% of maximum AUC), prefer the data-driven `Global_k1` candidate
  # (mean + 2*SD) which AUCell guidance recommends for sparse-signature
  # distributions (Aibar et al. 2017, Nat Methods). Falls through if
  # Global_k1 is unavailable.
  if (identical(method_name, "tenPercentOfMax")) {
    globk1_idx <- which(candidates$method == "Global_k1")
    if (length(globk1_idx)) {
      gk1_thr <- as.numeric(candidates$threshold[globk1_idx[1]])
      if (is.finite(gk1_thr)) {
        selected_value <- gk1_thr
        method_name <- "Global_k1"
      }
    }
  }

  list(
    selected = selected_value,
    method = method_name,
    candidates = candidates
  )
}

# Compute AUCell rankings + AUC matrix + threshold exploration for a
# named gene-set list, on a counts matrix subset to a single cell type.
# Returns AUC at two `aucMaxRank` settings: a default top window
# (`top_percent_default`, used for multi-gene module signatures and
# AUCell threshold exploration) and an extended top window
# (`top_percent_single`, used for single-gene cytokine analysis where
# extreme transcript dropout requires a wider rank window per AUCell
# vignette guidance for sparse-signature regimes).
#
# Memory-bounded mode: when the cell type has more cells than
# `chunk_ncells`, AUCell ranking is run in cell-chunked batches and the
# resulting per-chunk AUC matrices are concatenated column-wise. Per-cell
# rank is an intra-cell operation (depends only on that cell's count
# vector), so cell-chunking is mathematically equivalent to a single
# unchunked AUCell call: rank values, AUC values, and the multi-gene
# module-score AUC denominator (a function only of `gene_set_list` and
# `aucMaxRank`) are unchanged. Threshold exploration on the merged AUC
# matrix is performed via an `aucellResults` wrapper around the merged
# matrix; if the wrap is not supported by the installed AUCell version,
# `extract_aucell_threshold_info` falls back to the data-driven
# Global_k1 (mean + 2*SD) per gene set.
compute_aucell_for_celltype <- function(counts_mat, gene_set_list,
                                        top_percent_default = 0.05,
                                        top_percent_single  = 0.10,
                                        chunk_ncells        = 30000L) {
  if (!inherits(counts_mat, "dgCMatrix")) {
    counts_mat <- as(counts_mat, "CsparseMatrix")
  }

  n_genes_total <- nrow(counts_mat)
  n_cells_total <- ncol(counts_mat)

  K_top_default <- max(1L, ceiling(top_percent_default * n_genes_total))
  K_top_single  <- max(1L, ceiling(top_percent_single  * n_genes_total))

  # chunk_ncells floor is enforced by the caller (see chunk_ncells_use
  # at the top of the script); this branch only checks whether the
  # current cell type is large enough to need chunking.
  use_chunked <- isTRUE(n_cells_total > chunk_ncells)
  gene_set_names <- names(gene_set_list)

  build_one_block <- function(counts_block) {
    AUCell::AUCell_buildRankings(
      counts_block,
      plotStats = FALSE,
      verbose   = FALSE,
      BPPARAM   = BiocParallel::SerialParam()
    )
  }

  if (!use_chunked) {
    rankings <- build_one_block(counts_mat)
    auc_res_default <- AUCell::AUCell_calcAUC(
      geneSets   = gene_set_list,
      rankings   = rankings,
      aucMaxRank = K_top_default,
      verbose    = FALSE,
      nCores     = 1L
    )
    auc_res_single <- AUCell::AUCell_calcAUC(
      geneSets   = gene_set_list,
      rankings   = rankings,
      aucMaxRank = K_top_single,
      verbose    = FALSE,
      nCores     = 1L
    )
    auc_mat_default <- AUCell::getAUC(auc_res_default)
    auc_mat_single  <- AUCell::getAUC(auc_res_single)

    explore_res_default <- tryCatch(
      AUCell::AUCell_exploreThresholds(
        auc_res_default,
        plotHist = FALSE,
        smallestPopPercent = 0.25,
        assignCells = TRUE,
        verbose = FALSE,
        nCores = 1L
      ),
      error = function(e) {
        message("  AUCell_exploreThresholds failed (default window); per-gene set fallback to Global_k1: ",
                conditionMessage(e))
        NULL
      }
    )

    rm(rankings, auc_res_default, auc_res_single)
    invisible(gc(verbose = FALSE))
  } else {
    # Balanced split: produces n_chunks chunks of ~equal size (last chunk
    # is at most balanced_size; avoids tiny final chunks that could trip
    # AUCell edge cases).
    n_chunks       <- as.integer(ceiling(n_cells_total / chunk_ncells))
    balanced_size  <- as.integer(ceiling(n_cells_total / n_chunks))
    cell_chunks    <- split(
      seq_len(n_cells_total),
      ceiling(seq_len(n_cells_total) / balanced_size)
    )
    n_chunks_actual <- length(cell_chunks)
    message(sprintf(
      "  Chunked AUCell: %d cells -> %d chunks of <=%d cells (chunk_ncells=%d).",
      n_cells_total, n_chunks_actual, balanced_size, chunk_ncells
    ))

    auc_mat_default_chunks <- vector("list", n_chunks_actual)
    auc_mat_single_chunks  <- vector("list", n_chunks_actual)

    for (i in seq_along(cell_chunks)) {
      chunk_idx <- cell_chunks[[i]]
      message(sprintf("    Chunk %d/%d : %d cells", i, n_chunks_actual, length(chunk_idx)))
      counts_chunk <- counts_mat[, chunk_idx, drop = FALSE]
      rk_chunk <- build_one_block(counts_chunk)
      auc_d <- AUCell::AUCell_calcAUC(
        geneSets   = gene_set_list,
        rankings   = rk_chunk,
        aucMaxRank = K_top_default,
        verbose    = FALSE,
        nCores     = 1L
      )
      auc_s <- AUCell::AUCell_calcAUC(
        geneSets   = gene_set_list,
        rankings   = rk_chunk,
        aucMaxRank = K_top_single,
        verbose    = FALSE,
        nCores     = 1L
      )
      auc_mat_default_chunks[[i]] <- AUCell::getAUC(auc_d)
      auc_mat_single_chunks[[i]]  <- AUCell::getAUC(auc_s)
      rm(rk_chunk, auc_d, auc_s, counts_chunk)
      invisible(gc(verbose = FALSE))
    }

    # Concatenate AUC matrices column-wise. Defensive: enforce a
    # consistent gene-set ordering across chunks (in practice all chunks
    # share rownames since gene_set_list is identical, but we explicitly
    # align so any AUCell minor-version difference cannot shuffle rows).
    align_auc_chunk <- function(m, gs_names) {
      if (identical(rownames(m), gs_names)) return(m)
      out <- matrix(0, nrow = length(gs_names), ncol = ncol(m),
                    dimnames = list(gs_names, colnames(m)))
      present <- intersect(gs_names, rownames(m))
      if (length(present)) {
        out[present, ] <- m[present, , drop = FALSE]
      }
      out
    }
    auc_mat_default <- do.call(
      cbind,
      lapply(auc_mat_default_chunks, align_auc_chunk, gene_set_names)
    )
    auc_mat_single  <- do.call(
      cbind,
      lapply(auc_mat_single_chunks,  align_auc_chunk, gene_set_names)
    )
    rm(auc_mat_default_chunks, auc_mat_single_chunks)
    invisible(gc(verbose = FALSE))

    # Threshold exploration on merged default-window AUC. Wrap the merged
    # matrix in an aucellResults object (assay name "AUC", which is the
    # convention AUCell::AUCell_exploreThresholds dispatches on via
    # getAUC()). On any wrap failure (e.g. AUCell minor-version skew),
    # explore_res_default is set to NULL and the per-gene-set extractor
    # `extract_aucell_threshold_info` falls back to Global_k1
    # (mean + 2*SD). For single-gene cytokines this fallback has no
    # downstream effect because positivity uses rank-and-detection,
    # not the threshold; for the multi-gene module the Global_k1 path is
    # the documented fallback.
    explore_res_default <- tryCatch({
      auc_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(AUC = auc_mat_default)
      )
      aucell_results_fn <- tryCatch(
        get("aucellResults", envir = asNamespace("AUCell")),
        error = function(e) NULL
      )
      if (is.null(aucell_results_fn)) {
        stop("AUCell::aucellResults constructor not found in installed AUCell namespace.")
      }
      auc_obj <- aucell_results_fn(auc_se)
      AUCell::AUCell_exploreThresholds(
        auc_obj,
        plotHist = FALSE,
        smallestPopPercent = 0.25,
        assignCells = TRUE,
        verbose = FALSE,
        nCores = 1L
      )
    }, error = function(e) {
      message("  Chunked AUCell_exploreThresholds wrap failed; per-gene set fallback to Global_k1: ",
              conditionMessage(e))
      NULL
    })
  }

  list(
    auc_matrix_default = auc_mat_default,
    auc_matrix_single  = auc_mat_single,
    explore_res        = explore_res_default,
    K_top_default      = K_top_default,
    K_top_single       = K_top_single,
    n_genes_total      = n_genes_total,
    chunked            = use_chunked,
    chunk_ncells_used  = if (use_chunked) chunk_ncells else NA_integer_
  )
}

# =========================================================
# 5. Analysis settings
# =========================================================
top_percent_use         <- 0.05  # AUCell aucMaxRank fraction (default top window) for multi-gene module signatures
top_percent_use_single  <- 0.10  # extended top window for single-gene cytokine analysis (Aibar 2017 vignette: 1-20% allowed for sparse signatures)

msigdb_il6_jak_stat3_feature <- "MSigDB_Hallmark_IL6_JAK_STAT3_signaling"
msigdb_il6_jak_stat3_genes <- c(
  "A2M", "ACVR1B", "ACVRL1", "BAK1", "CBL", "CCL7", "CCR1", "CD14",
  "CD36", "CD38", "CD44", "CD9", "CNTFR", "CRLF2", "CSF1", "CSF2",
  "CSF2RA", "CSF2RB", "CSF3R", "CXCL1", "CXCL10", "CXCL11", "CXCL13", "CXCL3",
  "CXCL9", "DNTT", "EBI3", "FAS", "GRB2", "HAX1", "HMOX1", "IFNAR1",
  "IFNGR1", "IFNGR2", "IL10RB", "IL12RB1", "IL13RA1", "IL15RA", "IL17RA", "IL17RB",
  "IL18R1", "IL1B", "IL1R1", "IL1R2", "IL2RA", "IL2RG", "IL3RA", "IL4R",
  "IL6", "IL6ST", "IL7", "IL9R", "INHBE", "IRF1", "IRF9", "ITGA4",
  "ITGB3", "JUN", "LEPR", "LTB", "LTBR", "MAP3K8", "MYD88", "OSMR",
  "PDGFC", "PF4", "PIK3R5", "PIM1", "PLA2G2A", "PTPN1", "PTPN11", "PTPN2",
  "REG1A", "SOCS1", "SOCS3", "STAM2", "STAT1", "STAT2", "STAT3", "TGFB1",
  "TLR2", "TNF", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21", "TYK2"
)

feature_sets <- list(
  IL6 = list(
    genes = c("IL6"),
    celltypes = ct_order,
    single_gene = TRUE,
    minimum_gene_fraction = NA_real_,
    positive_rule = "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0"
  ),
  GFAP = list(
    genes = c("GFAP"),
    celltypes = c("Ast"),
    single_gene = TRUE,
    minimum_gene_fraction = NA_real_,
    positive_rule = "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0"
  ),
  AIF1 = list(
    genes = c("AIF1"),
    celltypes = c("Mic"),
    single_gene = TRUE,
    minimum_gene_fraction = NA_real_,
    positive_rule = "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0"
  )
)

feature_sets[[msigdb_il6_jak_stat3_feature]] <- list(
  genes = msigdb_il6_jak_stat3_genes,
  celltypes = main_module_cts,
  single_gene = FALSE,
  minimum_gene_fraction = 0.60,
  module_label = "MSigDB Hallmark IL6/JAK/STAT3 signaling",
  positive_rule = "AUC score > AUCell selected threshold (default top window aucMaxRank=0.05*G)"
)

score_tables <- list()
summary_tables <- list()
threshold_tables <- list()
stats_celllevel_tables <- list()
gene_mapping_tables <- list()
run_log <- list()

# =========================================================
# 6. Core function: one feature set in one cell type
# =========================================================
run_one_feature_set <- function(
  seu_ct, ct_name, feature_name, requested_genes, feature_info,
  meta_ct, aucell_pack
) {
  message("  Feature set: ", feature_name)

  feature_map <- resolve_features_ci(requested_genes, rownames(aucell_pack$auc_matrix_source))
  features_found <- feature_map$found

  if (length(features_found) == 0L) {
    warning("No features found for ", feature_name, " in ", ct_name, ". Skipping.")
    return(NULL)
  }

  expected_n <- length(requested_genes)
  found_n <- length(features_found)
  fraction_recovered <- found_n / expected_n
  missing_genes <- feature_map$missing
  gene_recovery_threshold <- if (!is.null(feature_info$minimum_gene_fraction) &&
                                 is.finite(feature_info$minimum_gene_fraction)) {
    feature_info$minimum_gene_fraction
  } else {
    0.60
  }

  if (feature_name == msigdb_il6_jak_stat3_feature) {
    if (fraction_recovered < gene_recovery_threshold) {
      warning(sprintf(
        "[%s/%s] Only %d/%d (%.1f%%) genes recovered. Skipping this combination. Missing: %s",
        ct_name, feature_name, found_n, expected_n, 100 * fraction_recovered,
        paste(missing_genes, collapse = ";")
      ))
      return(NULL)
    } else if (found_n != expected_n) {
      message(sprintf(
        "[%s/%s] %d/%d genes recovered; analysis proceeds with available subset. Missing: %s",
        ct_name, feature_name, found_n, expected_n, paste(missing_genes, collapse = ";")
      ))
    }
  }

  # AUC vectors at default and single-gene-extended top windows.
  is_single_gene <- length(features_found) == 1L
  use_detection_protection <- is_single_gene

  auc_mat_default <- aucell_pack$auc_matrix_default
  auc_mat_single  <- aucell_pack$auc_matrix_single
  if (!feature_name %in% rownames(auc_mat_default)) {
    warning("Gene set not present in AUC matrix: ", feature_name)
    return(NULL)
  }
  cells_in_order <- colnames(seu_ct)

  # Multi-gene module: report AUC at default window. Single-gene cytokine:
  # report AUC at extended window (the score that defines positivity); the
  # default-window AUC is also retained for transparency in the threshold
  # candidates table.
  score_vec_default <- as.numeric(auc_mat_default[feature_name, cells_in_order])
  names(score_vec_default) <- cells_in_order
  score_vec_single  <- as.numeric(auc_mat_single[feature_name, cells_in_order])
  names(score_vec_single) <- cells_in_order

  score_vec <- if (is_single_gene) score_vec_single else score_vec_default
  K_top_used <- if (is_single_gene) aucell_pack$K_top_single else aucell_pack$K_top_default

  # Threshold via AUCell default (with tenPercentOfMax -> Global_k1 fallback);
  # this threshold is applied to the default-window AUC and gates multi-gene
  # module positivity. For single-gene cytokines it is recorded for
  # transparency but not used as the binary positivity gate.
  thr_info <- extract_aucell_threshold_info(aucell_pack$explore_res, feature_name, score_vec_default)
  threshold_selected <- thr_info$selected
  threshold_method   <- thr_info$method
  candidates_df      <- thr_info$candidates

  aucell_mask <- score_vec_default > threshold_selected
  aucell_mask[!is.finite(score_vec_default)] <- FALSE

  if (use_detection_protection) {
    # Single-gene (IL6, GFAP, AIF1): rank-and-detection positivity at extended
    # top window (aucMaxRank = top_percent_use_single * G).
    #   positive = (AUC > 0 at extended window, equivalent to within-cell rank
    #               <= K_top_single) OR (raw transcript count > 0)
    count_row <- features_found[1]
    counts_src <- aucell_pack$auc_matrix_source
    if (count_row %in% rownames(counts_src)) {
      expr_count_proxy <- as.numeric(counts_src[count_row, cells_in_order, drop = TRUE])
    } else {
      expr_count_proxy <- rep(0, length(cells_in_order))
    }
    names(expr_count_proxy) <- cells_in_order
    detection_mask <- expr_count_proxy > 0
    detection_mask[!is.finite(expr_count_proxy)] <- FALSE

    rank_mask <- score_vec_single > 0
    rank_mask[!is.finite(score_vec_single)] <- FALSE

    positive_mask <- rank_mask | detection_mask
  } else {
    # Multi-gene module: AUCell-threshold positivity (default window).
    expr_count_proxy <- rep(NA_real_, length(cells_in_order))
    names(expr_count_proxy) <- cells_in_order
    detection_mask <- rep(FALSE, length(cells_in_order))
    rank_mask <- aucell_mask
    positive_mask <- aucell_mask
  }

  group_vec <- as.character(meta_ct[cells_in_order, "Group", drop = TRUE])

  positive_rule <- if (use_detection_protection) {
    "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0"
  } else {
    "AUC score > AUCell selected threshold (default top window aucMaxRank=0.05*G)"
  }

  sumtab <- summarize_aucell_positives(
    aucell_mask       = aucell_mask,
    rank_mask         = rank_mask,
    detection_mask    = detection_mask,
    positive_mask     = positive_mask,
    group_vec         = group_vec,
    group_levels      = display_levels,
    threshold_selected = threshold_selected,
    threshold_method  = threshold_method,
    positive_rule     = positive_rule
  )
  sumtab$Celltype_use   <- ct_name
  sumtab$Celltype_label <- unname(ct_label_map[ct_name])
  sumtab$feature_set    <- feature_name

  feature_label <- if (length(features_found) == 1L) features_found[1] else "featureSet"
  stats_cell <- compute_celllevel_stats(
    score_vec      = score_vec,
    group_vec      = group_vec,
    positive_mask  = positive_mask,
    feature_label  = feature_label,
    group_levels   = display_levels
  )
  stats_cell <- extract_stats_table(stats_cell, ct_name, feature_name, "celllevel")

  patient_vec <- if ("Patient" %in% colnames(meta_ct)) {
    as.character(meta_ct[cells_in_order, "Patient", drop = TRUE])
  } else if ("donor_id" %in% colnames(meta_ct)) {
    as.character(meta_ct[cells_in_order, "donor_id", drop = TRUE])
  } else if ("sample" %in% colnames(meta_ct)) {
    as.character(meta_ct[cells_in_order, "sample", drop = TRUE])
  } else if ("orig.ident" %in% colnames(meta_ct)) {
    as.character(meta_ct[cells_in_order, "orig.ident", drop = TRUE])
  } else {
    stop("None of 'Patient' / 'donor_id' / 'sample' / 'orig.ident' available in meta_ct; cannot resolve donor IDs.")
  }
  # Defensive: if patient_vec is degenerate (single literal value, e.g.
  # "Donor"), reconstruct from cell-barcode prefix.
  if (length(unique(stats::na.omit(patient_vec))) <= 1L) {
    bc_donor <- sub("-.*", "", cells_in_order)
    if (length(unique(bc_donor)) > 1L) {
      message(sprintf(
        "[%s/%s] Patient degenerate (n_unique=%d); reconstructed from cell-barcode prefix (n_unique=%d).",
        ct_name, feature_name,
        length(unique(stats::na.omit(patient_vec))), length(unique(bc_donor))
      ))
      patient_vec <- bc_donor
    }
  }

  score_df <- data.frame(
    cell = cells_in_order,
    Group = as.character(meta_ct[cells_in_order, "Group", drop = TRUE]),
    Patient = patient_vec,
    Celltypes_raw = as.character(meta_ct[cells_in_order, "Celltypes_raw", drop = TRUE]),
    Celltype_use = ct_name,
    Celltype_label = unname(ct_label_map[ct_name]),
    feature_set = feature_name,
    features_requested = paste(requested_genes, collapse = ";"),
    features_found = paste(features_found, collapse = ";"),
    n_features_found = length(features_found),
    score = score_vec,
    score_default_window = score_vec_default,
    expr_count_proxy = expr_count_proxy,
    threshold_selected = threshold_selected,
    threshold_method = threshold_method,
    K_top_used = K_top_used,
    aucell_mask = as.logical(aucell_mask),
    detection_mask = as.logical(detection_mask),
    rank_mask = as.logical(rank_mask),
    positive = as.logical(positive_mask),
    detection_used = use_detection_protection,
    stringsAsFactors = FALSE
  )

  threshold_df <- candidates_df
  if (!is.null(threshold_df)) {
    threshold_df$Celltype_use <- ct_name
    threshold_df$Celltype_label <- unname(ct_label_map[ct_name])
    threshold_df$feature_set <- feature_name
  }

  gene_map_df <- feature_map$mapping %>%
    mutate(
      Celltype_use = ct_name,
      Celltype_label = unname(ct_label_map[ct_name]),
      feature_set = feature_name,
      found_in_expr = !is.na(resolved),
      expected_n = expected_n,
      found_n = found_n,
      missing_n = expected_n - found_n,
      missing_genes = paste(missing_genes, collapse = ";"),
      fraction_recovered = fraction_recovered
    )

  ct_id <- safe_file_id(ct_name)
  feature_id <- safe_file_id(feature_name)

  write.csv(score_df, file.path(score_tab_dir, sprintf("%s_%s_scores_%s.csv", ct_id, feature_id, timestamp_tag2)), row.names = FALSE)
  write.csv(sumtab, file.path(summary_tab_dir, sprintf("%s_%s_summary_%s.csv", ct_id, feature_id, timestamp_tag2)), row.names = FALSE)

  if (!is.null(threshold_df)) {
    write.csv(threshold_df, file.path(threshold_tab_dir, sprintf("%s_%s_threshold_candidates_%s.csv", ct_id, feature_id, timestamp_tag2)), row.names = FALSE)
  }

  write.csv(gene_map_df, file.path(mapping_tab_dir, sprintf("%s_%s_gene_mapping_%s.csv", ct_id, feature_id, timestamp_tag2)), row.names = FALSE)

  if (!is.null(stats_cell)) {
    write.csv(stats_cell, file.path(stats_cell_tab_dir, sprintf("%s_%s_celllevel_stats_%s.csv", ct_id, feature_id, timestamp_tag2)), row.names = FALSE)
  }

  list(
    score_table = score_df,
    summary_table = sumtab,
    threshold_table = threshold_df,
    gene_mapping = gene_map_df,
    stats_celllevel = stats_cell,
    threshold_selected = threshold_selected,
    threshold_method = threshold_method,
    expected_n = expected_n,
    found_n = found_n,
    missing_n = expected_n - found_n,
    missing_genes = paste(missing_genes, collapse = ";"),
    fraction_recovered = fraction_recovered
  )
}

# =========================================================
# 7. Core function: one cell type
# =========================================================
run_one_ct <- function(seurat_obj, ct_name) {
  message("------------------------------------------------------------")
  message("Running CT: ", ct_name, " | ", unname(ct_label_map[ct_name]))

  cells_ct <- colnames(seurat_obj)[seurat_obj$Celltype_use == ct_name]

  if (length(cells_ct) == 0L) {
    warning("Skipping empty cell type: ", ct_name)
    return(NULL)
  }

  # Filter feature sets applicable to this cell type (preserves human-specific
  # restriction: GFAP only in Ast; AIF1 only in Mic; multi-gene module only in main four cell types).
  applicable_feature_names <- names(feature_sets)[
    vapply(feature_sets, function(fs) ct_name %in% fs$celltypes, logical(1))
  ]
  if (!length(applicable_feature_names)) {
    message("No feature sets applicable to ", ct_name, "; skipping.")
    return(NULL)
  }

  seu_ct <- subset(seurat_obj, cells = cells_ct)
  meta_ct <- seu_ct@meta.data

  message("n cells = ", ncol(seu_ct))
  print(table(seu_ct$Group))

  counts_ct <- seurat_get_layer(seu_ct, assay = "RNA", layer = "counts")

  message("Building AUCell rankings and AUC matrix...")
  gene_set_list <- lapply(feature_sets[applicable_feature_names], function(fs) {
    fm <- resolve_features_ci(fs$genes, rownames(counts_ct))
    fm$found
  })
  gene_set_list <- gene_set_list[vapply(gene_set_list, length, integer(1)) > 0L]

  if (!length(gene_set_list)) {
    warning("No gene sets resolved for CT: ", ct_name, "; skipping.")
    return(NULL)
  }

  aucell_pack <- compute_aucell_for_celltype(
    counts_mat = counts_ct,
    gene_set_list = gene_set_list,
    top_percent_default = top_percent_use,
    top_percent_single  = top_percent_use_single,
    chunk_ncells = chunk_ncells_use
  )
  # Attach the gene name pool that was used to resolve features (rownames of counts_ct).
  aucell_pack$auc_matrix_source <- counts_ct

  ct_results <- list()

  for (feature_name in applicable_feature_names) {
    ct_results[[feature_name]] <- run_one_feature_set(
      seu_ct = seu_ct,
      ct_name = ct_name,
      feature_name = feature_name,
      requested_genes = feature_sets[[feature_name]]$genes,
      feature_info = feature_sets[[feature_name]],
      meta_ct = meta_ct,
      aucell_pack = aucell_pack
    )
  }

  ct_bundle <- list(
    ct = ct_name,
    ct_label = unname(ct_label_map[ct_name]),
    cells = colnames(seu_ct),
    n_cells = ncol(seu_ct),
    group_table = table(seu_ct$Group),
    aucell_threshold = lapply(ct_results, function(rr) {
      if (is.null(rr)) return(NULL)
      list(
        selected = rr$threshold_selected,
        method = rr$threshold_method,
        candidates = rr$threshold_table
      )
    }),
    results = ct_results
  )

  saveRDS(
    ct_bundle,
    file.path(rds_dir, sprintf("hbcc_%s_%s.rds", safe_file_id(ct_name), timestamp_tag2))
  )

  rm(seu_ct, counts_ct, aucell_pack)
  gc(verbose = FALSE)

  ct_bundle
}

# =========================================================
# 8. Run all target cell types
# =========================================================
all_ct_results <- list()

for (ct in cts_present) {
  res_ct <- run_one_ct(obj, ct)
  all_ct_results[[ct]] <- res_ct

  if (!is.null(res_ct)) {
    for (feature_name in names(res_ct$results)) {
      rr <- res_ct$results[[feature_name]]
      if (is.null(rr)) next

      key <- paste(ct, feature_name, sep = "__")
      score_tables[[key]] <- rr$score_table
      summary_tables[[key]] <- rr$summary_table
      threshold_tables[[key]] <- rr$threshold_table
      gene_mapping_tables[[key]] <- rr$gene_mapping
      stats_celllevel_tables[[key]] <- rr$stats_celllevel
      run_log[[key]] <- data.frame(
        Celltype_use = ct,
        Celltype_label = unname(ct_label_map[ct]),
        feature_set = feature_name,
        n_cells = res_ct$n_cells,
        expected_n = rr$expected_n,
        found_n = rr$found_n,
        missing_n = rr$missing_n,
        missing_genes = rr$missing_genes,
        fraction_recovered = rr$fraction_recovered,
        stringsAsFactors = FALSE
      )
    }

    run_log[[paste(ct, "celltype", sep = "__")]] <- data.frame(
      Celltype_use = ct,
      Celltype_label = unname(ct_label_map[ct]),
      n_cells = res_ct$n_cells,
      stringsAsFactors = FALSE
    )
  }
}

# =========================================================
# 9. Export merged tables
# =========================================================
all_scores <- bind_rows(score_tables)
all_summaries <- bind_rows(summary_tables)
all_thresholds <- bind_rows(threshold_tables)
all_gene_mappings <- bind_rows(gene_mapping_tables)
all_stats_celllevel <- bind_rows(stats_celllevel_tables)
run_log_df <- bind_rows(run_log)

write.csv(all_scores, file.path(score_tab_dir, sprintf("%s_IL6_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_GFAP_scores_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)
write.csv(all_summaries, file.path(summary_tab_dir, sprintf("%s_IL6_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_GFAP_summary_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)
if (nrow(all_thresholds)) write.csv(all_thresholds, file.path(threshold_tab_dir, sprintf("%s_IL6_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_GFAP_threshold_candidates_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)
write.csv(all_gene_mappings, file.path(mapping_tab_dir, sprintf("%s_IL6_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_GFAP_gene_mapping_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)
if (nrow(all_stats_celllevel)) write.csv(all_stats_celllevel, file.path(stats_cell_tab_dir, sprintf("%s_IL6_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_GFAP_celllevel_stats_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)
write.csv(run_log_df, file.path(run_tab_dir, sprintf("%s_run_log_%s.csv", all_ct_tag, timestamp_tag2)), row.names = FALSE)

# Per-feature aggregated CSVs. CRITICAL: gate each write by `nrow(filter)>0`
# so that a step4 run for a CT that does NOT carry the feature does NOT create
# an empty shadow CSV with the same name as the per-CT/per-feature canonical
# output (which is written inside run_one_feature_set). Without this gate, a
# non-applicable run's empty file would be newer than the applicable run's
# data file and find_score_file in step5 (which picks newest by mtime) would
# read the empty file -> downstream plots break.
.scores_il6   <- all_scores    %>% filter(feature_set == "IL6")
.scores_mod   <- all_scores    %>% filter(feature_set == "MSigDB_Hallmark_IL6_JAK_STAT3_signaling")
.scores_gfap  <- all_scores    %>% filter(feature_set == "GFAP")
.scores_aif1  <- all_scores    %>% filter(feature_set == "AIF1")
.summ_il6     <- all_summaries %>% filter(feature_set == "IL6")
.summ_mod     <- all_summaries %>% filter(feature_set == "MSigDB_Hallmark_IL6_JAK_STAT3_signaling")
.summ_gfap    <- all_summaries %>% filter(feature_set == "GFAP")
.summ_aif1    <- all_summaries %>% filter(feature_set == "AIF1")
if (nrow(.scores_il6))  write.csv(.scores_il6,  file.path(score_tab_dir,   sprintf("%s_IL6_scores_%s.csv",                                                   all_ct_tag, timestamp_tag2)), row.names = FALSE)
if (nrow(.scores_mod))  write.csv(.scores_mod,  file.path(score_tab_dir,   sprintf("main%d_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_scores_%s.csv",            length(main_module_cts), timestamp_tag2)), row.names = FALSE)
if (nrow(.scores_gfap)) write.csv(.scores_gfap, file.path(score_tab_dir,   sprintf("Ast_GFAP_scores_%s.csv",                                                  timestamp_tag2)), row.names = FALSE)
if (nrow(.scores_aif1)) write.csv(.scores_aif1, file.path(score_tab_dir,   sprintf("Mic_AIF1_scores_%s.csv",                                                  timestamp_tag2)), row.names = FALSE)
if (nrow(.summ_il6))    write.csv(.summ_il6,    file.path(summary_tab_dir, sprintf("%s_IL6_summary_%s.csv",                                                   all_ct_tag, timestamp_tag2)), row.names = FALSE)
if (nrow(.summ_mod))    write.csv(.summ_mod,    file.path(summary_tab_dir, sprintf("main%d_MSigDB_Hallmark_IL6_JAK_STAT3_signaling_summary_%s.csv",           length(main_module_cts), timestamp_tag2)), row.names = FALSE)
if (nrow(.summ_gfap))   write.csv(.summ_gfap,   file.path(summary_tab_dir, sprintf("Ast_GFAP_summary_%s.csv",                                                 timestamp_tag2)), row.names = FALSE)
if (nrow(.summ_aif1))   write.csv(.summ_aif1,   file.path(summary_tab_dir, sprintf("Mic_AIF1_summary_%s.csv",                                                 timestamp_tag2)), row.names = FALSE)
rm(.scores_il6, .scores_mod, .scores_gfap, .scores_aif1,
   .summ_il6,   .summ_mod,   .summ_gfap,   .summ_aif1)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  xlsx_path <- file.path(run_tab_dir, sprintf("hbcc_AUCell_score_summary_%s.xlsx", timestamp_tag2))
  openxlsx::write.xlsx(
    list(
      scores_all = all_scores,
      summary_all = all_summaries,
      run_log = run_log_df,
      celltype_counts = celltype_count_table,
      gene_mapping = all_gene_mappings,
      stats_celllevel = if (!nrow(all_stats_celllevel)) data.frame() else all_stats_celllevel
    ),
    file = xlsx_path,
    overwrite = TRUE
  )
}

# =========================================================
# 10. Save lightweight run-level RDS
# =========================================================
run_bundle <- list(
  package_version = as.character(packageVersion("AUCell")),
  timestamp = timestamp_tag2,
  input_rds = hbcc_rds,
  output_root = out_root,
  params = list(
    top_percent_default_multi_gene = top_percent_use,
    top_percent_extended_single_gene = top_percent_use_single,
    aucMaxRank_fraction_multi_gene = top_percent_use,
    aucMaxRank_fraction_single_gene = top_percent_use_single,
    smallestPopPercent = 0.25,
    threshold_method_default = "AUCell_exploreThresholds_selected; tenPercentOfMax replaced by Global_k1 when selected",
    chunk_ncells_use = chunk_ncells_use,
    chunk_mode_note = "Cell-chunked AUCell ranking (memory-bounded). Per-cell rank is intra-cell; chunking is mathematically equivalent to a single AUCell call.",
    group_levels = display_levels,
    feature_sets = feature_sets,
    single_gene_positive_rule = "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0 (rank-and-detection)",
    multi_gene_positive_rule = "AUC score > AUCell selected threshold (default top window aucMaxRank=0.05*G)"
  ),
  celltype_order = ct_order,
  cts_present = cts_present,
  main4_celltypes = main_module_cts,
  celltype_label_map = ct_label_map,
  celltype_count_table = celltype_count_table,
  raw_mapping_check = raw_mapping_check,
  run_log = run_log_df,
  score_tables = score_tables,
  summary_tables = summary_tables,
  threshold_tables = threshold_tables,
  gene_mapping_tables = gene_mapping_tables,
  stats_celllevel_tables = stats_celllevel_tables
)

run_rds <- file.path(step4_dir, sprintf("hbcc_run_%s.rds", timestamp_tag2))
saveRDS(run_bundle, run_rds)

# =========================================================
# 11. Manifest and session information
# =========================================================
manifest_lines <- c(
  paste0("PACKAGE_VERSION: AUCell ", as.character(packageVersion("AUCell"))),
  paste0("INPUT_RDS: ", normalizePath(hbcc_rds, winslash = "/", mustWork = FALSE)),
  paste0("OUT_ROOT: ", normalizePath(out_root, winslash = "/", mustWork = FALSE)),
  paste0("TAB_DIR: ", normalizePath(tab_dir, winslash = "/", mustWork = FALSE)),
  paste0("RDS_DIR: ", normalizePath(rds_dir, winslash = "/", mustWork = FALSE)),
  paste0("RUN_RDS: ", normalizePath(run_rds, winslash = "/", mustWork = FALSE))
)

writeLines(
  manifest_lines,
  con = file.path(step4_dir, sprintf("exports_%s.txt", timestamp_tag2))
)

write.csv(
  data.frame(
    parameter = c(
      "AUCell_version",
      "top_percent_default_multi_gene", "top_percent_extended_single_gene",
      "aucMaxRank_fraction_multi_gene", "aucMaxRank_fraction_single_gene",
      "smallestPopPercent", "threshold_method_default",
      "threshold_method_fallback",
      "single_gene_positive_rule", "multi_gene_positive_rule",
      "chunk_ncells_use",
      "chunk_mode_note"
    ),
    value = c(
      as.character(packageVersion("AUCell")),
      top_percent_use, top_percent_use_single,
      top_percent_use, top_percent_use_single,
      0.25, "AUCell_exploreThresholds_selected",
      "tenPercentOfMax replaced by Global_k1 when AUCell selects it",
      "AUC score > 0 at extended top window (aucMaxRank=0.10*G) OR raw count > 0 (rank-and-detection)",
      "AUC score > AUCell selected threshold (default top window aucMaxRank=0.05*G)",
      as.character(chunk_ncells_use),
      "Cell-chunked AUCell ranking when n_cells > chunk_ncells_use; mathematically equivalent to a single unchunked call (per-cell rank is intra-cell)."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(run_tab_dir, sprintf("analysis_parameters_%s.csv", timestamp_tag2)),
  row.names = FALSE
)

sink(file.path(log_dir, sprintf("sessionInfo_%s.txt", timestamp_tag2)))
print(sessionInfo())
sink()

message("============================================================")
message("HBCC formal analysis (AUCell) finished successfully.")
message("Output root: ", normalizePath(out_root, winslash = "/", mustWork = FALSE))
message("Run RDS: ", normalizePath(run_rds, winslash = "/", mustWork = FALSE))
message("============================================================")
