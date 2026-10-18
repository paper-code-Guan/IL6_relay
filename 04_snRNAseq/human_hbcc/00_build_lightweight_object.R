## =========================================================
## 00_build_lightweight_object.R
##
## Combine 5 per-CT HBCC Seurat RDS files (Ast / ExcitNeu / InhibNeu / Mic /
## OL) into ONE lightweight Seurat object suitable for HBCC step5 in
## "5-CT light mode".
##
## Input data in this workflow:
##   - 5 per-CT RDS at HBCC layout (rownames already gene symbol; metadata
##     uses HBCC names: donor_id, group/Group, Celltype_use)
##   - 1 global UMAP table at HBCC_X_umap.csv.gz with columns
##     (cell, UMAP_1, UMAP_2) covering all ~1.49M cells
##
## Strategy:
##   1. Read each per-CT RDS, force-stamp Celltype_use from filename, copy
##      donor_id->Patient, group/Group->Group, subset features to KEEP_GENES
##      so the merged dgCMatrix never approaches the 2^31-1 NNZ ceiling.
##   2. merge() the 5 lightweight per-CT objects.
##   3. Read the external global UMAP CSV.gz, filter to cells in the merged
##      object, attach as the umap_35 reduction (step5 reads this name).
##   4. saveRDS to HBCC_LIGHTWEIGHT_OUTPUT.
##
## Required inputs are supplied through environment variables:
##   HBCC_RDS_AST
##   HBCC_RDS_EXCITNEU
##   HBCC_RDS_INHIBNEU
##   HBCC_RDS_MIC
##   HBCC_RDS_OL
##   HBCC_GLOBAL_UMAP
##   HBCC_LIGHTWEIGHT_OUTPUT
##   KEEP_GENES          (default: AIF1,GFAP)
##
## NOTE: step5 expects the UMAP reduction key to be 'umap_35' so we always
## attach under that name regardless of upstream naming.
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
set.seed(123)

required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop(name, " must be set.", call. = FALSE)
  }
  value
}

ct_inputs <- list(
  Ast      = required_env("HBCC_RDS_AST"),
  ExcitNeu = required_env("HBCC_RDS_EXCITNEU"),
  InhibNeu = required_env("HBCC_RDS_INHIBNEU"),
  Mic      = required_env("HBCC_RDS_MIC"),
  OL       = required_env("HBCC_RDS_OL")
)
global_umap_path <- required_env("HBCC_GLOBAL_UMAP")

out_path <- required_env("HBCC_LIGHTWEIGHT_OUTPUT")
out_path <- normalizePath(out_path, winslash = "/", mustWork = FALSE)

keep_genes_env <- Sys.getenv("KEEP_GENES", unset = "AIF1,GFAP")
keep_genes <- trimws(unlist(strsplit(keep_genes_env, ",")))
keep_genes <- keep_genes[nzchar(keep_genes)]

# ---------- Validate inputs ----------
nonexistent_rds <- names(ct_inputs)[!file.exists(unlist(ct_inputs))]
if (length(nonexistent_rds)) {
  stop("RDS file(s) missing on disk:\n  ",
       paste(sprintf("%s -> %s", nonexistent_rds, unlist(ct_inputs)[nonexistent_rds]),
             collapse = "\n  "))
}
if (!file.exists(global_umap_path)) {
  stop("Global UMAP file not found: ", global_umap_path)
}

cat("============================================================\n")
cat("00_build_lightweight_object.R\n")
cat("============================================================\n")
cat(sprintf("Output             : %s\n", out_path))
cat(sprintf("Global UMAP source : %s\n", global_umap_path))
cat(sprintf("KEEP_GENES         : %s\n", paste(keep_genes, collapse = ", ")))
cat("Per-CT inputs:\n")
for (ct in names(ct_inputs)) cat(sprintf("  %-10s : %s\n", ct, ct_inputs[[ct]]))
cat("\n")

# ---------- HBCC group canonicalization ----------
canonicalize_hbcc_group <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("Ctr", "CTR", "Ctrl", "CTRL", "Control", "CONTROL", "control",
                  "Normal", "NORMAL", "normal", "Healthy", "HEALTHY", "healthy")] <- "CTR"
  out[out %in% c("Scz", "SCZ", "Sz", "SZ", "sz", "Case", "CASE", "case",
                  "Schizophrenia", "schizophrenia", "SCHIZOPHRENIA")] <- "SCZ"
  out[!out %in% c("CTR", "SCZ")] <- NA_character_
  out
}

# ---------- Per-CT load + lightweight subset ----------
load_one_ct <- function(ct, rds_path) {
  cat(sprintf("[%s] reading %s\n", ct, rds_path))
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) stop(sprintf("[%s] not a Seurat object", ct))

  md <- obj@meta.data

  # --- donor_id -> Patient, group/Group canonicalize -> CTR/SCZ ---
  patient_src <- if ("Patient" %in% colnames(md)) md$Patient
                 else if ("donor_id" %in% colnames(md)) md$donor_id
                 else if ("orig.ident" %in% colnames(md)) md$orig.ident
                 else stop(sprintf("[%s] no Patient / donor_id / orig.ident column", ct))
  obj$Patient <- as.character(patient_src)

  group_src <- if ("Group" %in% colnames(md)) md$Group
               else if ("group" %in% colnames(md)) md$group
               else if ("Diagnosis" %in% colnames(md)) md$Diagnosis
               else stop(sprintf("[%s] no Group / group / Diagnosis column", ct))
  group_canonical <- canonicalize_hbcc_group(group_src)
  if (any(is.na(group_canonical))) {
    bad_levels <- unique(as.character(group_src)[is.na(group_canonical)])
    stop(sprintf("[%s] unrecognized Group values: %s", ct, paste(bad_levels, collapse = ", ")))
  }
  obj$Group <- factor(group_canonical, levels = c("CTR", "SCZ"))

  # --- Force CT id from filename (file-of-record convention) ---
  obj$Celltype_use <- ct
  obj$Celltypes_raw <- ct
  obj$Celltypes <- ct  # also for legacy/step4 compat

  # --- Subset features to keep_genes ---
  feat_present <- intersect(keep_genes, rownames(obj))
  feat_missing <- setdiff(keep_genes, rownames(obj))
  if (length(feat_missing)) {
    cat(sprintf("  [%s] WARNING: KEEP_GENES not found in this RDS: %s\n",
                ct, paste(feat_missing, collapse = ", ")))
  }
  if (!length(feat_present)) {
    feat_present <- rownames(obj)[1]
    cat(sprintf("  [%s] WARNING: zero KEEP_GENES present; placeholder gene '%s' kept\n",
                ct, feat_present))
  }
  obj_lite <- subset(obj, features = feat_present)

  cat(sprintf("  [%s] cells = %d  | CTR = %d  SCZ = %d  | features kept = %s\n",
              ct, ncol(obj_lite),
              sum(obj_lite$Group == "CTR"),
              sum(obj_lite$Group == "SCZ"),
              paste(feat_present, collapse = ",")))

  obj_lite
}

obj_list <- mapply(load_one_ct, names(ct_inputs), unlist(ct_inputs),
                   SIMPLIFY = FALSE, USE.NAMES = TRUE)

# ---------- Cell barcode uniqueness check ----------
all_bc <- unlist(lapply(obj_list, colnames), use.names = FALSE)
add_cell_ids <- NULL
if (length(all_bc) != length(unique(all_bc))) {
  cat(sprintf("[Cell barcodes] %d total, %d unique. Will prefix with CT name during merge.\n",
              length(all_bc), length(unique(all_bc))))
  add_cell_ids <- names(obj_list)
} else {
  cat(sprintf("[Cell barcodes] %d total, all unique across CTs.\n", length(all_bc)))
}

# ---------- Merge ----------
cat("\n[Merge] combining 5 lightweight per-CT obj...\n")
combined <- merge(obj_list[[1]], y = obj_list[-1],
                  add.cell.ids = add_cell_ids,
                  project = "hbcc_5CT_lightweight")

# Seurat 5: collapse multi-layer counts back to 1 layer
if (utils::packageVersion("Seurat") >= "5.0.0") {
  combined <- tryCatch(JoinLayers(combined, assay = "RNA"),
                       error = function(e) {
                         cat(sprintf("[JoinLayers] not applied: %s\n", conditionMessage(e)))
                         combined
                       })
}
cat(sprintf("  merged cells   : %d\n", ncol(combined)))
cat(sprintf("  merged features: %d\n", nrow(combined)))

# ---------- Read external global UMAP ----------
cat(sprintf("\n[UMAP] reading global UMAP from %s\n", global_umap_path))
umap_df <- read.csv(gzfile(global_umap_path), stringsAsFactors = FALSE,
                    check.names = FALSE)

# Tolerate alternate column names
cell_col_candidates <- c("cell", "Cell", "cell_id", "barcode", "Barcode")
u1_candidates <- c("UMAP_1", "umap1", "UMAP1", "X_umap1", "X1")
u2_candidates <- c("UMAP_2", "umap2", "UMAP2", "X_umap2", "X2")
cell_col <- intersect(cell_col_candidates, colnames(umap_df))[1]
u1_col   <- intersect(u1_candidates,        colnames(umap_df))[1]
u2_col   <- intersect(u2_candidates,        colnames(umap_df))[1]
if (any(is.na(c(cell_col, u1_col, u2_col)))) {
  stop(sprintf("Global UMAP CSV missing expected columns. Got: %s. Need cell + UMAP_1 + UMAP_2 (or aliases).",
               paste(colnames(umap_df), collapse = ", ")))
}
cat(sprintf("  cell column = %s, coord columns = %s + %s, rows = %d\n",
            cell_col, u1_col, u2_col, nrow(umap_df)))

# Build embedding matrix indexed by cell ID
emb_mat <- as.matrix(umap_df[, c(u1_col, u2_col), drop = FALSE])
rownames(emb_mat) <- as.character(umap_df[[cell_col]])
colnames(emb_mat) <- c("UMAP_1", "UMAP_2")

# If we prefixed barcodes during merge, the global UMAP IDs (which lack the
# prefix) won't directly index combined cell names. Strip the prefix back to
# look up the embedding row, then re-prefix as the row name.
if (!is.null(add_cell_ids)) {
  combined_bare <- sub("^[^_]+_", "", colnames(combined))
  emb_aligned <- emb_mat[combined_bare, , drop = FALSE]
  rownames(emb_aligned) <- colnames(combined)
} else {
  missing_in_global <- setdiff(colnames(combined), rownames(emb_mat))
  if (length(missing_in_global)) {
    stop(sprintf("Global UMAP missing %d cells from merged obj (e.g. %s ...). Cell IDs do not match.",
                 length(missing_in_global), paste(head(missing_in_global, 3), collapse = ", ")))
  }
  emb_aligned <- emb_mat[colnames(combined), , drop = FALSE]
}

# Always attach under the name 'umap_35' because HBCC step5 hard-codes that
# reduction name.
combined[["umap_35"]] <- CreateDimReducObject(
  embeddings = emb_aligned,
  key        = "UMAP_",
  assay      = DefaultAssay(combined)
)
cat("[UMAP] reduction 'umap_35' attached (key = 'UMAP_')\n")

# ---------- Diagnostic ----------
cat("\n[Metadata summary]\n")
md <- combined@meta.data
cat(sprintf("  total cells     : %d\n", nrow(md)))
cat(sprintf("  Group counts    : %s\n",
            paste(names(table(md$Group)), table(md$Group), sep = "=", collapse = ", ")))
cat(sprintf("  Patient (n_uniq): %d\n", length(unique(md$Patient))))
cat("  Celltype_use breakdown:\n")
ct_tab <- table(md$Celltype_use)
for (k in names(ct_tab)) cat(sprintf("    %-10s = %d\n", k, ct_tab[k]))

emb_check <- Embeddings(combined, "umap_35")
cat(sprintf("  UMAP range      : D1 [%.2f, %.2f]  D2 [%.2f, %.2f]\n",
            min(emb_check[,1]), max(emb_check[,1]),
            min(emb_check[,2]), max(emb_check[,2])))

nnz_counts <- nnzero(GetAssayData(combined, assay = "RNA", layer = "counts"))
cat(sprintf("  counts NNZ      : %s (dgCMatrix max = %s)\n",
            format(nnz_counts, big.mark = ","),
            format(.Machine$integer.max, big.mark = ",")))

# ---------- Save ----------
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(combined, out_path)

cat(sprintf("\n[Done] Saved lightweight obj to: %s\n", out_path))
cat(sprintf("       Size on disk: %.1f MB\n", file.info(out_path)$size / 1024^2))
cat("       Use this RDS as SEURAT_RDS_PATH for HBCC step5.\n")
cat("============================================================\n")
