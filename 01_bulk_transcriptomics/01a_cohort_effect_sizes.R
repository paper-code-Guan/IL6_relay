#!/usr/bin/env Rscript

# Stage 1 of the bulk DLPFC meta-analysis.
#
# Derives per-cohort IL6 effect sizes, standard errors and case-control P values,
# renders Supplementary Figure 1a-c, and writes the aggregate table consumed by
# 01b_meta_analysis.R.
#
# Standard errors are computed on the log-ratio scale (see log_ratio_se below).
# Select the estimator with --se-method:
#
#   delta      first-order delta method (default)
#   bootstrap  stratified nonparametric bootstrap, --bootstrap-reps controls
#              the number of replicates
#   legacy     the pre-correction pooled sd(expr_norm) / sqrt(N); retained only
#              to reproduce previously released numbers
#
# The controlled-access PsychENCODE count matrices and donor metadata are not
# redistributed with this archive. Two input modes are supported:
#
# 1. Raw-count mode reproduces the original preprocessing:
#    count matrix -> edgeR CPM -> log10(CPM + 1) -> IL6 extraction.
#    Each count matrix has genes in rows and donors in columns. Metadata contain
#    columns name and diagnosis, with diagnosis equal to Control or
#    Schizophrenia. A single metadata file may be supplied for all cohorts.
#
# 2. Donor-table mode starts from one table per cohort with columns:
#      group      character, SCZ or CTR
#      expr_norm  numeric, log10(CPM + 1)-transformed IL6 expression
#
# Raw-count mode:
#   Rscript 01a_cohort_effect_sizes.R \
#     --cmc-counts  <CMC_genecounts.matrix> \
#     --hbcc-counts <CMC_HBCC_genecounts.matrix> \
#     --libd-counts <LIBD_scControl_genecounts.matrix> \
#     --metadata    <group.csv> \
#     --output      ../../source_data/bulk_DLPFC_meta_analysis_input.csv \
#     --plot-output-dir output
#
# Donor-table mode:
#   Rscript 01a_cohort_effect_sizes.R \
#     --cmc  <CMC_IL6_donor_level.csv> \
#     --hbcc <HBCC_IL6_donor_level.csv> \
#     --libd <LIBD_IL6_donor_level.csv> \
#     --output ../../source_data/bulk_DLPFC_meta_analysis_input.csv \
#     --plot-output-dir output

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(metafor)
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

read_delimited <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Input table not found: %s", path), call. = FALSE)
  }
  sep <- if (grepl("\\.(tsv|txt|matrix)$", path, ignore.case = TRUE)) "\t" else ","
  read.delim(
    path,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
}

read_donor_table <- function(path) {
  tab <- read_delimited(path)
  missing <- setdiff(c("group", "expr_norm"), colnames(tab))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s",
              path, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  tab$group <- toupper(trimws(as.character(tab$group)))
  tab$expr_norm <- suppressWarnings(as.numeric(tab$expr_norm))
  tab <- tab[
    !is.na(tab$expr_norm) & is.finite(tab$expr_norm) &
      tab$group %in% c("SCZ", "CTR"),
    c("group", "expr_norm"),
    drop = FALSE
  ]
  if (!all(c("SCZ", "CTR") %in% tab$group)) {
    stop(sprintf("%s must contain both SCZ and CTR donors.", path),
         call. = FALSE)
  }
  tab
}

read_counts_and_extract_il6 <- function(count_path, metadata_path) {
  if (!file.exists(count_path)) {
    stop(sprintf("Count matrix not found: %s", count_path), call. = FALSE)
  }
  expr <- read.delim(
    count_path,
    row.names = 1L,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  if (nrow(expr) == 0L || ncol(expr) == 0L) {
    stop(sprintf("Count matrix is empty: %s", count_path), call. = FALSE)
  }
  if (anyDuplicated(colnames(expr))) {
    stop(sprintf("Count matrix has duplicated donor names: %s", count_path),
         call. = FALSE)
  }

  expr_matrix <- as.matrix(expr)
  storage.mode(expr_matrix) <- "numeric"
  if (anyNA(expr_matrix) || any(!is.finite(expr_matrix)) || any(expr_matrix < 0)) {
    stop(sprintf("Count matrix contains invalid numeric values: %s", count_path),
         call. = FALSE)
  }

  il6_id <- "ENSG00000136244"
  gene_ids <- sub("\\..*$", "", rownames(expr_matrix))
  il6_index <- which(gene_ids == il6_id)
  if (length(il6_index) != 1L) {
    stop(
      sprintf("Expected exactly one %s row in %s, found %d.",
              il6_id, count_path, length(il6_index)),
      call. = FALSE
    )
  }

  metadata <- read_delimited(metadata_path)
  missing <- setdiff(c("name", "diagnosis"), colnames(metadata))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s",
              metadata_path, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  metadata$name <- trimws(as.character(metadata$name))
  if (anyDuplicated(metadata$name)) {
    stop(sprintf("Metadata contain duplicated sample names: %s", metadata_path),
         call. = FALSE)
  }
  metadata <- metadata[match(colnames(expr_matrix), metadata$name), , drop = FALSE]
  diagnosis <- trimws(as.character(metadata$diagnosis))
  diagnosis_key <- tolower(diagnosis)
  keep <- !is.na(diagnosis_key) &
    diagnosis_key %in% c("control", "schizophrenia", "ctr", "scz")
  if (!any(keep)) {
    stop(sprintf("No schizophrenia or control donors remain for %s.", count_path),
         call. = FALSE)
  }

  expr_matrix <- expr_matrix[, keep, drop = FALSE]
  diagnosis_key <- diagnosis_key[keep]
  library_sizes <- colSums(expr_matrix)
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop(sprintf("At least one retained donor has a non-positive library size: %s",
                 count_path), call. = FALSE)
  }

  # Supplying full-library sizes to edgeR::cpm avoids constructing a second
  # full CPM matrix while remaining equivalent to cpm(expr)[IL6, ].
  il6_cpm <- edgeR::cpm(
    expr_matrix[il6_index, , drop = FALSE],
    lib.size = library_sizes
  )[1L, ]
  tab <- data.frame(
    group = ifelse(diagnosis_key %in% c("schizophrenia", "scz"), "SCZ", "CTR"),
    expr_norm = log10(as.numeric(il6_cpm) + 1),
    stringsAsFactors = FALSE
  )
  if (!all(c("SCZ", "CTR") %in% tab$group)) {
    stop(sprintf("%s must retain both schizophrenia and control donors.",
                 count_path), call. = FALSE)
  }
  tab
}

# Standard error of log2(mean(SCZ) / mean(CTR)).
#
# The effect size is a unitless log ratio, so its standard error has to be
# expressed on that scale. The pre-correction code used
#   sd(expr_norm) / sqrt(N)
# pooled over both diagnoses, which is the standard error of the grand mean on
# the log10(CPM + 1) scale. That quantity is smaller than the standard error of
# the log ratio by roughly a factor of ln(2) * mean(expr_norm), about twenty-fold
# for IL6, and it is not the standard error of any quantity entering the
# meta-analytic model. It is kept here only as method = "legacy".
log_ratio_se <- function(scz, ctr,
                         method = c("delta", "bootstrap", "legacy"),
                         reps = 2000L,
                         seed = 20260807L) {
  method <- match.arg(method)
  if (method == "legacy") {
    pooled <- c(scz, ctr)
    return(stats::sd(pooled) / sqrt(length(pooled)))
  }

  mean_scz <- mean(scz)
  mean_ctr <- mean(ctr)
  if (!is.finite(mean_scz) || !is.finite(mean_ctr) ||
      mean_scz <= 0 || mean_ctr <= 0) {
    stop(
      "Both group means must be positive to form a log ratio and its error.",
      call. = FALSE
    )
  }

  if (method == "delta") {
    return(
      (1 / log(2)) * sqrt(
        stats::var(scz) / (length(scz) * mean_scz^2) +
          stats::var(ctr) / (length(ctr) * mean_ctr^2)
      )
    )
  }

  set.seed(seed)
  draws <- vapply(seq_len(reps), function(i) {
    a <- mean(sample(scz, length(scz), replace = TRUE))
    b <- mean(sample(ctr, length(ctr), replace = TRUE))
    if (a <= 0 || b <= 0) NA_real_ else log2(a / b)
  }, numeric(1))
  draws <- draws[is.finite(draws)]
  if (length(draws) < 0.9 * reps) {
    warning(
      sprintf("Only %d of %d bootstrap replicates were usable.",
              length(draws), reps),
      call. = FALSE
    )
  }
  if (length(draws) < 2L) {
    stop("Too few usable bootstrap replicates to estimate a standard error.",
         call. = FALSE)
  }
  stats::sd(draws)
}

summarise_cohort <- function(tab,
                             se_method = "delta",
                             bootstrap_reps = 2000L) {
  scz <- tab$expr_norm[tab$group == "SCZ"]
  ctr <- tab$expr_norm[tab$group == "CTR"]
  effect <- log2(mean(scz) / mean(ctr))
  se <- log_ratio_se(scz, ctr, se_method, bootstrap_reps)
  data.frame(
    n_schizophrenia = length(scz),
    n_control = length(ctr),
    mean_schizophrenia = mean(scz),
    mean_control = mean(ctr),
    sd_schizophrenia = stats::sd(scz),
    sd_control = stats::sd(ctr),
    effect_size_log2_expression_ratio = effect,
    standard_error = se,
    standard_error_method = se_method,
    standard_error_legacy = log_ratio_se(scz, ctr, "legacy"),
    # Wald test implied by the effect size and its scale-correct standard error.
    # Reported for transparency only. Figure 1b labels significance on the pooled
    # estimate alone, using the random-effects P value computed in 01b. The
    # per-cohort Wilcoxon statistic below is what Supplementary Figure 1a-c and
    # Supplementary Table 1 report, and it suits the zero-inflated donor-level
    # distributions better than the Wald test.
    p_value_wald = 2 * stats::pnorm(-abs(effect / se)),
    p_value_exact = stats::wilcox.test(
      scz,
      ctr,
      alternative = "two.sided"
    )$p.value,
    stringsAsFactors = FALSE
  )
}

reported_p_value <- function(cohort, p_value) {
  if (cohort == "LIBD" && p_value < 2e-11) {
    return(c(value = 2e-11, bound = "upper"))
  }
  c(value = signif(p_value, 3L), bound = "exact")
}

p_plot_label <- function(cohort, p_value) {
  if (cohort == "LIBD" && p_value < 2e-11) {
    return("italic(P) < 2 %*% 10^-11")
  }
  if (p_value >= 0.001) {
    return(sprintf("italic(P) == %.4f", p_value))
  }
  exponent <- floor(log10(p_value))
  coefficient <- p_value / (10^exponent)
  sprintf("italic(P) == %.2f %%*%% 10^%d", coefficient, exponent)
}

# Plain-text counterpart of p_plot_label for the released statistics table.
p_text_label <- function(cohort, p_value) {
  if (cohort == "LIBD" && p_value < 2e-11) {
    return("P < 2e-11")
  }
  if (p_value >= 0.001) {
    return(sprintf("P = %.4f", p_value))
  }
  exponent <- floor(log10(p_value))
  sprintf("P = %.2fe%d", p_value / (10^exponent), exponent)
}

render_supplementary_panel <- function(tab, cohort, panel, p_value, output_dir) {
  font_family <- Sys.getenv("PLOT_FONT_FAMILY", unset = "Helvetica")
  figure <- ggplot(tab, aes(x = group, y = expr_norm, fill = group, colour = group)) +
    geom_violin(width = 0.9, alpha = 0.8, linewidth = 0.6, colour = NA) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0, seed = 20260807L),
      alpha = 1,
      size = 1.5,
      stroke = 0
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(tab$expr_norm) - 0.05 * diff(range(tab$expr_norm)),
      label = p_plot_label(cohort, p_value),
      parse = TRUE,
      family = font_family,
      size = 15 / ggplot2::.pt
    ) +
    scale_fill_manual(values = c(CTR = "#7FB2D5", SCZ = "#F47F72")) +
    scale_colour_manual(values = c(CTR = "#7FB2D5", SCZ = "#F47F72")) +
    scale_x_discrete(limits = c("CTR", "SCZ")) +
    labs(x = NULL, y = "IL6 Expression") +
    theme_classic(base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 1, colour = "black"),
      axis.title.y = element_text(size = 18, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 18, face = "bold", colour = "black"),
      axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
      axis.ticks = element_line(linewidth = 1, colour = "black"),
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank()
    )

  figure_file <- file.path(
    output_dir,
    sprintf("Supplementary_Fig1%s_%s.pdf", panel, cohort)
  )
  ggsave(
    figure_file,
    figure,
    width = 9,
    height = 9,
    units = "cm",
    device = grDevices::pdf,
    limitsize = FALSE,
    useDingbats = FALSE
  )
  write.csv(
    data.frame(
      cohort = cohort,
      n_schizophrenia = sum(tab$group == "SCZ"),
      n_control = sum(tab$group == "CTR"),
      p_value_exact = p_value,
      p_value_display = p_text_label(cohort, p_value),
      transformation = "log10(CPM + 1)",
      stringsAsFactors = FALSE
    ),
    file.path(
      output_dir,
      sprintf("Supplementary_Fig1%s_%s_statistics.csv", panel, cohort)
    ),
    row.names = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
output_file <- args[["output"]]
if (is.null(output_file)) {
  stop("--output is required.", call. = FALSE)
}

se_method <- args[["se-method"]]
if (is.null(se_method)) {
  se_method <- "delta"
}
se_method <- match.arg(se_method, c("delta", "bootstrap", "legacy"))

bootstrap_reps <- args[["bootstrap-reps"]]
bootstrap_reps <- if (is.null(bootstrap_reps)) {
  2000L
} else {
  as.integer(bootstrap_reps)
}
if (is.na(bootstrap_reps) || bootstrap_reps < 2L) {
  stop("--bootstrap-reps must be an integer of at least 2.", call. = FALSE)
}

donor_files <- list(
  CMC = args[["cmc"]],
  HBCC = args[["hbcc"]],
  LIBD = args[["libd"]]
)
count_files <- list(
  CMC = args[["cmc-counts"]],
  HBCC = args[["hbcc-counts"]],
  LIBD = args[["libd-counts"]]
)
has_donor_mode <- any(!vapply(donor_files, is.null, logical(1)))
has_count_mode <- any(!vapply(count_files, is.null, logical(1)))
if (has_donor_mode && has_count_mode) {
  stop("Do not mix donor-table and raw-count arguments.", call. = FALSE)
}
if (!has_donor_mode && !has_count_mode) {
  stop(
    paste(
      "Supply either --cmc/--hbcc/--libd donor tables or",
      "--cmc-counts/--hbcc-counts/--libd-counts matrices."
    ),
    call. = FALSE
  )
}

if (has_donor_mode) {
  if (any(vapply(donor_files, is.null, logical(1)))) {
    stop("--cmc, --hbcc and --libd are all required in donor-table mode.",
         call. = FALSE)
  }
  cohort_tables <- lapply(donor_files, read_donor_table)
} else {
  if (any(vapply(count_files, is.null, logical(1)))) {
    stop(
      "--cmc-counts, --hbcc-counts and --libd-counts are all required.",
      call. = FALSE
    )
  }
  metadata_files <- list(
    CMC = if (!is.null(args[["cmc-metadata"]])) {
      args[["cmc-metadata"]]
    } else {
      args[["metadata"]]
    },
    HBCC = if (!is.null(args[["hbcc-metadata"]])) {
      args[["hbcc-metadata"]]
    } else {
      args[["metadata"]]
    },
    LIBD = if (!is.null(args[["libd-metadata"]])) {
      args[["libd-metadata"]]
    } else {
      args[["metadata"]]
    }
  )
  if (any(vapply(metadata_files, is.null, logical(1)))) {
    stop(
      paste(
        "Raw-count mode requires --metadata or all three cohort-specific",
        "metadata arguments."
      ),
      call. = FALSE
    )
  }
  cohort_tables <- lapply(names(count_files), function(cohort) {
    read_counts_and_extract_il6(
      count_files[[cohort]],
      metadata_files[[cohort]]
    )
  })
  names(cohort_tables) <- names(count_files)
}

donor_output_dir <- args[["donor-output-dir"]]
if (!is.null(donor_output_dir)) {
  dir.create(donor_output_dir, recursive = TRUE, showWarnings = FALSE)
  for (cohort in names(cohort_tables)) {
    write.csv(
      cohort_tables[[cohort]],
      file.path(donor_output_dir, sprintf("%s_IL6_donor_level.csv", cohort)),
      row.names = FALSE
    )
  }
}

summary_rows <- lapply(names(cohort_tables), function(cohort) {
  row <- summarise_cohort(cohort_tables[[cohort]], se_method, bootstrap_reps)
  cbind(cohort = cohort, row, stringsAsFactors = FALSE)
})
study_data <- do.call(rbind, summary_rows)
rownames(study_data) <- NULL

plot_output_dir <- args[["plot-output-dir"]]
if (!is.null(plot_output_dir)) {
  dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
  panels <- c(CMC = "a", HBCC = "b", LIBD = "c")
  for (cohort in names(cohort_tables)) {
    render_supplementary_panel(
      cohort_tables[[cohort]],
      cohort,
      panels[[cohort]],
      study_data$p_value_exact[study_data$cohort == cohort],
      plot_output_dir
    )
  }
}

fit <- rma(
  yi = effect_size_log2_expression_ratio,
  sei = standard_error,
  method = "REML",
  data = study_data
)
study_data$ci_95_lower <- study_data$effect_size_log2_expression_ratio -
  1.96 * study_data$standard_error
study_data$ci_95_upper <- study_data$effect_size_log2_expression_ratio +
  1.96 * study_data$standard_error
study_data$weight_percent <- as.numeric(weights(fit))

reported <- lapply(seq_len(nrow(study_data)), function(index) {
  reported_p_value(study_data$cohort[[index]], study_data$p_value_exact[[index]])
})
study_data$p_value <- as.numeric(vapply(reported, `[[`, character(1), "value"))
study_data$p_value_bound <- vapply(reported, `[[`, character(1), "bound")

pooled_row <- data.frame(
  cohort = "Random-effects pooled",
  n_schizophrenia = sum(study_data$n_schizophrenia),
  n_control = sum(study_data$n_control),
  mean_schizophrenia = NA_real_,
  mean_control = NA_real_,
  sd_schizophrenia = NA_real_,
  sd_control = NA_real_,
  effect_size_log2_expression_ratio = as.numeric(fit$b),
  standard_error = fit$se,
  standard_error_method = se_method,
  standard_error_legacy = NA_real_,
  ci_95_lower = fit$ci.lb,
  ci_95_upper = fit$ci.ub,
  weight_percent = 100,
  p_value_wald = fit$pval,
  p_value = fit$pval,
  p_value_bound = "exact",
  stringsAsFactors = FALSE
)

column_order <- c(
  "cohort", "n_schizophrenia", "n_control",
  "mean_schizophrenia", "mean_control",
  "sd_schizophrenia", "sd_control",
  "effect_size_log2_expression_ratio", "standard_error",
  "standard_error_method", "standard_error_legacy",
  "ci_95_lower", "ci_95_upper", "weight_percent",
  "p_value_wald", "p_value", "p_value_bound"
)
output <- rbind(study_data[, column_order], pooled_row[, column_order])

dir.create(
  dirname(normalizePath(output_file, mustWork = FALSE)),
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(output, output_file, row.names = FALSE)

message(sprintf(
  paste(
    "Cohort effect sizes written to %s",
    "Standard-error method: %s",
    "Pooled estimate %.7f, 95%% CI %.7f to %.7f, P = %.7g,",
    "I-squared = %.2f%%, tau-squared = %.4f, Q = %.2f (df = %d), Q P = %.4g",
    sep = "\n"
  ),
  output_file, se_method, as.numeric(fit$b), fit$ci.lb, fit$ci.ub, fit$pval,
  fit$I2, fit$tau2, fit$QE, fit$k - 1L, fit$QEp
))
