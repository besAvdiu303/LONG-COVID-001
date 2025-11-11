#!/usr/bin/env Rscript
# ============================================================
# DESeq2_analysis.R
#
# Description:
# Differential gene expression analysis using DESeq2.
# This script loads raw count data and corresponding metadata,
# performs normalization and differential expression modeling
# while accounting for batch effects (GSE numbers),
# and exports normalized counts, contrasts, and VST-transformed data.
# ============================================================

suppressWarnings(suppressMessages({
    library(DESeq2)
    library(limma)
    library(Biobase)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(dplyr)
    library(stringr)
    library(tibble)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(matrixStats)
}))

# ---------------------- Paths ----------------------
datadir <- "../data"
results <- "../results"
resultsdir <- "../results/diff_exp_analysis"
if (!dir.exists(resultsdir)) dir.create(resultsdir, recursive = TRUE)

counts_path <- file.path(results, "count_matrix", "gene_countMatrix.csv")
metadata_path <- file.path(datadir, "metadata", "metadata_harmonized.csv")

# ---------------------- Load counts ----------------------
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

counts_matrix <- readr::read_csv(counts_path, show_col_types = FALSE) |>
    column_to_rownames("gene_id") |>
    as.matrix()

# Remove samples with zero total counts
counts_matrix <- counts_matrix[, colSums(counts_matrix) > 0, drop = FALSE]

# ---------------------- Load metadata ----------------------
if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path)

metadata <- readr::read_csv(metadata_path, col_names = TRUE, show_col_types = FALSE)

# Ensure required columns are present
required_cols <- c("sample_id", "gse_number", "condition")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
    stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Harmonize metadata and set factor levels (reference = healthy)
phenoData <- metadata |>
    mutate(
        sample_id  = trimws(sample_id),
        gse_number = trimws(gse_number),
        condition  = trimws(condition)
    ) |>
    mutate(
        condition = factor(
            condition,
            levels = c("healthy", "acuteCOVID", "longCOVID", "exposedNoSymptoms", "recovered")
        )
    ) |>
    mutate(rowname = sample_id) |>
    column_to_rownames("rowname")

# Remove samples with missing or invalid entries (avoids NA issues in DESeq2)
valid <- !is.na(phenoData$condition) & phenoData$condition != "" &
    !is.na(phenoData$gse_number) & phenoData$gse_number != ""
if (any(!valid)) {
    message(
        "Dropping ", sum(!valid), " sample(s) with missing condition or gse_number: ",
        paste(rownames(phenoData)[!valid], collapse = ", ")
    )
    phenoData <- phenoData[valid, , drop = FALSE]
}

# Display sample distribution
print(table(condition = phenoData$condition))
print(table(batch = phenoData$gse_number))

# ---------------------- Gene annotation (ENSEMBL) ----------------------
# Used for annotation of results but not for filtering.
database <- org.Hs.eg.db::org.Hs.eg.db
featureData <- AnnotationDbi::select(
    database,
    columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"),
    keys    = keys(database, keytype = "ENSEMBL"),
    keytype = "ENSEMBL"
) |>
    mutate(FEATUREID = ENSEMBL) |>
    group_by(FEATUREID) |>
    summarise(
        ENSEMBL  = paste(unique(ENSEMBL), collapse = ", "),
        ENTREZID = paste(unique(ENTREZID), collapse = ", "),
        SYMBOL   = paste(unique(SYMBOL), collapse = ", "),
        GENENAME = paste(unique(GENENAME), collapse = ", "),
        .groups  = "drop"
    ) |>
    mutate(rowname = FEATUREID) |>
    column_to_rownames("rowname")

# ---------------------- Match samples and genes ----------------------
# Ensure identical sample order between counts and metadata.
common_samples <- intersect(colnames(counts_matrix), rownames(phenoData))
if (length(common_samples) < 4) {
    stop("Too few overlapping samples between count matrix and metadata.")
}

ordered_samples <- colnames(counts_matrix)[colnames(counts_matrix) %in% common_samples]
counts_matrix <- counts_matrix[, ordered_samples, drop = FALSE]
phenoData <- phenoData[ordered_samples, , drop = FALSE]

if (nrow(counts_matrix) < 2) stop("Too few genes after filtering.")
if (nrow(phenoData) < 2) stop("Too few samples after alignment.")

# ---------------------- Build DESeq2 dataset ----------------------
phenoData$gse_number <- factor(phenoData$gse_number)
if ("healthy" %in% levels(phenoData$condition)) {
    phenoData$condition <- stats::relevel(phenoData$condition, ref = "healthy")
}

dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData   = phenoData,
    design    = ~ gse_number + condition
)

# ---------------------- Filtering and model fitting ----------------------
# Remove genes with very low counts across all samples
keep_features <- rowSums(counts(dds)) >= 10
dds <- dds[keep_features, ]

# Fit the DESeq2 model
dds <- DESeq(dds)

# ---------------------- Export normalized counts ----------------------
norm_counts <- counts(dds, normalized = TRUE) |>
    as.data.frame() |>
    tibble::rownames_to_column("ENSEMBL")

readr::write_csv(norm_counts, file.path(resultsdir, "normalized_counts.csv"))

# Save DESeq2 object for downstream use
saveRDS(dds, file = file.path(resultsdir, "DDS_subset.RDS"))

# ---------------------- Differential expression contrasts ----------------------
# Define pairwise comparisons between biological conditions
desired_levels <- c("healthy", "acuteCOVID", "longCOVID", "exposedNoSymptoms", "recovered")
present_levels <- intersect(desired_levels, levels(colData(dds)$condition))

contrast_pairs <- list(
    acuteCOVID_vs_healthy           = c("condition", "acuteCOVID", "healthy"),
    longCOVID_vs_healthy            = c("condition", "longCOVID", "healthy"),
    exposedNoSymptoms_vs_healthy    = c("condition", "exposedNoSymptoms", "healthy"),
    recovered_vs_healthy            = c("condition", "recovered", "healthy"),
    longCOVID_vs_acuteCOVID         = c("condition", "longCOVID", "acuteCOVID"),
    exposedNoSymptoms_vs_acuteCOVID = c("condition", "exposedNoSymptoms", "acuteCOVID"),
    recovered_vs_acuteCOVID         = c("condition", "recovered", "acuteCOVID"),
    longCOVID_vs_exposedNoSymptoms  = c("condition", "longCOVID", "exposedNoSymptoms"),
    recovered_vs_exposedNoSymptoms  = c("condition", "recovered", "exposedNoSymptoms"),
    longCOVID_vs_recovered          = c("condition", "longCOVID", "recovered"),
    acuteCOVID_vs_recovered         = c("condition", "acuteCOVID", "recovered")
)

# Keep only contrasts where both conditions exist in the dataset
contrast_pairs <- Filter(function(v) all(v[2:3] %in% present_levels), contrast_pairs)

# Compute results for all valid contrasts
results_list <- lapply(contrast_pairs, function(ctr) {
    res <- results(dds, contrast = ctr)
    res
})

# Save all DESeq2 contrast results
saveRDS(results_list, file = file.path(resultsdir, "DESeq2_results.RDS"))

# ---------------------- Variance Stabilizing Transformation (VST) ----------------------
# Calculate and export VST-transformed counts for downstream analyses (e.g., PCA, clustering)
message("Performing variance stabilizing transformation (VST)...")
vst_data <- vst(dds, blind = FALSE)

# Save VST object for downstream analysis and visualization
saveRDS(vst_data, file = file.path(resultsdir, "VST_data.RDS"))

# Export VST expression matrix
vst_mat <- assay(vst_data) |>
    as.data.frame() |>
    tibble::rownames_to_column("ENSEMBL")

readr::write_csv(vst_mat, file.path(resultsdir, "vst_matrix_raw.csv"))
message("Saved raw VST-transformed matrix (unadjusted for batch).")

# ============================================================
#  limma batch correction
# ============================================================
message("Applying limma batch correction on VST matrix...")

vsd_mat <- assay(vst_data)
batch <- vst_data$gse_number
design <- model.matrix(~ 0 + vst_data$condition)

vsd_limma <- removeBatchEffect(vsd_mat, batch = batch, design = design)

# Save limma-corrected matrix (same naming as ComBat version)
saveRDS(vsd_limma, file = file.path(resultsdir, "VST_data_limma_corrected.RDS"))
readr::write_csv(
    as.data.frame(vsd_limma) |>
        tibble::rownames_to_column("ENSEMBL"),
    file.path(resultsdir, "vst_matrix_limma_corrected.csv")
)

message("Saved limma batch-corrected VST matrix.")

message("Differential expression analysis completed successfully.")
message("Results saved to: ", normalizePath(resultsdir))
