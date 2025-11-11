#!/usr/bin/env Rscript
# WGCNA_modules.R
#
# Description:
# Performs Weighted Gene Co-expression Network Analysis (WGCNA)
# on VST-normalized expression data (no batch correction or outlier filtering).
#
# Inputs:
#   - results/diff_exp_analysis/VST_data.RDS  (DESeq2 output)
#   - data/metadata/metadata_harmonized.csv
#
# Outputs:
#   - results/wgcna/
#   - plots/wgcna/
#
# Notes:
#   - Grey modules are hidden in plots but retained in tables.
#   - The workflow includes: soft-thresholding, network construction,
#     module detection, eigengene extraction, and module–trait correlations.

options(stringsAsFactors = FALSE)
set.seed(123)

suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(here)
    library(Biobase)
    library(WGCNA)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(matrixStats)
})

# ---------------------- Paths and setup ----------------------
here::i_am("scripts/WGCNA_modules.R")
allowWGCNAThreads(nThreads = max(1, parallel::detectCores() - 1))
has_limma <- requireNamespace("limma", quietly = TRUE)

# Input files (from DESeq2)
VST_RDS <- here::here("results", "diff_exp_analysis", "VST_data_limma_corrected.RDS")
META_CSV <- here::here("data", "metadata", "metadata_harmonized.csv")

# Output directories
OUTDIR <- here::here("results", "wgcna")
OUTDIR_plots <- here::here("plots", "wgcna")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTDIR_plots, recursive = TRUE, showWarnings = FALSE)

# ---------------------- Load data and align ----------------------
stopifnot(file.exists(VST_RDS), file.exists(META_CSV))

vsd <- readRDS(VST_RDS)

if (inherits(vsd, "SummarizedExperiment") || inherits(vsd, "DESeqTransform")) {
    expr <- SummarizedExperiment::assay(vsd)
} else if (is.matrix(vsd)) {
    expr <- vsd
} else {
    stop("Unknown object type in VST RDS file: ", class(vsd))
}

expr <- as.matrix(expr)
rownames(expr) <- sub("\\..*$", "", rownames(expr))

meta <- read_csv(META_CSV, show_col_types = FALSE)
stopifnot(all(c("sample_id", "condition") %in% names(meta)))
meta$condition <- factor(meta$condition,
    levels = c("healthy", "acuteCOVID", "recovered", "longCOVID", "exposedNoSymptoms")
)

# Align samples between expression and metadata
common <- intersect(colnames(expr), meta$sample_id)
if (length(common) < 4) stop("Too few overlapping samples between expression and metadata.")
expr <- expr[, common, drop = FALSE]
meta <- meta |>
    filter(sample_id %in% common) |>
    arrange(match(sample_id, colnames(expr)))
stopifnot(identical(colnames(expr), meta$sample_id))

message("Expression matrix: ", nrow(expr), " genes × ", ncol(expr), " samples")

# ---------------------- Gene filtering (variance-based) ----------------------
targetGenes <- 8000
if (nrow(expr) > targetGenes) {
    v <- matrixStats::rowVars(expr)
    expr <- expr[names(sort(v, decreasing = TRUE))[1:targetGenes], , drop = FALSE]
}
saveRDS(expr, file.path(OUTDIR, "WGCNA_input_matrix.rds"))

# ---------------------- Sample clustering (QC) ----------------------
datExpr0 <- t(expr)
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]

pdf(file.path(OUTDIR_plots, "WGCNA_sampleClustering_rawExpr.pdf"), width = 16, height = 6)
sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering (raw VST expression)", xlab = "", sub = "", cex = 0.5)
dev.off()

# ---------------------- Soft-threshold power selection ----------------------
powers <- c(1:20, seq(22, 30, 2))
sft <- pickSoftThreshold(datExpr0,
    powerVector = powers, networkType = "signed",
    corFnc = "bicor", corOptions = list(maxPOutliers = 0.05), verbose = 5
)

pdf(file.path(OUTDIR_plots, "WGCNA_softThresholding_powers.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)", ylab = "Scale-Free Topology Model Fit (signed R²)",
    type = "n", main = "Scale independence"
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers, cex = 0.9, col = "blue"
)
abline(h = 0.80, col = "red", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)", ylab = "Mean connectivity",
    type = "n", main = "Mean connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "blue")
dev.off()

r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
cand <- powers[r2 >= 0.80]
pickedPower <- if (length(cand) > 0) min(cand) else 6
message("Selected soft-threshold power: ", pickedPower)

# ---------------------- Network construction and module detection ----------------------
net <- blockwiseModules(
    datExpr0,
    power = pickedPower,
    networkType = "signed",
    TOMType = "signed",
    corType = "bicor",
    maxBlockSize = 50000,
    deepSplit = 2,
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3
)

moduleColors <- net$colors
MEs <- orderMEs(net$MEs)

# ---------------------- Gene dendrogram and module colors ----------------------
pdf(file.path(OUTDIR_plots, "WGCNA_geneDendrogram_moduleColors.pdf"), width = 12, height = 6)
plotDendroAndColors(
    net$dendrograms[[1]], labels2colors(moduleColors),
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
)
dev.off()

# ---------------------- Save module assignments ----------------------
ens_ids <- colnames(datExpr0)
ann <- AnnotationDbi::select(org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL", columns = "SYMBOL"
) |> distinct(ENSEMBL, .keep_all = TRUE)
symb <- ann$SYMBOL[match(ens_ids, ann$ENSEMBL)]

tab_assign <- tibble(ENSEMBL = ens_ids, SYMBOL = symb, module = labels2colors(moduleColors))
write_csv(tab_assign, file.path(OUTDIR, "WGCNA_gene_module_assignment.csv"))

mod_sizes <- sort(table(labels2colors(moduleColors)), decreasing = TRUE)
mod_sizes_df <- data.frame(module = names(mod_sizes), size = as.integer(mod_sizes))
write_csv(mod_sizes_df, file.path(OUTDIR, "WGCNA_module_sizes.csv"))

# ---------------------- Module eigengenes and trait correlations ----------------------
MEs_noGrey <- MEs[, !grepl("^MEgrey$", colnames(MEs)), drop = FALSE]

# Export module eigengenes
MEs_mat <- as.data.frame(MEs_noGrey) %>%
    tibble::rownames_to_column("sample_id")
readr::write_csv(MEs_mat, file.path(OUTDIR, "WGCNA_module_eigengenes_samples.csv"))
message("[INFO] Wrote module eigengenes: ", file.path(OUTDIR, "WGCNA_module_eigengenes_samples.csv"))

# Build trait matrix (categorical + numeric)
traits <- model.matrix(~ 0 + condition, data = meta)
colnames(traits) <- gsub("^condition", "", colnames(traits))

# Add optional covariates
if ("age" %in% names(meta)) {
    traits <- cbind(traits, age = as.numeric(meta$age))
}
if ("gender" %in% names(meta)) {
    g <- tolower(meta$gender)
    g <- ifelse(g %in% c("male", "m"), 1,
        ifelse(g %in% c("female", "f"), 0, NA)
    )
    traits <- cbind(traits, gender = g)
}
traits <- traits[, apply(traits, 2, function(x) !all(is.na(x)) && var(x, na.rm = TRUE) > 0), drop = FALSE]

corME <- bicor(MEs_noGrey, traits, use = "pairwise.complete.obs", maxPOutliers = 0.05)
pME <- corPvalueStudent(corME, nSamples = nrow(MEs_noGrey))

pdf(file.path(OUTDIR_plots, "WGCNA_moduleTrait_heatmap_with_covariates.pdf"), width = 12, height = 8)
labeledHeatmap(
    Matrix = corME,
    xLabels = colnames(traits),
    yLabels = gsub("^ME", "", colnames(MEs_noGrey)),
    ySymbols = gsub("^ME", "", colnames(MEs_noGrey)),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = signif(corME, 2),
    colorBarLabelPosition = 0.5
)
dev.off()

write_csv(
    as.data.frame(corME) |> rownames_to_column("ME"),
    file.path(OUTDIR, "WGCNA_ME_trait_correlations_with_covariates.csv")
)
write_csv(
    as.data.frame(pME) |> rownames_to_column("ME"),
    file.path(OUTDIR, "WGCNA_ME_trait_pvalues_with_covariates.csv")
)

# ---------------------- Hub gene analysis (kME) ----------------------
kME <- bicor(datExpr0, MEs_noGrey, use = "pairwise.complete.obs", maxPOutliers = 0.05)
kME_df <- as.data.frame(kME) |>
    rownames_to_column("ENSEMBL") |>
    left_join(tab_assign, by = "ENSEMBL")
write_csv(kME_df, file.path(OUTDIR, "WGCNA_geneModuleMembership_kME.csv"))

get_top_hubs <- function(kME_df, topN = 50) {
    out <- list()
    mods <- setdiff(unique(kME_df$module), "grey")
    for (mod in mods) {
        col_me <- paste0("ME", mod)
        if (!col_me %in% names(kME_df)) next
        sub <- kME_df |>
            filter(module == mod) |>
            mutate(kME_abs = abs(.data[[col_me]])) |>
            arrange(desc(kME_abs)) |>
            slice_head(n = topN)
        out[[mod]] <- sub
    }
    bind_rows(out, .id = "module_color")
}
top_hubs <- get_top_hubs(kME_df)
if (nrow(top_hubs) > 0) {
    write_csv(top_hubs, file.path(OUTDIR, "WGCNA_topHubGenes_perModule.csv"))
}

# ---------------------- Session info ----------------------
sink(file.path(OUTDIR, "sessionInfo.txt"))
sessionInfo()
sink()

message("WGCNA completed successfully.")
message("Results: ", normalizePath(OUTDIR))
message("Plots:   ", normalizePath(OUTDIR_plots))
