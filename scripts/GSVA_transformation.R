#!/usr/bin/env Rscript
# ============================================================
# GSVA_transformation.R
#
# Description:
# Performs Gene Set Variation Analysis (GSVA) on the
# variance-stabilized (VST) counts after limma batch correction
# (no outlier filtering), followed by differential pathway analysis
# using limma. Generates publication-quality volcano plots,
# combined 2×2 volcano layouts, and QC visualizations
# (heatmap, PCA, optional UMAP).
# ============================================================

set.seed(123)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(here)
  library(Biobase)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(GSVA)
  library(BiocParallel)
  library(limma)
  library(EnhancedVolcano)
  library(patchwork)
  library(ggrepel)
  library(cowplot)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
})

# ---------------------- Paths ----------------------
here::i_am("scripts/GSVA_transformation.R")
message("here() root: ", here::here())

# Input files (no batch correction)
VST_RDS <- here::here("results", "diff_exp_analysis", "VST_data_limma_corrected.RDS")
META_CSV <- here::here("data", "metadata", "metadata_harmonized.csv")

# Output directories
RES_DIR <- here::here("results", "gsva")
PLOT_DIR <- here::here("plots", "gsva")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------- Load VST data ----------------------
if (!file.exists(VST_RDS)) stop("Missing VST_data.RDS at: ", VST_RDS)
vsd <- readRDS(VST_RDS)

if (inherits(vsd, "SummarizedExperiment") || inherits(vsd, "DESeqTransform")) {
  expr <- SummarizedExperiment::assay(vsd)
} else if (is.matrix(vsd)) {
  expr <- vsd
} else {
  stop("Unknown object type in VST RDS file: ", class(vsd))
}

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"
rownames(expr) <- sub("\\..*$", "", rownames(expr))
expr <- expr[!rowSums(is.na(expr)) > 0, , drop = FALSE]

message("Loaded VST expression matrix: ", nrow(expr), " genes × ", ncol(expr), " samples")

# ---------------------- Load metadata and align samples ----------------------
meta <- readr::read_csv(META_CSV, show_col_types = FALSE)
stopifnot(all(c("sample_id", "condition") %in% names(meta)))
if (!("gse_number" %in% names(meta))) meta$gse_number <- NA_character_

meta$condition <- factor(
  meta$condition,
  levels = c("healthy", "acuteCOVID", "recovered", "longCOVID", "exposedNoSymptoms")
)

common <- intersect(colnames(expr), meta$sample_id)
if (length(common) < 4) stop("Too few overlapping samples between expression and metadata.")
expr <- expr[, common, drop = FALSE]
meta <- meta |>
  filter(sample_id %in% common) |>
  arrange(match(sample_id, colnames(expr)))
stopifnot(identical(colnames(expr), meta$sample_id))

# ---------------------- WikiPathways gene sets ----------------------
species <- "Homo sapiens"

get_wikipathways <- function(species = "Homo sapiens") {
  w <- tryCatch(
    msigdbr::msigdbr(species = species, collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
    error = function(e) NULL
  )
  if (!is.null(w) && nrow(w) > 0) {
    return(w)
  }
  all_df <- msigdbr::msigdbr(species = species)
  dplyr::filter(all_df, grepl("WIKIPATHWAYS", gs_name, ignore.case = TRUE))
}

msig_wp <- get_wikipathways(species)
stopifnot(nrow(msig_wp) > 0)
message("WikiPathways entries: ", nrow(msig_wp))

genes0 <- rownames(expr)
prop_ens <- mean(grepl("^ENS[A-Z]*\\d", genes0))
ID_SPACE <- if (is.finite(prop_ens) && prop_ens > 0.6) "ENSEMBL" else "SYMBOL"
message("Detected ID space: ", ID_SPACE)

# Build gene set list
build_sets <- function(df, id_space, expr_ids, keep_min = 10, keep_max = 1000) {
  if (id_space == "ENSEMBL") {
    if ("ensembl_gene" %in% names(df)) {
      lst <- split(sub("\\..*$", "", df$ensembl_gene), df$gs_name)
    } else {
      map <- AnnotationDbi::select(org.Hs.eg.db,
        keys = unique(df$gene_symbol),
        keytype = "SYMBOL", columns = "ENSEMBL"
      ) |>
        distinct(SYMBOL, .keep_all = TRUE)
      lst <- split(map$ENSEMBL[match(df$gene_symbol, map$SYMBOL)], df$gs_name)
    }
  } else {
    lst <- split(df$gene_symbol, df$gs_name)
  }
  lst <- lapply(lst, function(v) unique(na.omit(intersect(v, expr_ids))))
  lst <- lst[vapply(lst, function(v) length(v) >= keep_min && length(v) <= keep_max, logical(1))]
  lst
}

gene_sets <- build_sets(msig_wp, ID_SPACE, expr_ids = genes0)
if (length(gene_sets) == 0) stop("No valid gene sets retained.")
message("Gene sets retained: ", length(gene_sets))

# ---------------------- GSVA computation ----------------------
BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))
param <- GSVA::gsvaParam(exprData = expr, geneSets = gene_sets, kcdf = "Gaussian")
gsva_mat <- GSVA::gsva(param = param)
if (inherits(gsva_mat, "SummarizedExperiment")) gsva_mat <- SummarizedExperiment::assay(gsva_mat)
stopifnot(is.matrix(gsva_mat))

# Save GSVA results
gsva_eset <- ExpressionSet(assayData = gsva_mat)
pd <- meta |> column_to_rownames("sample_id")
pData(gsva_eset) <- pd
out_rds <- file.path(RES_DIR, "GSVA_scores.ExpressionSet.RDS")
saveRDS(gsva_eset, out_rds)
write_csv(
  as.data.frame(gsva_mat) |> rownames_to_column("pathway"),
  file.path(RES_DIR, "GSVA_scores.csv")
)
message("GSVA scores saved to: ", normalizePath(out_rds))


# Differential Pathway Activity (limma, robust)


message("Running limma differential pathway analysis...")

cond <- meta$condition
valid <- !is.na(cond) & cond != ""
expr_ok <- gsva_mat[, valid, drop = FALSE]
cond_ok <- droplevels(cond[valid])

# Ensure enough replicates per group
tab <- table(cond_ok)
message("Samples per condition:\n", paste(names(tab), tab, sep = "=", collapse = ", "))
small_groups <- names(tab)[tab < 2]
if (length(small_groups) > 0) {
  warning("These groups have <2 samples and will be dropped: ", paste(small_groups, collapse = ", "))
  keep <- !(cond_ok %in% small_groups)
  expr_ok <- expr_ok[, keep, drop = FALSE]
  cond_ok <- droplevels(cond_ok[keep])
}

# Build design matrix
X <- model.matrix(~ 0 + cond_ok)
colnames(X) <- gsub("^cond_ok", "", colnames(X))

fit <- lmFit(expr_ok, design = X)

# Define contrasts that actually exist in data
pairs <- list(
  acuteCOVID_vs_healthy   = c("acuteCOVID", "healthy"),
  acuteCOVID_vs_recovered = c("acuteCOVID", "recovered"),
  recovered_vs_healthy    = c("recovered", "healthy"),
  longCOVID_vs_recovered  = c("longCOVID", "recovered")
)
pairs_keep <- pairs[vapply(pairs, function(ab) all(ab %in% colnames(X)), logical(1))]
if (length(pairs_keep) == 0) stop("No valid contrasts available after filtering.")
contr_defs <- vapply(pairs_keep, function(ab) sprintf("%s - %s", ab[1], ab[2]), character(1))

Cmat <- makeContrasts(contrasts = contr_defs, levels = X)
fit2 <- contrasts.fit(fit, contrasts = Cmat)
fit2 <- eBayes(fit2, trend = TRUE)

message("Limma model fitted with ", ncol(Cmat), " valid contrasts.")


# GSVA Volcano Plots


message("Generating GSVA volcano plots...")

p_cut <- 0.05
fc_cut <- 0.3
volcano_pdf_dir <- file.path(PLOT_DIR, "volcano")
dir.create(volcano_pdf_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function: clean pathway names
clean_label <- function(x) {
  x <- gsub("^WP_", "", x)
  gsub("_", " ", x)
}

# Function to create individual volcano plots
make_volcano_gsva <- function(df, title_text) {
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title_text, "(no data)")))
  }
  
  df$adj.P.Val_capped <- pmax(df$adj.P.Val, 1e-40)
  df$LABEL <- clean_label(df$pathway)
  df$log10P <- -log10(df$adj.P.Val_capped)
  
  df <- df %>%
    mutate(
      significance = case_when(
        abs(logFC) >= fc_cut & adj.P.Val < p_cut ~ "Pass p-value and GSVA score cutoffs",
        abs(logFC) >= fc_cut                     ~ "Pass GSVA score cutoff",
        adj.P.Val < p_cut                        ~ "Pass p-value cutoff",
        TRUE                                     ~ "Not significant"
      ),
      significance = factor(
        significance,
        levels = c(
          "Not significant",
          "Pass p-value cutoff",
          "Pass p-value and GSVA score cutoffs"
        )
      )
    )
  
  label_df <- df %>%
    filter(significance == "Pass p-value and GSVA score cutoffs") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = 5)
  
  ggplot(df, aes(x = logFC, y = log10P)) +
    geom_point(aes(color = significance), alpha = 0.9, size = 1.8) +
    scale_color_manual(
      values = c(
        "Not significant"                        = "grey30",
        "Pass p-value cutoff"                    = "royalblue",
        "Pass p-value and GSVA score cutoffs"    = "red2"
      ),
      drop = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = LABEL),
      size = 3.6, color = "black",
      max.overlaps = Inf,
      box.padding = grid::unit(0.6, "lines"),
      point.padding = grid::unit(0.5, "lines"),
      segment.color = "grey50", segment.size = 0.4,
      min.segment.length = 0, show.legend = FALSE,
      force = 12, force_pull = 0.002, max.time = 10, max.iter = 10000,
      nudge_x = 0.1, nudge_y = 0.2, direction = "both"
    ) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey50") +
    coord_cartesian(xlim = c(-4, 4), ylim = c(0, 40)) +
    labs(
      title = gsub("_", " ", title_text),
      x = "GSVA score",
      y = expression(-log[10] ~ adjusted ~ italic(P)),
      color = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

# Generate volcano plots
ttables <- list()
plots_list <- list()

for (cn in colnames(fit2$coefficients)) {
  cname <- gsub(" - ", "_vs_", cn)
  message("Processing contrast: ", cname)

  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P") |>
    tibble::rownames_to_column("pathway")
  
  tt_sig <- tt |>
    dplyr::filter(abs(logFC) >= fc_cut & adj.P.Val < p_cut)
  
  if (nrow(tt_sig) == 0) {
    message("No significant pathways for: ", cname)
  } else {
    message(nrow(tt_sig), " significant pathways for: ", cname)
  }
  
  ttables[[cname]] <- tt_sig
  readr::write_csv(tt_sig, file.path(RES_DIR, paste0("GSVA_limma_", cname, "_significant.csv")))
  
  p <- make_volcano_gsva(tt, gsub("_", " ", cname))
  plots_list[[cname]] <- p
  
  ggsave(file.path(volcano_pdf_dir, paste0("volcano_", cname, ".pdf")),
         p, width = 7, height = 6, device = cairo_pdf)
  
  message("Saved volcano plot for: ", cname)
}


# Combined 2×2 Volcano Layout


want_keys <- c(
  "acuteCOVID_vs_healthy",
  "recovered_vs_healthy",
  "acuteCOVID_vs_recovered",
  "longCOVID_vs_recovered"
)

get_or_blank <- function(key) {
  if (!is.null(plots_list[[key]])) {
    plots_list[[key]] +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  } else {
    ggplot() +
      theme_void() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  }
}

combo <- ((get_or_blank("acuteCOVID_vs_healthy") | get_or_blank("acuteCOVID_vs_recovered")) /
  (get_or_blank("recovered_vs_healthy") | get_or_blank("longCOVID_vs_recovered"))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ")",
    caption = NULL
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.tag = element_text(size = 28, face = "bold"),
    plot.tag.position = c(0.1, 0.99)
  ) +
    theme(plot.background = element_blank())

ggsave(file.path(volcano_pdf_dir, "GSVA_combined_2x2_volcano.pdf"),
  combo,
  width = 14, height = 10, device = cairo_pdf
)

message("Saved combined GSVA 2×2 volcano plot.")


# GSVA Heatmaps


gsva_mat <- as.matrix(gsva_mat)
storage.mode(gsva_mat) <- "numeric"
gsva_mat[!is.finite(gsva_mat)] <- 0

rownames(gsva_mat) <- gsub("^WP_", "", rownames(gsva_mat))
rownames(gsva_mat) <- gsub("_", " ", rownames(gsva_mat))

meta_df <- as.data.frame(meta)
stopifnot(all(c("sample_id", "condition", "gse_number") %in% names(meta_df)))

common <- intersect(colnames(gsva_mat), meta_df$sample_id)
if (length(common) < 3) stop("Not enough overlapping samples between GSVA and metadata.")

gsva_mat <- gsva_mat[, common, drop = FALSE]
meta_df <- meta_df[match(common, meta_df$sample_id), , drop = FALSE]
stopifnot(identical(colnames(gsva_mat), meta_df$sample_id))

ann_col_df <- data.frame(
  Condition = factor(meta_df$condition),
  GSE = factor(meta_df$gse_number)
)
rownames(ann_col_df) <- meta_df$sample_id

cond_levels <- levels(ann_col_df$Condition)
ncond <- length(cond_levels)
if (ncond <= 1) {
  cond_colors <- setNames("grey70", cond_levels)
} else if (ncond == 2) {
  cond_colors <- setNames(c("#1f78b4", "#e31a1c"), cond_levels)
} else {
  cond_colors <- setNames(RColorBrewer::brewer.pal(max(3, ncond), "Set2")[1:ncond], cond_levels)
}

gse_levels <- levels(ann_col_df$GSE)
ngse <- length(gse_levels)
if (ngse == 1) {
  gse_colors <- setNames("grey60", gse_levels)
} else {
  pal_gse <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(ngse)
  gse_colors <- setNames(pal_gse, gse_levels)
}

ann_colors <- list(
  Condition = cond_colors,
  GSE = gse_colors
)

v <- matrixStats::rowVars(gsva_mat)
topk <- head(order(v, decreasing = TRUE), 50)
z <- t(scale(t(gsva_mat[topk, , drop = FALSE])))
z[!is.finite(z)] <- 0

ord <- order(ann_col_df$Condition, ann_col_df$GSE)
z <- z[, ord, drop = FALSE]
ann_col_df <- ann_col_df[ord, , drop = FALSE]

pheatmap::pheatmap(
  z,
  annotation_col    = ann_col_df,
  annotation_colors = ann_colors,
  filename          = file.path(PLOT_DIR, "GSVA_top50_heatmap.pdf"),
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  fontsize_row      = 7,
  fontsize_col      = 8,
  border_color      = NA,
  treeheight_row    = 25,
  annotation_legend = TRUE,
  width             = 10,
  height            = 8
)

message("Saved GSVA_top50_heatmap.pdf")


# QC: UMAP and PCA on GSVA scores


if (requireNamespace("uwot", quietly = TRUE)) {
  library(uwot)
  library(ggplot2)
  library(ggrepel)

  set.seed(123)

  umap_input <- t(gsva_mat)
  umap_input <- scale(umap_input)
  umap_input <- umap_input[, apply(umap_input, 2, var) > 0, drop = FALSE]

  emb <- uwot::umap(
    umap_input,
    n_neighbors = 15,
    min_dist = 0.3,
    metric = "euclidean"
  )

  df_umap <- as.data.frame(emb)
  colnames(df_umap) <- c("UMAP1", "UMAP2")
  df_umap$sample_id <- rownames(umap_input)

  df_umap <- df_umap |>
    dplyr::left_join(meta, by = "sample_id")

  if ("gse_number" %in% colnames(meta)) {
    gp2 <- ggplot(df_umap, aes(UMAP1, UMAP2, color = gse_number, shape = condition)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal(base_size = 16) +
      labs(
        x = "UMAP1", y = "UMAP2"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)
      )

    ggsave(
      file.path(PLOT_DIR, "GSVA_UMAP_GSE.pdf"),
      gp2,
      width = 8, height = 6
    )
  }

  message("Saved GSVA UMAP plots (colored by condition and optionally by GSE).")
} else {
  message("UMAP skipped — install with: install.packages('uwot')")
}

pca <- prcomp(t(gsva_mat), center = TRUE, scale. = TRUE)
ve <- (pca$sdev^2) / sum(pca$sdev^2) * 100

df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("sample_id") |>
  dplyr::left_join(meta[, c("sample_id", "condition")], by = "sample_id")

p_pca <- ggplot(df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_minimal(base_size = 16) +
  labs(
    x = paste0("PC1 (", round(ve[1], 1), "% variance)"),
    y = paste0("PC2 (", round(ve[2], 1), "% variance)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank()
  )

ggsave(file.path(PLOT_DIR, "GSVA_PCA.pdf"), p_pca, width = 8, height = 6)
message("Saved PCA plot on GSVA scores.")

# ---------------------- Session info ----------------------
sink(file.path(RES_DIR, "sessionInfo.txt"))
sessionInfo()
sink()
if (requireNamespace("devtools", quietly = TRUE)) {
  writeLines(
    capture.output(devtools::session_info()),
    file.path(RES_DIR, "sessionInfo_full.txt")
  )
}

message("GSVA transformation complete.")
message("Results: ", normalizePath(RES_DIR))
message("Plots:   ", normalizePath(PLOT_DIR))
