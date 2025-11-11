#!/usr/bin/env Rscript
# ============================================================
# DESeq2_visualization.R
#
# Description:
#   Generates QC and result visualizations for DESeq2 analysis:
#   - Volcano plots for all contrasts
#   - Combined volcano overview (2x2)
#   - MA plot for main contrast
#   - PCA on VST-transformed data
#
#   Modified so that volcano_data_<contrast>.csv files
#   include only genes that pass BOTH padj & log2FC cutoffs.
# ============================================================

# -------------------- Load libraries --------------------
# Load plotting, annotation, and DESeq2-related packages
suppressWarnings(suppressMessages({
  library(EnhancedVolcano)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
}))

# -------------------- Config --------------------
# Define significance thresholds and default main contrast
P_CUTOFF <- 0.05
FC_CUTOFF <- 2.0
MAIN_CONTRAST <- "acuteCOVID_vs_healthy"

# -------------------- Paths --------------------
# Set up input/output directories and verify required files
results_root <- "../results"
deseq2_dir <- file.path(results_root, "diff_exp_analysis")
plot_root <- "../plots/diff_exp"
volcano_dir <- file.path(plot_root, "volcano")
extras_dir <- file.path(plot_root, "extras")

dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extras_dir, recursive = TRUE, showWarnings = FALSE)

res_rds <- file.path(deseq2_dir, "DESeq2_results.RDS")
vsd_rds <- file.path(deseq2_dir, "VST_data.RDS")

if (!file.exists(res_rds)) stop("Missing DESeq2 results: ", res_rds)
if (!file.exists(vsd_rds)) stop("Missing VST_data.RDS: ", vsd_rds)

res_list <- readRDS(res_rds)
vsd <- readRDS(vsd_rds)

# -------------------- Annotation helper --------------------
# Function to map ENSEMBL IDs to gene symbols
ens_to_symbol <- function(ens_ids) {
  AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ens_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
}

# -------------------- Data preparation --------------------
# Prepare DESeq2 results for plotting and export
prep_df <- function(res_obj) {
  df <- as.data.frame(res_obj)
  df$ENSEMBL <- rownames(df)
  df$SYMBOL <- ens_to_symbol(df$ENSEMBL)

  df <- df %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(SYMBOL) & SYMBOL != "")

  # Cap extreme values for better visualization
  df$padj_capped <- pmax(df$padj, 1e-40)
  df$log2FoldChange_capped <- pmax(pmin(df$log2FoldChange, 20), -20)

  df
}

# -------------------- Volcano plots --------------------
# Function to generate a single volcano plot
make_volcano <- function(df, title_text) {
  EnhancedVolcano(
    toptable = df,
    lab = df$SYMBOL,
    x = "log2FoldChange_capped",
    y = "padj_capped",
    pCutoff = P_CUTOFF,
    FCcutoff = FC_CUTOFF,
    title = title_text,
    subtitle = NULL,
    caption = NULL,
    xlab = bquote(~ Log[2] ~ "fold change"),
    ylab = bquote(~ -Log[10] ~ adjusted ~ italic(P)),

    # --- Custom colors and legend ---
    col = c("grey30", "forestgreen", "royalblue", "red2"),
    legendLabels = c(
      "Not significant",
      "Pass log2FC cutoff",
      "Pass p-value cutoff",
      "Significant (both)"
    ),
    legendPosition = "right",
    legendLabSize = 16,

    # --- Visual settings ---
    pointSize = 1.8,
    colAlpha = 0.9,
    labSize = 4.5
  ) +
    coord_cartesian(xlim = c(-25, 25), ylim = c(0, 40)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank()
    ) +
    guides(color = guide_legend(title = NULL))
}

# -------------------- Generate individual volcano plots --------------------
# Loop through all contrasts, generate volcano plots, and export CSVs of significant genes
plots_list <- list()
for (nm in names(res_list)) {
  df <- prep_df(res_list[[nm]])
  p <- make_volcano(df, gsub("_", " ", nm))
  plots_list[[nm]] <- p

  # Save volcano plot (PDF)
  ggsave(
    file.path(volcano_dir, paste0("volcano_", nm, ".pdf")),
    p,
    width = 7, height = 6, device = cairo_pdf
  )
  message("Saved volcano plot for: ", nm)

  # Save CSV containing only significant genes
  df_sig <- df %>%
    filter(
      !is.na(padj),
      padj < P_CUTOFF,
      abs(log2FoldChange) > log2(FC_CUTOFF)
    )

  csv_path <- file.path(deseq2_dir, paste0("volcano_data_", nm, ".csv"))
  write_csv(df_sig, csv_path)
  message(
    "Saved filtered volcano data CSV (significant only) for: ", nm,
    " -> ", csv_path, " (", nrow(df_sig), " genes)"
  )
}

# -------------------- Combined 2x2 Volcano layout --------------------
# Combine four selected contrasts into a single 2x2 volcano panel
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
    tag_prefix = "",
    caption = NULL
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 18),
    plot.tag = element_text(size = 28, face = "bold"),
    plot.tag.position = c(0.1, 0.99)
  ) +
    theme(plot.background = element_blank())

ggsave(file.path(volcano_dir, "combined_2x2_volcano.pdf"),
  combo,
  width = 14, height = 10, device = cairo_pdf
)
message("Saved combined volcano 2x2 plot.")

# -------------------- MA plot --------------------
# Generate MA plot for the main contrast defined above
main_key <- if (MAIN_CONTRAST %in% names(res_list)) MAIN_CONTRAST else names(res_list)[1]
res_main <- res_list[[main_key]]

pdf(file.path(extras_dir, paste0("MA_plot_", main_key, ".pdf")), width = 7, height = 6)
plotMA(res_main, ylim = c(-2, 2))
dev.off()
message("Saved simple MA plot for: ", main_key)

# -------------------- PCA comparison: Raw vs. corrected --------------------
# Create PCA plots for raw VST data and limma batch-corrected VST matrix
# (skipped if corrected file not found or metadata missing)
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(sva)))

vsd_rds_corrected <- file.path(deseq2_dir, "VST_data_limma_corrected.RDS")

SHOW_LABELS <- FALSE
X_LIMIT <- c(-150, 150)
Y_LIMIT <- c(-150, 150)

if (file.exists(vsd_rds_corrected) && all(c("condition", "gse_number") %in% colnames(colData(vsd)))) {
  message("Generating side-by-side PCA plots (Raw vs. corrected)...")

  # ---- Helper: PCA for DESeq2 object ----
  # Computes PCA coordinates using DESeq2 plotPCA()
  make_pca_deseq <- function(vsd_obj, color_var, shape_var = NULL, title_text) {
    pcaData <- plotPCA(vsd_obj, intgroup = unique(na.omit(c(color_var, shape_var))), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaData$sample <- rownames(pcaData)

    pcaData$PC1 <- pmax(pmin(pcaData$PC1, X_LIMIT[2]), X_LIMIT[1])
    pcaData$PC2 <- pmax(pmin(pcaData$PC2, Y_LIMIT[2]), Y_LIMIT[1])

    aes_args <- aes_string(x = "PC1", y = "PC2", color = color_var)
    if (!is.null(shape_var)) {
      aes_args <- aes_string(x = "PC1", y = "PC2", color = color_var, shape = shape_var)
    }

    p <- ggplot(pcaData, aes_args) +
      geom_point(size = 5, alpha = 0.9) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_cartesian(xlim = X_LIMIT, ylim = Y_LIMIT) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      ggtitle(title_text) +
      guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

    if (SHOW_LABELS) {
      p <- p + geom_text_repel(aes(label = sample), size = 3.5, max.overlaps = 20)
    }
    return(p)
  }

  # ---- Helper: PCA for matrix ----
  # Runs PCA manually on a matrix (limma-corrected data)
  make_pca_matrix <- function(mat, meta, color_var, shape_var = NULL, title_text) {
    pca_res <- prcomp(t(mat), scale. = TRUE)
    percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2])
    pca_df <- data.frame(pca_res$x, meta)

    pca_df$PC1 <- pmax(pmin(pca_df$PC1, X_LIMIT[2]), X_LIMIT[1])
    pca_df$PC2 <- pmax(pmin(pca_df$PC2, Y_LIMIT[2]), Y_LIMIT[1])

    aes_args <- aes_string(x = "PC1", y = "PC2", color = color_var)
    if (!is.null(shape_var)) {
      aes_args <- aes_string(x = "PC1", y = "PC2", color = color_var, shape = shape_var)
    }

    p <- ggplot(pca_df, aes_args) +
      geom_point(size = 5, alpha = 0.9) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_cartesian(xlim = X_LIMIT, ylim = Y_LIMIT) +
      theme_minimal(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      ggtitle(title_text) +
      guides(color = guide_legend(title = NULL), shape = guide_legend(title = NULL))

    if (SHOW_LABELS) {
      p <- p + geom_text_repel(aes(label = rownames(meta)), size = 3.5, max.overlaps = 20)
    }
    return(p)
  }

  # Load limma-corrected VST matrix and corresponding metadata
  vsd_corrected_obj <- readRDS(vsd_rds_corrected)

  if (inherits(vsd_corrected_obj, "ExpressionSet")) {
    vsd_corrected <- Biobase::exprs(vsd_corrected_obj)
    meta <- as.data.frame(Biobase::pData(vsd_corrected_obj))
  } else if (is.matrix(vsd_corrected_obj)) {
    vsd_corrected <- vsd_corrected_obj
    meta <- as.data.frame(colData(vsd))
  } else {
    stop("Unknown format of limma-corrected object: expected matrix or ExpressionSet.")
  }

  # Generate and save PCA plots (color by GSE and condition)
  p_raw_gse <- make_pca_deseq(vsd, "gse_number", "condition", "Raw VST")
  p_corrected_gse <- make_pca_matrix(vsd_corrected, meta, "gse_number", "condition", "Corrected VST")

  p_combined_gse <- (p_raw_gse | p_corrected_gse) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A", tag_suffix = ")") &
    theme(
      legend.position = "right",
      plot.tag = element_text(face = "bold", size = rel(1.2))
    )

  ggsave(file.path(extras_dir, "PCA_raw_vs_corrected_colorGSE.pdf"),
    p_combined_gse,
    width = 14, height = 8, device = cairo_pdf
  )
  message("Saved PCA comparison (color = GSE, shape = condition).")

  p_raw_cond <- make_pca_deseq(vsd, "condition", NULL, "Raw VST")
  p_corrected_cond <- make_pca_matrix(vsd_corrected, meta, "condition", NULL, "Batch corrected VST")

  p_combined_cond <- (p_raw_cond | p_corrected_cond) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right")

  ggsave(file.path(extras_dir, "PCA_raw_vs_corrected_colorCondition.pdf"),
    p_combined_cond,
    width = 14, height = 8, device = cairo_pdf
  )
  message("Saved PCA comparison (color = condition).")
} else if (!file.exists(vsd_rds_corrected)) {
  message("No corrected VST found — generating only raw PCA plot.")
} else {
  message("Missing 'condition' or 'gse_number' in colData(vsd); skipping PCA plots.")
}

# -------------------- Done --------------------
# End of visualization pipeline
message("DESeq2 visualization complete. All plots saved under: ", normalizePath(plot_root))
