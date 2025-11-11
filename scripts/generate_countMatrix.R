#!/usr/bin/env Rscript
# ============================================================
# generate_countMatrix.R
#
# Description:
#   Imports transcript-level Kallisto quantifications, aggregates
#   to gene level using tximport, and produces a gene-level count
#   matrix. Only samples listed in metadata_harmonized.csv are
#   included. It also performs sequencing-depth diagnostics,
#   merges with sample metadata, generates summary plots, and
#   runs a preliminary PCA for quality control.
# ============================================================

# -------------------- Load libraries --------------------
# Load required packages for import, aggregation, and visualization
suppressWarnings(suppressMessages({
  library(tximport)
  library(readr)
  library(GenomicFeatures)
  library(dplyr)
  library(data.table)
  library(cli)
  library(AnnotationDbi)
  library(stringr)
  library(ggplot2)
  library(tibble)
}))

# -------------------- Define paths --------------------
# Define input/output file paths and directories
kallisto_dir <- "/home/bavdiu/LONG_COVID_001/data/01-dataflow/01-kallisto_alignment"
count_dir <- "../results/count_matrix"
plot_dir <- "../plots/sequencing_depth"
metadata_file <- "../data/metadata/metadata_harmonized.csv"
gtf_path <- "/srv/opt/genomes/hs/GRCh38_107/Homo_sapiens.GRCh38.107.gtf.gz"

# -------------------- Load metadata --------------------
# Import metadata and verify required columns
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

metadata_raw <- fread(metadata_file)
needed <- c("sample_id", "gse_number")
missing_cols <- setdiff(needed, names(metadata_raw))
if (length(missing_cols) > 0) {
  stop("Metadata is missing columns: ", paste(missing_cols, collapse = ", "))
}

metadata <- metadata_raw |>
  transmute(
    sample_id = trimws(sample_id),
    gse_group = trimws(gse_number)
  ) |>
  filter(!is.na(sample_id), sample_id != "") |>
  distinct() |>
  mutate(gse_group = ifelse(is.na(gse_group) | gse_group == "", "Unknown", gse_group))

cli_alert_info("Loaded metadata with {nrow(metadata)} unique samples.")

# -------------------- Identify sample directories --------------------
# Identify Kallisto output folders matching metadata samples
if (!dir.exists(kallisto_dir)) {
  stop("Kallisto directory not found: ", kallisto_dir)
}

cand <- list.files(kallisto_dir, pattern = "^SRR", full.names = FALSE)
samples_all <- cand[dir.exists(file.path(kallisto_dir, cand))]

samples <- intersect(samples_all, metadata$sample_id)
excluded <- setdiff(samples_all, samples)

if (length(samples) == 0) {
  stop("No matching Kallisto outputs found for metadata samples.")
}
if (length(excluded) > 0) {
  cli_alert_warning("Excluded {length(excluded)} Kallisto directories not listed in metadata.")
}

cli_alert_info("Using {length(samples)} samples that match metadata.")

# -------------------- Build transcript-to-gene mapping --------------------
# Create mapping table from transcripts to genes using GTF annotation
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = keys(txdb, keytype = "TXNAME"),
  columns = "GENEID",
  keytype = "TXNAME"
) |>
  rename(transcript_id = TXNAME, gene_id = GENEID)

# -------------------- Locate Kallisto abundance files --------------------
# Locate all abundance.tsv files for selected samples
files <- file.path(kallisto_dir, samples, "abundance.tsv")
names(files) <- samples

# -------------------- Verify input files --------------------
# Check that all expected abundance.tsv files exist
cli_progress_step("Checking abundance files...")
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  cli_alert_danger("Missing {length(missing_files)} abundance.tsv file(s), e.g.:")
  print(utils::head(missing_files))
  stop("Aborting due to missing files.")
}

# -------------------- Import quantifications --------------------
# Import transcript-level quantifications and summarize to gene level
cli_progress_step("Importing {length(samples)} samples with tximport...")
txi <- tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "no"
)

# -------------------- Save gene-level count matrix --------------------
# Save aggregated gene-level count matrix for downstream analysis
if (!dir.exists(count_dir)) dir.create(count_dir, recursive = TRUE)

count_matrix <- txi$counts |>
  round() |>
  as.data.frame() |>
  data.table::setDT(keep.rownames = "gene_id")

output_file <- file.path(count_dir, "gene_countMatrix.csv")
fwrite(count_matrix, output_file)
cli_alert_success("Saved filtered count matrix to {output_file}")

# -------------------- Sequencing-depth analysis --------------------
# Compute sequencing depth per sample and summarize per GSE group
cli_alert_info("Performing depth analysis...")

sample_depths <- colSums(txi$counts) |>
  tibble::enframe(name = "sample", value = "depth")

sample_depths <- sample_depths |>
  left_join(metadata, by = c("sample" = "sample_id")) |>
  mutate(gse_group = ifelse(is.na(gse_group) | gse_group == "", "Unknown", gse_group))

group_depths <- sample_depths |>
  group_by(gse_group) |>
  summarise(
    total_depth = sum(depth),
    mean_depth = mean(depth),
    median_depth = median(depth),
    sample_count = n(),
    .groups = "drop"
  ) |>
  arrange(desc(total_depth))

fwrite(sample_depths, file.path(count_dir, "sample_depths.csv"))
fwrite(group_depths, file.path(count_dir, "group_depths.csv"))

# -------------------- Generate depth plots --------------------
# Create and export sequencing depth summary plots
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

# Plot 1: Per-sample depth by GSE group
p1 <- ggplot(
  sample_depths,
  aes(
    x = factor(gse_group, levels = unique(group_depths$gse_group)),
    y = depth / 1e6
  )
) +
  geom_boxplot(aes(fill = gse_group), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  scale_y_log10() +
  labs(
    title = "Sequencing Depth per Sample by GSE Group",
    x = "GSE Group",
    y = "Total Reads (Millions, log10)",
    caption = paste("Total samples:", nrow(sample_depths))
  ) +
  plot_theme

# Plot 2: Total depth per GSE group
p2 <- ggplot(
  group_depths,
  aes(
    x = reorder(gse_group, -total_depth),
    y = total_depth / 1e6, fill = gse_group
  )
) +
  geom_col(alpha = 0.8) +
  geom_text(
    aes(label = paste0(round(total_depth / 1e6, 1), "M\n(n=", sample_count, ")")),
    vjust = -0.3, size = 3, lineheight = 0.8
  ) +
  labs(
    title = "Total Sequencing Depth by GSE Group",
    x = "GSE Group",
    y = "Total Reads (Millions)"
  ) +
  plot_theme +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Plot 3: Depth per sample (ordered)
sample_depths <- sample_depths |>
  mutate(sample_index = dplyr::row_number())

p3 <- ggplot(sample_depths, aes(x = sample_index, y = depth / 1e6)) +
  geom_col(fill = "#2c7fb8", alpha = 0.8) +
  labs(
    title = "Sequencing Depth per Sample",
    x = "Sample Index",
    y = "Total Reads (Millions)",
    caption = paste("Total samples:", nrow(sample_depths))
  ) +
  plot_theme

# -------------------- PCA diagnostic --------------------
# Perform PCA on log-transformed counts for quality control
log_counts <- log2(txi$counts + 1)
gene_var <- matrixStats::rowVars(as.matrix(log_counts))
keep <- gene_var > 0
removed_const <- sum(!keep)
log_counts_f <- log_counts[keep, , drop = FALSE]

N_top <- min(5000L, nrow(log_counts_f))
if (nrow(log_counts_f) > N_top) {
  top_idx <- order(gene_var[keep], decreasing = TRUE)[seq_len(N_top)]
  log_counts_f <- log_counts_f[top_idx, , drop = FALSE]
}

if (ncol(log_counts_f) < 2L) {
  stop("PCA requires at least two samples.")
}
if (nrow(log_counts_f) < 2L) {
  stop("PCA requires at least two variable genes.")
}

pca <- prcomp(t(log_counts_f), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:5, drop = FALSE]) |>
  tibble::rownames_to_column("sample") |>
  left_join(dplyr::select(sample_depths, sample, gse_group), by = "sample")

var_expl <- (pca$sdev^2) / sum(pca$sdev^2) * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = gse_group)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = paste0(
      "PCA of Samples (colored by GSE Group)\n",
      removed_const, " constant genes removed; ",
      nrow(log_counts_f), " genes used"
    ),
    x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
    y = paste0("PC2 (", round(var_expl[2], 1), "%)")
  ) +
  theme_minimal()

# -------------------- Save plots --------------------
# Export all diagnostic and depth plots as PDF files
ggsave(file.path(plot_dir, "depth_by_gse_group.pdf"), p1, width = 10, height = 7)
ggsave(file.path(plot_dir, "total_depth_by_gse_group.pdf"), p2, width = 10, height = 7)
ggsave(file.path(plot_dir, "depth_per_sample.pdf"), p3, width = 10, height = 6)
ggsave(file.path(plot_dir, "pca_by_gse_group.pdf"), p_pca, width = 8, height = 6)

cli_alert_success("Saved plots to {plot_dir}")
cli_alert_info("Sample depth summary (millions of reads):")
print(summary(sample_depths$depth / 1e6))
