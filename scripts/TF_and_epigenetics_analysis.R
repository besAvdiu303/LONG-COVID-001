#!/usr/bin/env Rscript
# TF_and_epigenetics_analysis.R
#
# Description:
# Performs transcription factor (TF) activity inference using
# DoRothEA and VIPER on VST-normalized gene expression data.
# Includes pairwise comparisons, limma modeling, WikiPathways
# enrichment, TF–pathway correlation analysis, and visualization
# of master regulators and epigenetic modules.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(Biobase)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dorothea)
  library(viper)
  library(ggplot2)
  library(broom)
  library(limma)
  library(msigdbr)
  library(clusterProfiler)
  library(ggpubr)
  library(patchwork)
  library(stringr)
  library(pheatmap)
  library(igraph)
  library(ggraph)
  library(dplyr)
  library(rstatix)
  library(tidyr)
  library(tibble)
  library(grid)
})

# ---------------------- Helper functions ----------------------
check_alignment <- function(expr, meta) {
  stopifnot(all(meta$sample_id %in% colnames(expr)))
  meta <- meta[match(colnames(expr), meta$sample_id), ]
  stopifnot(all(meta$sample_id == colnames(expr)))
  return(meta)
}

# ---------------------- Paths and directories ----------------------
VST_RDS <- here("results", "diff_exp_analysis", "VST_data_limma_corrected.RDS")
META_FILE <- here("data", "metadata", "metadata_harmonized.csv")

OUT_DIR <- here("results", "DoRothEA_VIPER")
PLOT_DIR <- here("plots", "DoRothEA_VIPER")
ENRICH_DIR <- file.path(OUT_DIR, "enrichment")
LIMMA_DIR <- file.path(OUT_DIR, "limma")
PLOT_ENRICH_DIR <- file.path(PLOT_DIR, "enrichment")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LIMMA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------- Load expression data ----------------------
message("[INFO] Loading VST expression data...")
vsd <- readRDS(VST_RDS)

if (inherits(vsd, "SummarizedExperiment") || inherits(vsd, "DESeqTransform")) {
  expr <- SummarizedExperiment::assay(vsd)
} else if (is.matrix(vsd)) {
  expr <- vsd
} else {
  stop("Unknown object type in VST RDS file: ", class(vsd))
}

expr <- as.matrix(expr)
mode(expr) <- "numeric"
message("[INFO] Expression matrix dimensions: ", paste(dim(expr), collapse = " × "))

# ---------------------- Map Ensembl IDs to gene symbols ----------------------
message("[INFO] Mapping Ensembl IDs to gene symbols...")
ensg_ids <- sub("\\..*$", "", rownames(expr))
sym_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensg_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

expr <- expr[!is.na(sym_map), ]
rownames(expr) <- sym_map[!is.na(sym_map)]

# Collapse duplicated gene symbols
expr <- as.data.frame(expr) %>%
  rownames_to_column("Gene") %>%
  group_by(Gene) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

message("[INFO] Retained ", nrow(expr), " unique genes after mapping.")

# ---------------------- Load metadata ----------------------
message("[INFO] Loading metadata...")
meta <- read_csv(META_FILE, show_col_types = FALSE)
stopifnot(all(c("sample_id", "condition") %in% colnames(meta)))

meta <- meta %>%
  filter(sample_id %in% colnames(expr)) %>%
  mutate(condition = factor(
    condition,
    levels = c("healthy", "acuteCOVID", "recovered", "longCOVID", "exposedNoSymptoms")
  ))

meta <- check_alignment(expr, meta)
message("[INFO] Metadata loaded: ", nrow(meta), " samples.")

# ---------------------- Load DoRothEA regulons and run VIPER ----------------------
message("[INFO] Loading DoRothEA regulons (confidence A–C)...")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
regulon_list <- dorothea::df2regulon(regulons)
message("[INFO] Regulons loaded: ", length(regulon_list), " TFs")

message("[INFO] Running VIPER to infer TF activities...")
tf_activities <- viper(expr, regulon_list, verbose = TRUE)
tf_activities <- as.data.frame(tf_activities)
message("[INFO] VIPER matrix: ", nrow(tf_activities), " TFs × ", ncol(tf_activities), " samples")

saveRDS(tf_activities, file.path(OUT_DIR, "VIPER_TF_activity.RDS"))

# ---------------------- Pairwise t-tests and barplots ----------------------
comparisons <- list(
  acute_vs_healthy = c("acuteCOVID", "healthy"),
  acute_vs_recovered = c("acuteCOVID", "recovered"),
  recovered_vs_healthy = c("recovered", "healthy"),
  long_vs_recovered = c("longCOVID", "recovered")
)

for (cmp in names(comparisons)) {
  groups <- comparisons[[cmp]]
  g1 <- groups[1]
  g2 <- groups[2]
  message("[INFO] Comparison: ", g1, " vs ", g2)

  idx <- meta$condition %in% groups
  dat <- tf_activities[, idx, drop = FALSE]
  cond <- meta$condition[idx]

  res <- apply(dat, 1, function(x) {
    tt <- try(t.test(x[cond == g1], x[cond == g2]), silent = TRUE)
    if (inherits(tt, "try-error")) {
      return(NULL)
    }
    broom::tidy(tt) %>%
      dplyr::select(estimate1, estimate2, statistic, p.value)
  })
  res <- res[!sapply(res, is.null)]
  res <- bind_rows(res, .id = "TF") %>%
    mutate(logFC = estimate1 - estimate2, p_adj = p.adjust(p.value, "fdr")) %>%
    arrange(p_adj)

  out_csv <- file.path(OUT_DIR, paste0("DoRothEA_VIPER_", cmp, ".csv"))
  write_csv(res, out_csv)

  top_res <- res %>% slice_min(p_adj, n = 20)
  p <- ggplot(top_res, aes(x = reorder(TF, logFC), y = logFC, fill = -log10(p_adj))) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(
      title = paste("Top TF activities:", g1, "vs", g2),
      x = "Transcription Factor",
      y = "Delta Activity (VIPER)"
    )

  ggsave(file.path(PLOT_DIR, paste0("DoRothEA_VIPER_", cmp, ".pdf")), p, width = 8, height = 6)
}

message("[INFO] Pairwise VIPER comparisons complete.")

# ---------------------- Combine pairwise barplots (2×2 grid) ----------------------
plots_list <- list()
for (cmp in names(comparisons)) {
  csv_path <- file.path(OUT_DIR, paste0("DoRothEA_VIPER_", cmp, ".csv"))
  if (!file.exists(csv_path)) next
  res <- read_csv(csv_path, show_col_types = FALSE)
  top_res <- res %>% slice_min(p_adj, n = 20)
  g1 <- comparisons[[cmp]][1]
  g2 <- comparisons[[cmp]][2]

  p <- ggplot(top_res, aes(x = reorder(TF, logFC), y = logFC, fill = -log10(p_adj))) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = paste(g1, "vs", g2),
      x = "Transcription Factor",
      y = "logFC Activity (VIPER)",
      fill = expression(-log[10] ~ adjusted ~ italic(P))
    )

  plots_list[[cmp]] <- p
}

combo <- (
  (plots_list$acute_vs_healthy | plots_list$acute_vs_recovered) /
    (plots_list$recovered_vs_healthy | plots_list$long_vs_recovered)
) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 24, face = "bold"))

ggsave(file.path(PLOT_DIR, "DoRothEA_VIPER_combined_2x2_individual_legends.pdf"),
  combo,
  width = 14, height = 10, device = cairo_pdf
)

message("[INFO] Combined VIPER plot saved.")

# ---------------------- Limma and WikiPathways enrichment ----------------------
message("[INFO] Starting extended analysis with limma and WikiPathways enrichment...")

meta <- check_alignment(tf_activities, meta)
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)

fit <- lmFit(tf_activities, design)
contrasts <- makeContrasts(
  acuteCOVID - healthy,
  acuteCOVID - recovered,
  recovered - healthy,
  longCOVID - recovered,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contrasts))

contrast_names <- colnames(fit2$coefficients)
limma_tables <- lapply(contrast_names, function(cn) {
  df <- topTable(fit2, coef = cn, number = Inf, sort.by = "P") %>%
    rownames_to_column("TF") %>%
    rename(p_adj = adj.P.Val, p_value = P.Value)
  return(df)
})
names(limma_tables) <- contrast_names

wp_df <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:WIKIPATHWAYS"
)

for (cn in names(limma_tables)) {
  df <- limma_tables[[cn]]
  sig_TFs <- df %>%
    filter(p_adj < 0.05) %>%
    pull(TF)
  if (length(sig_TFs) < 5) next

  ego <- enricher(
    gene = sig_TFs,
    TERM2GENE = wp_df %>% select(gs_name, gene_symbol),
    universe = rownames(tf_activities),
    pAdjustMethod = "fdr"
  )
  if (is.null(ego) || nrow(ego@result) == 0) next

  enrich_tbl <- ego@result %>%
    arrange(p.adjust) %>%
    rename(Pathway = Description, FDR = p.adjust, Count = Count) %>%
    mutate(
      Pathway = gsub("^WP_", "", Pathway),
      Pathway = gsub("_", " ", Pathway)
    )

  write_csv(enrich_tbl, file.path(ENRICH_DIR, paste0("WikiPathways_enrichment_", cn, ".csv")))

  top_terms <- enrich_tbl %>% slice_min(FDR, n = 20)
  clean_title <- gsub("_", " vs ", cn)

  p_enrich <- ggplot(top_terms, aes(x = reorder(Pathway, -log10(FDR)), y = -log10(FDR), fill = Count)) +
    geom_col(width = 0.7) +
    coord_flip() +
    theme_minimal(base_size = 18) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
    labs(title = clean_title, x = NULL, y = expression(-log[10] ~ FDR))

  ggsave(file.path(PLOT_ENRICH_DIR, paste0("DoRothEA_WikiPathways_", cn, ".pdf")),
    p_enrich,
    width = 15, height = 13
  )
}

# ---------------------- Master regulators and TF–pathway networks ----------------------
message("[INFO] Computing TF–pathway correlations and epigenetic modules...")

GSVA_FILE <- here("results", "gsva", "GSVA_scores.csv")
if (!file.exists(GSVA_FILE)) stop("GSVA file not found: ", GSVA_FILE)

gsva_mat <- read_csv(GSVA_FILE, show_col_types = FALSE) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

common_samples <- intersect(colnames(tf_activities), colnames(gsva_mat))
stopifnot(length(common_samples) >= 4)

tf_activities <- tf_activities[, common_samples, drop = FALSE]
gsva_mat <- gsva_mat[, common_samples, drop = FALSE]

cor_mat <- cor(t(tf_activities), t(gsva_mat), method = "spearman", use = "pairwise.complete.obs")

corr_threshold <- 0.4
tf_degree <- apply(cor_mat, 1, function(x) sum(abs(x) >= corr_threshold))
master_regulators <- sort(tf_degree, decreasing = TRUE)
mr_table <- tibble(TF = names(master_regulators), n_correlated_pathways = master_regulators)
write_csv(mr_table, file.path(OUT_DIR, "DoRothEA_VIPER_master_regulators.csv"))

top_mr <- head(mr_table, 15)
p <- ggplot(top_mr, aes(x = reorder(TF, n_correlated_pathways), y = n_correlated_pathways)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)) +
  labs(title = "Top Master Regulators (TFs)", x = "Transcription Factor", y = "Correlated Pathways >= |0.4|")
ggsave(file.path(PLOT_DIR, "DoRothEA_master_regulators.pdf"), p, width = 7, height = 6)

# ---------------------- Individual TF–pathway subnetworks ----------------------
message("[INFO] Building individual TF–Pathway subnetworks...")

NETWORK_DIR <- file.path(PLOT_DIR, "networks")
dir.create(NETWORK_DIR, showWarnings = FALSE, recursive = TRUE)

corr_threshold <- 0.5
n_top_tfs <- 7
n_top_paths <- 10

top_tfs <- head(names(sort(tf_degree, decreasing = TRUE)), n_top_tfs)

for (tf in top_tfs) {
  cor_vec <- cor_mat[tf, ]
  cor_df <- tibble(Pathway = names(cor_vec), cor = cor_vec) %>%
    filter(abs(cor) >= corr_threshold) %>%
    arrange(desc(abs(cor))) %>%
    slice_head(n = n_top_paths) %>%
    mutate(
      Pathway = gsub("^WP_", "", Pathway),
      Pathway = gsub("_", " ", Pathway),
      Pathway = stringr::str_wrap(Pathway, 30)
    )
  if (nrow(cor_df) == 0) next

  edges <- cor_df %>%
    mutate(TF = tf) %>%
    select(TF, Pathway, cor)

  nodes <- tibble(
    name = c(tf, edges$Pathway),
    type = c("TF", rep("Pathway", nrow(edges)))
  )

  g <- igraph::graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  V(g)$degree <- igraph::degree(g)

  p_net <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(color = cor, alpha = abs(cor))) +
    geom_node_point(aes(color = type, shape = type, size = degree)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.8) +
    scale_color_manual(values = c(TF = "firebrick3", Pathway = "steelblue")) +
    scale_shape_manual(values = c(TF = 21, Pathway = 22)) +
    scale_edge_color_gradient2(low = "blue", high = "red", mid = "grey90", midpoint = 0) +
    scale_size_continuous(range = c(3, 8)) +
    theme_void() +
    ggtitle(paste0(tf, " – correlated pathways (|ρ| ≥ ", corr_threshold, ")"))

  ggsave(file.path(NETWORK_DIR, paste0("TF_network_", tf, ".pdf")),
    p_net,
    width = 10, height = 8
  )
}

message("[INFO] Saved ", n_top_tfs, " individual TF subnetworks in: ", normalizePath(NETWORK_DIR))

# ---------------------- Epigenetic GSVA subset and plots ----------------------
epi_terms <- c(
  "ACETYLATION", "UBIQUITINATION", "METHYLATION", "HISTONE", "CHROMATIN", "EPIGENETIC", "NUCLEOSOME", "KDM", "KMT", "PRMT", "demethylase",
  "COHESIN", "SIRTUIN", "POLYCOMB", "DNA_REPAIR", "PRC2", "ASSEMLBY", "methyltransferase", "Oncogenic", "H1-5"
)

epi_idx <- grepl(paste(epi_terms, collapse = "|"), rownames(gsva_mat), ignore.case = TRUE)
gsva_epi <- gsva_mat[epi_idx, , drop = FALSE]
message("[INFO] Found ", sum(epi_idx), " epigenetic-related pathways.")

df_epi <- as.data.frame(t(gsva_epi)) |>
  rownames_to_column("sample_id") |>
  left_join(meta[, c("sample_id", "condition")], by = "sample_id") |>
  pivot_longer(-c(sample_id, condition), names_to = "Pathway", values_to = "Score") |>
  mutate(
    Pathway = gsub("^WP_", "", Pathway),
    Pathway = gsub("_", " ", Pathway),
    Pathway = str_wrap(Pathway, 25)
  )

df_epi <- df_epi %>%
  filter(condition != "exposedNoSymptoms") %>%
  droplevels()

my_comparisons <- list(
  c("healthy", "acuteCOVID"),
  c("recovered", "acuteCOVID"),
  c("recovered", "longCOVID")
)

stat.test <- df_epi %>%
  group_by(Pathway) %>%
  group_modify(~ compare_means(
    Score ~ condition,
    data = .x,
    method = "wilcox.test",
    comparisons = my_comparisons
  )) %>%
  ungroup() %>%
  filter(p < 0.05) %>%
  mutate(
    p.label = paste0("p = ", formatC(p, format = "e", digits = 2)),
    y.position = ave(p, Pathway, FUN = function(x)
      seq(
        from = max(df_epi$Score, na.rm = TRUE) * 1.05,
        length.out = length(x),
        by = 0.1
      )
    )
  )

p2 <- ggplot(df_epi, aes(condition, Score, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 4) +
  stat_pvalue_manual(
    stat.test,
    label = "p.label",
    size = 7,
    tip.length = 0.01,
    bracket.nudge.y = 0.15
  ) +
  theme_minimal(base_size = 20) +
  labs(x = NULL, y = "GSVA Score") +
  theme(
    axis.text.x = element_text(size = 32, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 36, face = "bold"),
    legend.position = "none",
    strip.text = element_text(size = 30, face = "bold", lineheight = 0.9),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

ggsave(
  file.path(PLOT_DIR, "GSVA_epigenetic_pathways_significant_only.pdf"),
  p2,
  width = 32,
  height = 28
)

cor_epi <- cor(t(tf_activities), t(gsva_epi), method = "spearman")
cor_epi_df <- as.data.frame(cor_epi) %>% rownames_to_column("TF")
write_csv(cor_epi_df, file.path(OUT_DIR, "DoRothEA_VIPER_TF_vs_EpigeneticPathways_cor.csv"))

# ---------------------- Completion ----------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("[INFO] DoRothEA + VIPER full analysis complete.")
message("Results: ", normalizePath(OUT_DIR))
message("Plots:   ", normalizePath(PLOT_DIR))
