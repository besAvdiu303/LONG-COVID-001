#!/usr/bin/env Rscript
# WGCNA_module_enrichment_wikipathways.R
#
# Description:
# Functional enrichment for WGCNA modules using WikiPathways
# (via msigdbr + clusterProfiler).
#
# Inputs:
#   - results/wgcna/WGCNA_gene_module_assignment.csv
#   - results/wgcna/WGCNA_ME_trait_correlations.csv
#   - results/wgcna/WGCNA_ME_trait_pvalues.csv
#   - results/wgcna/WGCNA_input_matrix.rds
#
# Outputs:
#   - results/wgcna/enrichment_wikipathways/*.csv
#   - plots/wgcna/enrichment_wikipathways/*.pdf
#   - results/wgcna/enrichment_wikipathways/enrichment_WikiPathways_allModules_top.csv
#
# Notes:
#   - Modules selected by |r| ≥ 0.25 and FDR < 0.05
#   - Automatically detects msigdbr column names (≥ v10)
#   - Pathway names cleaned ("WP"/IDs removed)

set.seed(123)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(here)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
})

# ---------------------- Paths and I/O ----------------------
here::i_am("scripts/WGCNA_module_enrichment.R")
message("here() root: ", here::here())

WGCNA_DIR <- here::here("results", "wgcna")
ASSIGN_CSV <- file.path(WGCNA_DIR, "WGCNA_gene_module_assignment.csv")
COR_CSV <- file.path(WGCNA_DIR, "WGCNA_ME_trait_correlations_with_covariates.csv")
PVAL_CSV <- file.path(WGCNA_DIR, "WGCNA_ME_trait_pvalues_with_covariates.csv")
UNIV_RDS <- file.path(WGCNA_DIR, "WGCNA_input_matrix.rds")

ENR_DIR <- file.path(WGCNA_DIR, "enrichment_wikipathways")
PLOTS_DIR <- here::here("plots", "wgcna", "enrichment_wikipathways")
dir.create(ENR_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(file.exists(ASSIGN_CSV, COR_CSV, PVAL_CSV, UNIV_RDS)))

# ---------------------- Load data ----------------------
assign <- read_csv(ASSIGN_CSV, show_col_types = FALSE)
corME <- read_csv(COR_CSV, show_col_types = FALSE)
pME <- read_csv(PVAL_CSV, show_col_types = FALSE)
expr_univ <- readRDS(UNIV_RDS)
univ_ens <- sub("\\..*$", "", rownames(expr_univ))

# ---------------------- Universe mapping ----------------------
ens_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(univ_ens),
    keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID")
) |> distinct(ENSEMBL, .keep_all = TRUE)

univ_tbl <- tibble(ENSEMBL = univ_ens) |> left_join(ens_map, by = "ENSEMBL")
univ_entrez <- unique(na.omit(univ_tbl$ENTREZID))
univ_symbol <- unique(na.omit(univ_tbl$SYMBOL))
message(sprintf(
    "Universe sizes — ENSEMBL=%d, ENTREZ=%d, SYMBOL=%d",
    length(univ_ens), length(univ_entrez), length(univ_symbol)
))

# ---------------------- Select modules by trait association ----------------------
cor_long <- pivot_longer(corME, -ME, names_to = "trait", values_to = "cor")
p_long <- pivot_longer(pME, -ME, names_to = "trait", values_to = "pval")
stats_long <- left_join(cor_long, p_long, by = c("ME", "trait")) |>
    mutate(FDR = p.adjust(pval, "BH"))

abs_cor_thr <- 0.25
fdr_thr <- 0.05
hit_long <- stats_long |> filter(abs(cor) >= abs_cor_thr, FDR < fdr_thr)
hit_long$module <- sub("^ME", "", hit_long$ME)
mods <- sort(unique(hit_long$module[hit_long$module != "grey"]))
message(
    "Modules selected (|r|≥", abs_cor_thr, ", FDR<", fdr_thr, "): ",
    ifelse(length(mods) > 0, paste(mods, collapse = ", "), "<none>")
)

if (length(mods) == 0) {
    message("No modules meet selection criteria — exiting gracefully.")
    quit(save = "no")
}

# ---------------------- Load WikiPathways via msigdbr ----------------------
msig_wp <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
stopifnot(nrow(msig_wp) > 0)
message("WikiPathways entries: ", nrow(msig_wp))

cols <- names(msig_wp)
gene_col <- case_when(
    "entrez_gene" %in% cols ~ "entrez_gene",
    "entrez_gene_id" %in% cols ~ "entrez_gene_id",
    "gene_symbol" %in% cols ~ "gene_symbol",
    "human_gene_symbol" %in% cols ~ "human_gene_symbol",
    "ensembl_gene" %in% cols ~ "ensembl_gene",
    "ensembl_gene_id" %in% cols ~ "ensembl_gene_id",
    TRUE ~ NA_character_
)
stopifnot(!is.na(gene_col))

# ---------------------- Prepare TERM2GENE and TERM2NAME ----------------------
clean_pathway <- function(x) {
    x <- gsub("\\s*\\(WP\\d+\\)$", "", x)
    x <- gsub("^WP_?", "", x)
    trimws(x)
}

TERM2GENE <- msig_wp |>
    select(term = gs_name, gene = all_of(gene_col)) |>
    mutate(gene = sub("\\..*$", "", as.character(gene))) |>
    distinct() |>
    filter(!is.na(gene) & gene != "")

TERM2NAME <- msig_wp |>
    select(term = gs_name, name = gs_name) |>
    mutate(name = clean_pathway(name)) |>
    distinct()

# ---------------------- Choose best ID space ----------------------
id_spaces <- list(
    ENSEMBL = intersect(TERM2GENE$gene, univ_ens),
    SYMBOL  = intersect(TERM2GENE$gene, univ_symbol),
    ENTREZ  = intersect(TERM2GENE$gene, univ_entrez)
)
lens <- vapply(id_spaces, length, integer(1))
best <- names(which.max(lens))
universe_cp <- id_spaces[[best]]
message("Using ID space: ", best, " (universe size = ", length(universe_cp), ")")

# ---------------------- Helper functions ----------------------
assign$ENSEMBL <- sub("\\..*$", "", assign$ENSEMBL)
genes_for_module <- function(df, id_space = best) {
    if (id_space == "SYMBOL") {
        unique(na.omit(df$SYMBOL))
    } else if (id_space == "ENSEMBL") {
        unique(na.omit(df$ENSEMBL))
    } else {
        m <- AnnotationDbi::select(org.Hs.eg.db,
            keys = unique(df$ENSEMBL), keytype = "ENSEMBL", columns = "ENTREZID"
        ) |> distinct(ENSEMBL, .keep_all = TRUE)
        unique(na.omit(m$ENTREZID))
    }
}

calc_height <- function(n, base = 3, per = 0.45, min_h = 6) {
    max(min_h, base + per * n)
}

# ---------------------- Enrichment per module ----------------------
mod_lists <- split(assign, assign$module)
mod_lists <- mod_lists[names(mod_lists) %in% mods]

for (mod in names(mod_lists)) {
    df <- mod_lists[[mod]]
    genes <- intersect(genes_for_module(df, best), universe_cp)
    if (length(genes) < 5) {
        message("Module ", mod, ": too few genes after filtering.")
        next
    }

    ego <- enricher(
        gene = genes,
        universe = universe_cp,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME,
        pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 5, maxGSSize = 1000
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
        message("No enrichment hits for module ", mod)
        next
    }

    ego@result$Description <- clean_pathway(ego@result$Description)
    res_df <- as.data.frame(ego) |> mutate(module = mod, ID_SPACE = best)

    out_csv <- file.path(ENR_DIR, paste0("enrichment_WikiPathways_", mod, ".csv"))
    write_csv(res_df, out_csv)
    message("Saved: ", out_csv)

    # Save top N pathways for combined summary
    n_show <- min(20, nrow(res_df))
    topN_df <- res_df |>
        arrange(p.adjust) |>
        slice_head(n = n_show) |>
        mutate(module = mod, rank_within_module = row_number())

    if (!exists("all_top_list")) all_top_list <- list()
    all_top_list[[mod]] <- topN_df

    # Visualization
    h_dot <- calc_height(n_show)

    pdf(file.path(PLOTS_DIR, paste0("WikiPathways_", mod, "_dotplot.pdf")), width = 8, height = h_dot)
    p_dot <- dotplot(ego, showCategory = n_show) +
        ggtitle(paste0("Module ", mod)) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    print(p_dot)
    dev.off()

    assign(paste0("p_", mod), p_dot)

    ego2 <- enrichplot::pairwise_termsim(ego)
    pdf(file.path(PLOTS_DIR, paste0("WikiPathways_", mod, "_emapplot.pdf")), width = 9, height = 8)
    print(emapplot(ego2, showCategory = min(15, nrow(res_df))))
    dev.off()

    pdf(file.path(PLOTS_DIR, paste0("WikiPathways_", mod, "_cnetplot.pdf")), width = 9, height = 8)
    print(cnetplot(ego, showCategory = min(8, nrow(res_df)), circular = FALSE, node_label = "all"))
    dev.off()
}

# ---------------------- Combine all modules ----------------------
alls <- list.files(ENR_DIR, "^enrichment_WikiPathways_.*\\.csv$", full.names = TRUE)
if (length(alls) > 0) {
    all_df <- bind_rows(lapply(alls, read_csv, show_col_types = FALSE))
    write_csv(all_df, file.path(ENR_DIR, "enrichment_WikiPathways_allModules_combined.csv"))
    message("Saved all enrichment results to enrichment_WikiPathways_allModules_combined.csv")
}

if (exists("all_top_list") && length(all_top_list) > 0) {
    all_top_df <- bind_rows(all_top_list)
    write_csv(all_top_df, file.path(ENR_DIR, "enrichment_WikiPathways_allModules_top.csv"))
    message("Saved top-N (<=20) pathways per module to enrichment_WikiPathways_allModules_top.csv")
}

# ---------------------- Combine turquoise and yellow dotplots ----------------------
suppressPackageStartupMessages({
    library(patchwork)
})

if (exists("p_turquoise") && exists("p_yellow")) {
    p_combined <- (p_turquoise | p_yellow) +
        plot_annotation(tag_levels = "A", tag_suffix = ")") &
        theme(plot.tag = element_text(size = 20, face = "bold"))

    file_out <- file.path(PLOTS_DIR, "WikiPathways_combined_turquoise_yellow_dotplots.pdf")
    ggsave(file_out, plot = p_combined, width = 16, height = 14, device = cairo_pdf)
    message("Combined plot saved to: ", file_out)
} else {
    warning("Required plots (p_turquoise, p_yellow) not found — skipping combination.")
}

# ---------------------- Session info ----------------------
sink(file.path(ENR_DIR, "sessionInfo.txt"))
sessionInfo()
sink()

message("WikiPathways enrichment completed successfully.")
message("Results: ", normalizePath(ENR_DIR))
message("Plots:   ", normalizePath(PLOTS_DIR))
