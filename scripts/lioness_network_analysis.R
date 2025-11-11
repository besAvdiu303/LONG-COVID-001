#!/usr/bin/env Rscript
# lioness_network_analysis.R
#
# Description:
#   This script applies the LIONESS (Linear Interpolation to Obtain Network Estimates 
#   for Single Samples) framework to infer sample-specific gene co-expression networks 
#   for a set of biologically relevant WikiPathways. 
#
#   The workflow performs the following:
#     1. Loads VST-normalized, limma-corrected expression data and associated metadata.
#     2. Selects user-defined WikiPathways of interest via the graphite package.
#     3. Converts pathway graphs to gene symbol space and extracts intersecting genes 
#        from the expression matrix.
#     4. Computes pairwise gene–gene correlations using WGCNA’s bicor and applies LIONESS 
#        to estimate single-sample edge weights.
#     5. Performs differential edge-weight analysis between biological conditions 
#        (acuteCOVID, longCOVID, recovered, healthy) using limma linear modeling.
#     6. Identifies significantly rewired edges (P < 0.05) and visualizes:
#           - Edge-wise differential networks using tidygraph + ggraph.
#           - Volcano plots of edge-weight changes via EnhancedVolcano + ggrepel.
#     7. Combines comparisons (2×2 grid) into composite network and volcano plots 
#        for each pathway.
#     8. Exports per-pathway CSVs of significant edges and all plots to 
#        results/lioness_networks and plots/lioness_networks.
#



# --- Load packages -----------------------------------------------------------
library(graphite)
library(tidygraph)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(WGCNA)
library(lionessR)
library(igraph)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(rio)
library(here)
library(org.Hs.eg.db)
library(patchwork)
library(ggraph)
library(ggplot2)
library(ggrepel)

# --- Define directories ------------------------------------------------------
datadir <- here::here("data")
resultsdir <- here::here("results")
res_lioness <- here::here(resultsdir, "lioness_networks")
plotdir <- here::here("plots", "lioness_networks")

dir.create(res_lioness, showWarnings = FALSE, recursive = TRUE)
dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)

matrix_file <- here::here(resultsdir, "diff_exp_analysis", "VST_data_limma_corrected.RDS")
metadata_file <- here::here(datadir, "metadata", "metadata_harmonized.csv")

# --- Load expression and metadata -------------------------------------------
counts_data <- rio::import(matrix_file)
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
metadata <- metadata %>% mutate(across(everything(), \(x) trimws(as.character(x))))
cat("Matrix dimensions:", dim(counts_data), "\n")
cat("Metadata loaded with", nrow(metadata), "rows and", ncol(metadata), "columns.\n")

# Match by sample_id
common_samples <- intersect(colnames(counts_data), metadata$sample_id)
metadata <- metadata %>% filter(sample_id %in% common_samples)
counts_data <- counts_data[, common_samples]
cat("Common samples:", length(common_samples), "\n")

# --- Pathways of interest ----------------------------------------------------
# pathways_of_interest <- c("WP3678", "WP5117", "WP4946", "WP3996", "WP5483", "WP2369", "WP2916", "WP3630")
pathways_of_interest <- c("WP2369", "WP2916", "WP1455")


all_pathways <- graphite::pathways("hsapiens", "wikipathways")
selected_paths <- Filter(function(p) graphite::pathwayId(p) %in% pathways_of_interest, as.list(all_pathways))
cat("Retrieved", length(selected_paths), "WikiPathways of interest.\n")

# --- Gene ID mapping ---------------------------------------------------------
ensembl_ids <- rownames(counts_data)
symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL", keytype = "ENSEMBL",
  multiVals = "first"
)
valid <- !is.na(symbol_map) & symbol_map != ""
counts_matrix_symbol <- counts_data[valid, ]
rownames(counts_matrix_symbol) <- symbol_map[valid]
cat("Mapped genes:", length(unique(rownames(counts_matrix_symbol))), "\n")

# --- Enable multithreading ---------------------------------------------------
WGCNA::enableWGCNAThreads(nThreads = 10)

# --- Comparisons to perform --------------------------------------------------
comparisons <- list(
  acuteCOVID_vs_healthy   = c("acuteCOVID", "healthy"),
  recovered_vs_healthy    = c("recovered", "healthy"),
  acuteCOVID_vs_recovered = c("acuteCOVID", "recovered"),
  longCOVID_vs_recovered  = c("longCOVID", "recovered")
)

# --- Fisher-Z transform ------------------------------------------------------
fisher_z <- function(x) {
  x[x > 0.9999] <- 0.9999
  x[x < -0.9999] <- -0.9999
  0.5 * log((1 + x) / (1 - x))
}

# --- Main loop ---------------------------------------------------------------
for (p in selected_paths) {
  pw_id <- graphite::pathwayId(p)
  pw_name <- graphite::pathwayTitle(p)

  cat("\n============================\nAnalyzing Pathway:", pw_name, "\n============================\n")

  # --- Prepare pathway ------------------------------------------------------
  p_sym <- graphite::convertIdentifiers(p, to = "SYMBOL")
  g <- tryCatch(graph_from_graphnel(graphite::pathwayGraph(p_sym)), error = function(e) NULL)
  if (is.null(g)) next
  POI <- as_tbl_graph(g)

  pathway_genes <- POI %>%
    activate(nodes) %>%
    as.data.frame() %>%
    pull(name) %>%
    str_replace("^SYMBOL:", "") %>%
    unique()

  selected_genes <- intersect(pathway_genes, rownames(counts_matrix_symbol))
  cat("Selected genes in pathway:", length(selected_genes), "\n")
  if (length(selected_genes) < 5) {
    warning(paste("Too few overlapping genes — skipping pathway", pw_name))
    next
  }

  counts_matrix_vst_selected <- counts_matrix_symbol[selected_genes, ]

  # --- LIONESS single-sample correlation ------------------------------------
  funWGCNA <- function(x) WGCNA::cor(t(x), nThreads = 10)
  lioness_result <- lionessR::lioness(counts_matrix_vst_selected, f = funWGCNA)

  lioness_mat <- fisher_z(SummarizedExperiment::assay(lioness_result))
  lioness_mat[!is.finite(lioness_mat)] <- NA
  lioness_mat <- lioness_mat[rowSums(is.na(lioness_mat)) == 0, , drop = FALSE]
  if (nrow(lioness_mat) == 0) {
    warning(paste("All correlations invalid — skipping pathway", pw_name))
    next
  }

  # --- Extract edges safely -------------------------------------------------
  edge_mat <- as.data.frame(as_edgelist(g, names = TRUE))
  if (ncol(edge_mat) == 2) colnames(edge_mat) <- c("reg", "tar")
  edge_list_df <- edge_mat %>%
    mutate(
      reg = str_replace(reg, "^SYMBOL:", ""),
      tar = str_replace(tar, "^SYMBOL:", ""),
      edgename = paste(reg, tar, sep = "_")
    )

  selected_edges <- intersect(rownames(lioness_mat), edge_list_df$edgename)
  cat("Selected edges in pathway:", length(selected_edges), "\n")
  if (length(selected_edges) < 5) {
    warning(paste("Too few valid edges — skipping pathway", pw_name))
    next
  }

  lioness_wp <- lioness_mat[selected_edges, , drop = FALSE]
  rownames(metadata) <- metadata$sample_id
  lioness_wp <- lioness_wp[, rownames(metadata), drop = FALSE]

  # --- Initialize plot lists -------------------------------------------------
  plots_list <- list()
  sig_results_list <- list()
  net_plots_list <- list()

  # --- Differential analysis per comparison ---------------------------------
  plot_volcano <- function(res, title_text) {
    p <- EnhancedVolcano(
      res,
      lab = res$label,
      x = "logFC",
      y = "P.Value",
      subtitle = NULL,
      caption = NULL,
      pCutoff = 0.05,
      FCcutoff = 0,
      boxedLabels = FALSE,
      drawConnectors = TRUE,
      widthConnectors = 0.6,
      colConnectors = "black",
      labSize = 0, # disable internal labels
      pointSize = 4,
      colAlpha = 0.9,
      title = title_text,
      xlab = "Edge Weight Difference",
      ylab = bquote(~ -Log[10] ~ adjusted ~ italic(P)),
      col = c("grey70", "forestgreen", "royalblue", "red2"),
      legendLabels = c("Log2FC", "Not significant", "P-value", "Significant"),
      legendPosition = "bottom",
      legendLabSize = 14,
      legendIconSize = 4.5,
      borderWidth = 0.8,
      borderColour = "black"
    ) +
      coord_cartesian(xlim = c(-7, 7), ylim = c(0, 20)) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "bottom",
        legend.justification = "center",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        plot.background = element_rect(fill = NA, color = NA),
        panel.background = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = NA, color = NA)
      )
    # Add repelled labels AFTER EnhancedVolcano
    p <- p + geom_text_repel(
      data = res[res$label != "", ],
      aes(x = logFC, y = labelY, label = label),
      size = 4,
      max.overlaps = Inf,
      nudge_y = 0.8,
      box.padding = 0.4,
      point.padding = 0.3,
      force = 3,
      max.iter = 5000,
      segment.size = 0.3,
      segment.alpha = 0.8
    )
    return(p)
  }


  for (cmp_name in names(comparisons)) {
    grp <- comparisons[[cmp_name]]
    sub_meta <- metadata %>% filter(condition %in% grp)
    if (length(unique(sub_meta$condition)) < 2) next

    sub_data <- lioness_wp[, sub_meta$sample_id, drop = FALSE]
    group <- factor(sub_meta$condition, levels = grp)
    design <- model.matrix(~ 0 + group)
    colnames(design) <- grp
    contrast <- makeContrasts(contrasts = paste0(grp[1], " - ", grp[2]), levels = design)
    fit <- lmFit(sub_data, design) |>
      contrasts.fit(contrast) |>
      eBayes()
    res <- topTable(fit, coef = 1, number = Inf)

    # ---- Store significant results for later combination ----
    res_tbl <- res %>%
      tibble::rownames_to_column("edgename") %>%
      tidyr::separate(edgename, into = c("gene1", "gene2"), sep = "_", remove = FALSE) %>%
      dplyr::filter(P.Value < 0.05) %>%
      dplyr::mutate(comparison = cmp_name)

    sig_results_list[[cmp_name]] <- res_tbl

    message(pw_id, " | ", cmp_name, " | Significant edges: ", nrow(res_tbl))


    # --- Network with tidygraph --------------------------------------------
    res_tbl <- res %>% tibble::rownames_to_column("edgename")
    edges_tbl <- edge_list_df %>%
      left_join(res_tbl, by = "edgename") %>%
      mutate(
        weight = ifelse(is.na(logFC), 0, logFC),
        dir = ifelse(weight >= 0, "up", "down"),
        sig = P.Value < 0.05,
        weight_abs = abs(weight),
        weight_pos = ifelse(!is.finite(weight_abs) | weight_abs <= 0, 0.001, weight_abs)
      )

    net_tidy <- tbl_graph(
      nodes = tibble(name = unique(c(edges_tbl$reg, edges_tbl$tar))),
      edges = edges_tbl,
      directed = FALSE
    ) %>%
      activate(edges) %>%
      mutate(weight = weight_pos) %>%
      activate(nodes) %>%
      mutate(
        degree = centrality_degree(mode = "all"),
        betweenness = centrality_betweenness()
      )

    # --- Safe layout creation ----------------------------------------------
    layout_safe <- tryCatch(
      create_layout(net_tidy, layout = "fr", weights = net_tidy %E>% pull(weight)),
      error = function(e) {
        message("Layout weights invalid, using unweighted layout instead.")
        create_layout(net_tidy, layout = "fr")
      }
    )
    layout_df <- create_layout(net_tidy, layout = "fr")
    # --- Build network plot -------------------------------------------------
    p_net <- ggraph(layout_df) +
      geom_edge_link(aes(color = dir, width = weight_pos, alpha = sig)) +
      geom_node_point(aes(size = degree), color = "grey20") +
      ggrepel::geom_text_repel(
        data = layout_df,
        aes(x = x, y = y, label = name),
        size = 5,
        max.overlaps = Inf,
        box.padding = 0.6,
        point.padding = 0.6,
        force = 2,
        max.iter = 5000,
        segment.size = 0.3,
        segment.alpha = 0.6
      ) +
      scale_edge_color_manual(values = c("down" = "blue", "up" = "red")) +
      scale_edge_width(range = c(0.3, 2)) +
      scale_edge_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3), guide = "none") +
      theme_void() +
      ggtitle(gsub("_vs_", " vs ", cmp_name))


    net_plots_list[[cmp_name]] <- p_net

    # --- Volcano plot -------------------------------------------------------
    title_text <- gsub("_vs_", " vs ", cmp_name)

    # --- Volcano prep: choose labels & repel positions ----
    res$neglogP <- -log10(res$P.Value)

    # pick top 7 to label
    top_labels <- rownames(res)[order(res$neglogP, decreasing = TRUE)][1:7]
    res$label <- ifelse(rownames(res) %in% top_labels, rownames(res), "")

    # replace underscores with arrow
    res$label <- gsub("_", " \u279C ", res$label)

    # visual cap only — do NOT alter P.Value
    res$labelY <- pmin(res$neglogP, 15)


    plots_list[[cmp_name]] <- plot_volcano(res, gsub("_vs_", " vs ", cmp_name))
  }

  # --- After finishing all comparisons for this pathway ---
  sig_combined <- dplyr::bind_rows(sig_results_list)

  if (nrow(sig_combined) > 0) {
    out_file <- file.path(res_lioness, paste0("LIONESS_", pw_id, "_SIGNIFICANT_ALL_CONTRASTS.csv"))
    rio::export(sig_combined, out_file)
    message("Combined significant edges saved: ", out_file)
  } else {
    message("No significant edges for pathway ", pw_id, " — no file written.")
  }

  # --- Combined 2x2 Network plot -------------------------------------------
  want_keys <- c(
    "acuteCOVID_vs_healthy", "recovered_vs_healthy",
    "acuteCOVID_vs_recovered", "longCOVID_vs_recovered"
  )

  get_or_blank <- function(list_obj, key) {
    if (!is.null(list_obj[[key]])) {
      list_obj[[key]] +
        theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none" # hide individual legends for combined layout
        )
    } else {
      ggplot() +
        theme_void() +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
    }
  }

  # Combine 2×2 grid of network plots with spacing and larger labels
  combo_net <- ((get_or_blank(net_plots_list, "acuteCOVID_vs_healthy") |
    get_or_blank(net_plots_list, "acuteCOVID_vs_recovered")) /
    plot_spacer() / # vertical gap between rows
    (get_or_blank(net_plots_list, "recovered_vs_healthy") |
      get_or_blank(net_plots_list, "longCOVID_vs_recovered"))) +
    plot_layout(guides = "collect", heights = c(1, 0.05, 1)) +
    plot_annotation(
      title = paste0(pw_name, " (", pw_id, ")"),
      tag_levels = "A", # subplot letters A–D
      tag_suffix = ")"
    ) &
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.tag = element_text(size = 24, face = "bold"), # A–D labels +1 pt
      legend.position = "bottom", # global legend below
      legend.direction = "horizontal",
      legend.spacing.x = unit(1.2, "cm"), # space between legends
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      plot.margin = margin(10, 10, 30, 10)
    ) &
    guides(
      edge_width = "none",
      edge_alpha = "none",
      edge_color = guide_legend(title = "Regulation"), 
      size = guide_legend(title = "Degree") 
    )

  # Save the combined figure
  ggsave(
    filename = file.path(plotdir, paste0("Network_combined_", pw_id, ".pdf")),
    plot = combo_net, width = 14, height = 14, device = cairo_pdf
  )


  # --- Combined 2x2 Volcano plot -------------------------------------------
  combo_volcano <- ((get_or_blank(plots_list, "acuteCOVID_vs_healthy") |
    get_or_blank(plots_list, "acuteCOVID_vs_recovered")) /
    (get_or_blank(plots_list, "recovered_vs_healthy") |
      get_or_blank(plots_list, "longCOVID_vs_recovered"))) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(pw_name, " (", pw_id, ")"),
      tag_levels = "A",
      tag_suffix = ")"
    ) &
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.title.position = "plot",

      # Bigger & bold subplot tag letters (A,B,C,D)
      plot.tag = element_text(size = 28, face = "bold"),
      plot.tag.position = c(0.02, 0.98),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 18),
      plot.background = element_blank()
    )

  ggsave(
    filename = file.path(plotdir, paste0("Volcano_combined_", pw_id, ".pdf")),
    plot = combo_volcano, width = 14, height = 10, device = cairo_pdf
  )

  cat("Completed:", pw_name, "\n")
}
cat("\nAll selected pathways analyzed and combined volcano/network plots exported.\n")
