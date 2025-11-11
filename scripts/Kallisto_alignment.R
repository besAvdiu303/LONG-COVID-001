#!/usr/bin/env Rscript
# ============================================================
# Kallisto_alignment.R
#
# Description:
#   Automates Kallisto quantification for all paired-end FASTQ
#   files in the specified directory. Pairs reads automatically
#   based on filename conventions and creates per-sample output
#   directories for each SRR/sample ID.
# ============================================================

# ---------- Load libraries ----------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(here)
})

# ---------- Define paths ----------
fastq_dir <- "/srv/storage/RNASeq/LONGCOVID/paired_end"
output_dir <- "../data/kallisto_alignment"
genome_idx <- "/srv/opt/genomes/hs/GRCh38_107/Homo_sapiens.GRCh38.cdna.all.idx"
kallisto_bin <- "/opt/kallisto/kallisto"

# ---------- Validate inputs ----------
if (!dir.exists(fastq_dir)) {
  stop("FASTQ directory not found: ", fastq_dir)
}
if (!file.exists(genome_idx)) {
  stop("Kallisto index file not found: ", genome_idx)
}
if (!file.exists(kallisto_bin)) {
  stop("Kallisto binary not found at: ", kallisto_bin)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---------- Identify and organize FASTQ pairs ----------
workflow <- data.frame(
  FILES = dir(
    fastq_dir,
    pattern = "_[12]\\.fastq\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
) %>%
  mutate(SAMPLES = basename(FILES)) %>%
  mutate(SAMPLES = str_replace(SAMPLES, "_[12]\\.fastq\\.gz$", "")) %>%
  group_by(SAMPLES) %>%
  summarise(
    FILES = paste(FILES, collapse = " "), # combine R1 and R2
    OUTPUT = unique(SAMPLES),
    .groups = "drop"
  ) %>%
  mutate(OUTPUT = here::here(output_dir, OUTPUT)) %>%
  mutate(
    command = paste(
      kallisto_bin, "quant -i", genome_idx,
      "-t 10 -o", shQuote(OUTPUT), FILES
    )
  )

# ---------- Create output directories ----------
walk(workflow$OUTPUT, ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# ---------- Execute Kallisto commands ----------
cat("Running Kallisto quantification for", nrow(workflow), "samples...\n")

walk2(workflow$SAMPLES, workflow$command, function(samp, cmd) {
  cat(" → Processing sample:", samp, "\n")
  system(cmd)
})

cat("\n✅ Kallisto alignment completed for all samples.\n")
cat("Results stored in:", normalizePath(output_dir), "\n")
