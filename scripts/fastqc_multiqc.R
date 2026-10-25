# ==============================================
# File: fastqc_multiqc.R
# Description:
#   Runs FastQC and MultiQC on FASTQ files that match
#   SRR numbers in the metadata file.
# ==============================================

# --- Load libraries ---
library(readxl)
library(dplyr)
library(stringr)

# --- Define paths (relative to project root) ---
metadata_path <- "../data/02-metadata/metadata_PBMC_subset2.xlsx"
fastq_dir <- "../data/01-raw/" # adjust if needed
qc_output <- "../results/fastqc/" # where results go

# --- Create output directory if missing ---
if (!dir.exists(qc_output)) dir.create(qc_output, recursive = TRUE)

# --- 1. Read metadata and extract SRR IDs ---
metadata <- read_excel(metadata_path)
metadata <- metadata %>%
    rename_with(~ tolower(gsub("[^A-Za-z0-9]", "_", .x)))

# try to identify the SRR column automatically
srr_col <- names(metadata)[grepl("srr", names(metadata), ignore.case = TRUE)]
if (length(srr_col) == 0) {
    stop("❌ No SRR column found in metadata. Please check the column names.")
}

srr_ids <- unique(metadata[[srr_col[1]]])
srr_ids <- srr_ids[!is.na(srr_ids)]
cat("✅ Found", length(srr_ids), "SRR IDs in metadata.\n")

# --- 2. List all FASTQ files and filter ---
all_fastqs <- list.files(fastq_dir, pattern = "\\.fastq(\\.gz)?$", full.names = TRUE)
matched_fastqs <- all_fastqs[
    sapply(srr_ids, function(id) any(grepl(id, all_fastqs)))
]

if (length(matched_fastqs) == 0) {
    stop("❌ No FASTQ files matching SRR IDs found.")
}

cat("✅ Found", length(matched_fastqs), "FASTQ files matching metadata SRR IDs.\n")

# --- 3. Run FastQC on matched files ---
# (Make sure fastqc is in your PATH)
cat("🚀 Running FastQC...\n")
fastqc_cmd <- paste(
    "fastqc",
    paste(shQuote(matched_fastqs), collapse = " "),
    "-o", shQuote(qc_output),
    "--quiet"
)
system(fastqc_cmd)

cat("✅ FastQC completed.\n")

# --- 4. Run MultiQC ---
cat("📊 Running MultiQC...\n")
multiqc_cmd <- paste("multiqc", shQuote(qc_output), "-o", shQuote(qc_output))
system(multiqc_cmd)

cat("🎉 MultiQC report generated in:", qc_output, "\n")
