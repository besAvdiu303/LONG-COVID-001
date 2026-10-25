# LONG-COVID-001

Analysis code for my master's thesis *"Systematic Characterization of the Transcriptome of COVID-19 and Long COVID Patients Compared to Healthy Individuals"* (MSc Bio Data Science, FH Wiener Neustadt, 2025). The work was done at the Institute of Analytical Chemistry, University of Vienna, supervised by Dr. Thomas Mohr.

The question was whether long COVID shows up as its own stable state in the blood transcriptome, or whether it looks like leftover inflammation from the acute infection. I combined public bulk RNA-seq data of PBMCs from 10 GEO studies and compared five groups.

## Data

515 PBMC samples after filtering. I started from 28 GEO/SRA projects and kept only PBMC samples from people without known pre-existing conditions and with clear clinical metadata.

| Condition | n | GEO series |
|---|---|---|
| healthy | 46 | GSE179627, GSE206263, GSE224615 |
| exposedNoSymptoms | 79 | GSE162562, GSE206263 |
| acuteCOVID | 306 | GSE162562, GSE174621, GSE178824, GSE179627, GSE190680, GSE201530, GSE205244, GSE206263 |
| recovered | 61 | GSE179627, GSE190747 |
| longCOVID | 23 | GSE224615 |

`data/metadata/metadata_harmonized.csv` has one row per sample: study, GSM and SRR accession, age, gender, harmonized condition and the original label from GEO. The raw FASTQ files are not in this repository. They can be downloaded from SRA with the SRR accessions in the `sample_id` column (I used `fasterq-dump` from the SRA Toolkit).

## Pipeline

The scripts are in `scripts/` and run in this order:

1. `fastqc_multiqc.R` runs FastQC and MultiQC on the FASTQ files.
2. `Kallisto_alignment.R` quantifies each sample with kallisto against the Ensembl GRCh38 (release 107) cDNA index.
3. `generate_countMatrix.R` summarizes transcripts to genes with tximport and writes the gene count matrix. It also checks sequencing depth and makes a first PCA.
4. `DESeq2_analysis.R` fits DESeq2 with `design = ~ gse_number + condition`, so the study is a covariate. Genes with fewer than 10 counts in total are dropped. It writes all pairwise contrasts, the VST matrix, and a batch-corrected VST matrix (limma `removeBatchEffect` on the GSE number). The corrected matrix is the input for steps 6 to 9.
5. `DESeq2_visualization.R` makes volcano plots, MA plots and PCA before and after batch correction.
6. `GSVA_transformation.R` computes GSVA scores for WikiPathways gene sets (from MSigDB via msigdbr) and tests them between groups with limma.
7. `WGCNA_modules.R` builds a signed co-expression network (bicor) from the 8,000 most variable genes and correlates module eigengenes with condition, age and gender. `WGCNA_module_enrichment.R` then runs WikiPathways enrichment for modules with |r| ≥ 0.25 and FDR < 0.05.
8. `TF_and_epigenetics_analysis.R` estimates transcription factor activity with DoRothEA regulons (confidence A to C) and VIPER, compares groups, and correlates TF activity with the GSVA scores of epigenetic pathways. It needs the GSVA output from step 6.
9. `lioness_network_analysis.R` builds single-sample networks with LIONESS for selected WikiPathways and tests edge weights between groups with limma. The pathway IDs are set at the top of the script (`pathways_of_interest`).

Results go to `results/<step>/` and plots to `plots/<step>/`.

## Running it

I used R 4.5.1 on a Linux server. Main packages: DESeq2 1.48.2, limma 3.64.3, tximport 1.36.1, GSVA 2.2.0, msigdbr 25.1.1, WGCNA 1.73, clusterProfiler 4.16.0, dorothea, viper 1.32.0, lionessR 1.22.0, graphite 1.54.0, igraph, tidygraph, ggraph, EnhancedVolcano 1.26.0, org.Hs.eg.db 3.21.0, here, tidyverse. FastQC 0.12.1, MultiQC 1.18 and kallisto 0.50.1 have to be installed separately.

Run the scripts from inside `scripts/`, because some use paths like `../results` and others use `here::here()`. `Kallisto_alignment.R` and `generate_countMatrix.R` still contain the server paths I worked with (FASTQ folder, kallisto index, GTF file). Change those to your own paths first.

## What came out

In short, long COVID samples did not look like acute COVID and also not like recovered samples. The strongest signals were in chromatin regulation (for example EHMT2/G9a, PHF10, MLLT6 up, ARID2 down), in serotonin signalling, and in mitochondrial and NAD+ metabolism. TF activity differed for FOXP1, TP53 and RELA, among others. Many of these patterns are also described for ME/CFS. Details and figures are in the thesis.

## Limitations

- All 23 long COVID samples come from one study (GSE224615). Study and condition are therefore partly mixed up for this group, and the batch correction cannot fully separate them.
- After `removeBatchEffect` the samples still cluster partly by study in the PCA. The correction helps, but it does not remove study effects completely.
- The groups differ in age. The recovered group is much older (mean 73 years) than the others.
- The results are associations from bulk RNA-seq. They would need to be checked in a new cohort, ideally with single-cell or epigenetic data.

## Citation

Avdiu, B. (2025). *Systematic Characterization of the Transcriptome of COVID-19 and Long COVID Patients Compared to Healthy Individuals.* Master's thesis, FH Wiener Neustadt / University of Vienna.
