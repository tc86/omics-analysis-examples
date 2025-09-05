#!/usr/bin/env Rscript

# ------------------------------- #
#   DESeq2 Differential Analysis  #
#   Author: Tamrin Chowdhury      
#   Date: 07/13/2019
#   Usage:
#   Rscript deseq2_run.R \
#     --counts counts.tsv \
#     --meta meta.tsv \
#     --condition_col condition \
#     --ref Telomerase \
#     --outdir results_deseq2
# ------------------------------- #

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(ggplot2)
})

# 1) Read metadata (first column = sample IDs)
meta <- read_tsv("meta.tsv", col_types = cols())
rownames(meta) <- meta[[1]]
meta <- meta[, -1, drop = FALSE]  # keep only metadata columns

# 2) Read counts (first column = gene IDs)
cts <- read_tsv("counts.tsv", col_types = cols())

# If a column named 'Gene' exists, use it as rownames; otherwise use the first column
if ("Gene" %in% colnames(cts)) {
  rownames(cts) <- make.names(cts$Gene, unique = TRUE)
  cts <- cts %>% select(-Gene)
} else {
  rownames(cts) <- make.names(pull(cts, 1), unique = TRUE)
  cts <- cts[, -1, drop = FALSE]
}

# 3) Align samples between counts and metadata
common <- intersect(colnames(cts), rownames(meta))
if (length(common) < 2) stop("Need at least 2 overlapping samples between counts and meta.")
cts  <- as.matrix(cts[, common, drop = FALSE])
meta <- meta[common, , drop = FALSE]

# --- 4) Build DESeq2 object ---
# Ensure 'condition' exists and set reference to 'Telomerase'
if (!"condition" %in% colnames(meta)) stop("Metadata must contain a 'condition' column.")
meta$condition <- factor(meta$condition)
if (!("Telomerase" %in% levels(meta$condition))) {
  stop("Reference level 'Telomerase' not found in 'condition'.")
}
meta$condition <- relevel(meta$condition, ref = "Telomerase")

dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData   = meta,
                              design    = ~ condition)

# 5) Filter low-count genes (same spirit as your rowSums >= 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 6) Run DESeq2
dds <- DESeq(dds)

# Optional: see available result names
print(resultsNames(dds))

# 7) Extract results for ALT vs Telomerase and Neither vs Telomerase ---
AvsT <- results(dds, contrast = c("condition", "ALT", "Telomerase"))
NvsT <- results(dds, contrast = c("condition", "Neither", "Telomerase"))

# Order by adjusted p-value, then raw p-value (safer than p-value only)
AvsTordered <- as.data.frame(AvsT) %>%
  rownames_to_column("gene") %>%
  arrange(padj, pvalue)

NvsTordered <- as.data.frame(NvsT) %>%
  rownames_to_column("gene") %>%
  arrange(padj, pvalue)

# --- 8) Write results ---
write.table(AvsTordered, file = "ALT_vs_Telomerase.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(NvsTordered, file = "Neither_vs_Telomerase.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --- 9) Save normalized counts (handy for QC/plots) ---
norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene")
write.table(norm_counts, file = "normalized_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# --- 10) Quick MA plot (optional) ---
png("MA_plot.png", width = 900, height = 650, res = 120)
plotMA(dds, ylim = c(-4, 4))
dev.off()

message("Done. Files written: ALT_vs_Telomerase.tsv, Neither_vs_Telomerase.tsv, normalized_counts.tsv, MA_plot.png")
