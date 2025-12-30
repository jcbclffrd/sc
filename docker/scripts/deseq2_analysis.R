#!/usr/bin/env Rscript
# deseq2_analysis.R - Differential TE analysis with DESeq2

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

# Get environment variables
count_file <- Sys.getenv("COUNT_MATRIX", "/data/te_counts_matrix.tsv")
metadata_file <- Sys.getenv("METADATA", "/metadata/samples.csv")
output_dir <- Sys.getenv("OUTPUT_DIR", "/output")
condition_col <- Sys.getenv("CONDITION_COL", "condition")
control_name <- Sys.getenv("CONTROL", "Control")
treatment_name <- Sys.getenv("TREATMENT", "AD")

cat("==========================================\n")
cat("DESeq2 Differential Analysis\n")
cat("==========================================\n")
cat("Count matrix:", count_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Condition column:", condition_col, "\n")
cat("Comparison:", treatment_name, "vs", control_name, "\n")
cat("==========================================\n\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read count matrix
counts <- read.delim(count_file, header = TRUE, row.names = 1)
cat("Loaded counts:", nrow(counts), "TEs ×", ncol(counts), "samples\n")

# Read metadata
metadata <- read.csv(metadata_file, header = TRUE)
rownames(metadata) <- metadata$sample_id

# Match samples
common_samples <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_samples]
metadata <- metadata[common_samples, ]

cat("Matched samples:", length(common_samples), "\n\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = as.formula(paste0("~ ", condition_col))
)

# Set reference level
dds[[condition_col]] <- relevel(dds[[condition_col]], ref = control_name)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
cat("After filtering:", nrow(dds), "TEs retained\n\n")

# Run DESeq2
cat("Running DESeq2...\n")
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c(condition_col, treatment_name, control_name))
res_ordered <- res[order(res$padj), ]

# Save results
write.csv(as.data.frame(res_ordered), 
          file = file.path(output_dir, "deseq2_results.csv"))

# Summary
cat("\nResults Summary:\n")
cat("Total TEs:", nrow(res), "\n")
cat("Upregulated (padj < 0.05, log2FC > 1):", 
    sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Downregulated (padj < 0.05, log2FC < -1):", 
    sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE), "\n")

# Volcano plot
pdf(file.path(output_dir, "volcano_plot.pdf"), width = 8, height = 6)
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, 
                              "Significant", "Not significant")
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_bw() +
  labs(title = paste(treatment_name, "vs", control_name),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
dev.off()

# MA plot
pdf(file.path(output_dir, "ma_plot.pdf"), width = 8, height = 6)
plotMA(res, main = paste(treatment_name, "vs", control_name))
dev.off()

# Heatmap of top 50 TEs
vsd <- vst(dds, blind = FALSE)
top_tes <- head(rownames(res_ordered), 50)
pdf(file.path(output_dir, "heatmap_top50.pdf"), width = 10, height = 12)
pheatmap(assay(vsd)[top_tes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = metadata[, condition_col, drop = FALSE],
         main = "Top 50 Differential TEs")
dev.off()

# PCA plot
pdf(file.path(output_dir, "pca_plot.pdf"), width = 8, height = 6)
plotPCA(vsd, intgroup = condition_col) +
  theme_bw() +
  ggtitle("PCA - All Samples")
dev.off()

cat("\n✓ Analysis complete!\n")
cat("Results saved to:", output_dir, "\n")
