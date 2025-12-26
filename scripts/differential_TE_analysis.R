#!/usr/bin/env Rscript

# Differential TE Expression Analysis: AD vs Control
# Using DESeq2 for robust differential expression analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

# Set working directory
setwd("~/sc/tetranscripts_bulk")

# Load count matrix
cat("Loading TE count matrix...\n")
te_counts <- read.table("combined_counts_matrix_TEs_only.tsv", 
                        header = TRUE, 
                        row.names = 1, 
                        sep = "\t",
                        check.names = FALSE)

cat(sprintf("Loaded %d TEs across %d samples\n", nrow(te_counts), ncol(te_counts)))

# Extract sample IDs and patient numbers from column names
sample_info <- data.frame(
  sample_id = colnames(te_counts),
  srr = sapply(strsplit(colnames(te_counts), "_"), `[`, 1),
  patient = sapply(strsplit(colnames(te_counts), "_"), `[`, 2),
  stringsAsFactors = FALSE
)

# Load patient metadata with diagnosis
metadata <- read.csv("../data/patient_multiomics_mapping.csv", stringsAsFactors = FALSE)

# Merge with diagnosis information
sample_info <- merge(sample_info, 
                     metadata[, c("patient_sample_id", "diagnosis", "age", "sex", "tangle_stage", "plaque_stage")],
                     by.x = "patient", 
                     by.y = "patient_sample_id",
                     all.x = TRUE)

# Reorder to match count matrix columns
sample_info <- sample_info[match(colnames(te_counts), sample_info$sample_id), ]
rownames(sample_info) <- sample_info$sample_id

# Ensure diagnosis is a factor with Control as reference
sample_info$diagnosis <- factor(sample_info$diagnosis, levels = c("Control", "AD"))

cat("\nSample distribution:\n")
print(table(sample_info$diagnosis))

# Filter out lowly expressed TEs (present in at least 25% of samples with >= 10 reads)
min_samples <- ceiling(ncol(te_counts) * 0.25)
keep <- rowSums(te_counts >= 10) >= min_samples
te_counts_filtered <- te_counts[keep, ]

cat(sprintf("\nFiltered to %d TEs (from %d) expressed in >= %d samples\n", 
            sum(keep), nrow(te_counts), min_samples))

# Create DESeq2 dataset
cat("\nCreating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = te_counts_filtered,
  colData = sample_info,
  design = ~ diagnosis
)

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get results
cat("Extracting results...\n")
res <- results(dds, contrast = c("diagnosis", "AD", "Control"), alpha = 0.05)

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Convert to dataframe and add TE classification
res_df <- as.data.frame(res_ordered)
res_df$TE <- rownames(res_df)

# Parse TE family and class from name (format: family:subfamily:class)
te_info <- strsplit(res_df$TE, ":")
res_df$TE_family <- sapply(te_info, function(x) if(length(x) >= 2) x[2] else x[1])
res_df$TE_class <- sapply(te_info, function(x) if(length(x) >= 3) x[3] else "Unknown")

# Add significance categories
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Upregulated in AD"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Downregulated in AD"
res_df$significance <- factor(res_df$significance, 
                              levels = c("Upregulated in AD", "Downregulated in AD", "Not Significant"))

# Save full results
write.csv(res_df, "differential_TE_results_AD_vs_Control.csv", row.names = FALSE)

cat("\n=============================================================\n")
cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("=============================================================\n")
cat(sprintf("Total TEs tested: %d\n", nrow(res_df)))
cat(sprintf("Significantly DE (padj < 0.05): %d\n", sum(res_df$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  - Upregulated in AD: %d\n", sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)))
cat(sprintf("  - Downregulated in AD: %d\n", sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)))
cat("=============================================================\n\n")

# Show top upregulated TEs in AD
cat("Top 20 TEs UPREGULATED in AD:\n")
top_up <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange > 0, ] %>%
  arrange(padj) %>%
  head(20)
if(nrow(top_up) > 0) {
  print(top_up[, c("TE", "baseMean", "log2FoldChange", "padj", "TE_class")])
} else {
  cat("  None\n")
}

cat("\nTop 20 TEs DOWNREGULATED in AD:\n")
top_down <- res_df[res_df$padj < 0.05 & res_df$log2FoldChange < 0, ] %>%
  arrange(padj) %>%
  head(20)
if(nrow(top_down) > 0) {
  print(top_down[, c("TE", "baseMean", "log2FoldChange", "padj", "TE_class")])
} else {
  cat("  None\n")
}

# Create visualizations
cat("\nGenerating visualizations...\n")

# 1. Volcano plot
pdf("volcano_plot_TE_AD_vs_Control.pdf", width = 10, height = 8)
p1 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated in AD" = "red", 
                                 "Downregulated in AD" = "blue", 
                                 "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(title = "Differential TE Expression: AD vs Control",
       x = "Log2 Fold Change (AD / Control)",
       y = "-Log10(Adjusted P-value)",
       color = "Significance") +
  theme(legend.position = "bottom")
print(p1)
dev.off()

# 2. MA plot
pdf("MA_plot_TE_AD_vs_Control.pdf", width = 10, height = 8)
p2 <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated in AD" = "red", 
                                 "Downregulated in AD" = "blue", 
                                 "Not Significant" = "grey")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  theme_bw() +
  labs(title = "MA Plot: TE Expression AD vs Control",
       x = "Log10(Mean Expression)",
       y = "Log2 Fold Change (AD / Control)",
       color = "Significance") +
  theme(legend.position = "bottom")
print(p2)
dev.off()

# 3. TE class enrichment
if(sum(res_df$padj < 0.05, na.rm = TRUE) > 0) {
  class_summary <- res_df %>%
    filter(padj < 0.05) %>%
    group_by(TE_class, significance) %>%
    summarise(count = n(), .groups = "drop")
  
  pdf("TE_class_enrichment_AD_vs_Control.pdf", width = 10, height = 6)
  p3 <- ggplot(class_summary, aes(x = reorder(TE_class, count), y = count, fill = significance)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Upregulated in AD" = "red", 
                                  "Downregulated in AD" = "blue")) +
    coord_flip() +
    theme_bw() +
    labs(title = "TE Class Distribution of Differentially Expressed TEs",
         x = "TE Class",
         y = "Number of DE TEs",
         fill = "Direction")
  print(p3)
  dev.off()
}

# 4. Heatmap of top DE TEs
if(sum(res_df$padj < 0.05, na.rm = TRUE) > 0) {
  top_tes <- res_df %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(50) %>%
    pull(TE)
  
  if(length(top_tes) > 0) {
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    top_counts <- norm_counts[top_tes, ]
    
    # Z-score transformation
    top_counts_scaled <- t(scale(t(log2(top_counts + 1))))
    
    # Annotation
    annotation_col <- data.frame(
      Diagnosis = sample_info$diagnosis,
      row.names = colnames(top_counts_scaled)
    )
    
    ann_colors <- list(
      Diagnosis = c("Control" = "#4DAF4A", "AD" = "#E41A1C")
    )
    
    pdf("heatmap_top50_DE_TEs.pdf", width = 12, height = 14)
    pheatmap(top_counts_scaled,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
             breaks = seq(-3, 3, length.out = 101),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = FALSE,
             main = "Top 50 Differentially Expressed TEs (AD vs Control)",
             fontsize_row = 6)
    dev.off()
  }
}

# 5. PCA plot
vsd <- vst(dds, blind = FALSE)
pdf("PCA_plot_TE_expression.pdf", width = 10, height = 8)
pcaData <- plotPCA(vsd, intgroup = "diagnosis", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p5 <- ggplot(pcaData, aes(PC1, PC2, color = diagnosis, label = name)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, vjust = 0, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Control" = "#4DAF4A", "AD" = "#E41A1C")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  labs(title = "PCA of TE Expression: AD vs Control") +
  theme(legend.position = "bottom")
print(p5)
dev.off()

cat("\n=============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================\n")
cat("Output files:\n")
cat("  - differential_TE_results_AD_vs_Control.csv\n")
cat("  - volcano_plot_TE_AD_vs_Control.pdf\n")
cat("  - MA_plot_TE_AD_vs_Control.pdf\n")
cat("  - TE_class_enrichment_AD_vs_Control.pdf\n")
cat("  - heatmap_top50_DE_TEs.pdf\n")
cat("  - PCA_plot_TE_expression.pdf\n")
cat("=============================================================\n")
