#!/usr/bin/env Rscript

# Differential Expression Analysis for Bulk RNA-seq TEs using DESeq2
# AD vs Control comparison

invisible({
  library(DESeq2)
  library(ggplot2)
})

cat("=================================================================\n")
cat("DIFFERENTIAL TE EXPRESSION ANALYSIS: AD vs CONTROL\n")
cat("=================================================================\n\n")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript differential_TE_analysis_bulk.R <count_matrix> <metadata> <output_dir>\n\n")
  cat("Arguments:\n")
  cat("  count_matrix  : Path to combined counts matrix (TSV format)\n")
  cat("  metadata      : Path to metadata CSV file\n")
  cat("  output_dir    : Directory to save results\n\n")
  quit(status = 1)
}

count_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

# Validate inputs
if (!file.exists(count_file)) {
  cat("Error: Count matrix file not found:", count_file, "\n")
  quit(status = 1)
}

if (!file.exists(metadata_file)) {
  cat("Error: Metadata file not found:", metadata_file, "\n")
  quit(status = 1)
}

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Loading count matrix...\n")
counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat(sprintf("  %d features x %d samples\n", nrow(counts), ncol(counts)))

cat("\nLoading metadata...\n")
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
cat(sprintf("  %d samples in metadata\n", nrow(metadata)))

# Filter metadata to included samples
if ("included_in_de_analysis" %in% colnames(metadata)) {
  metadata <- metadata[metadata$included_in_de_analysis == "Yes", ]
  cat(sprintf("  Using %d samples for DE analysis\n", nrow(metadata)))
}

# Match samples between counts and metadata
common_samples <- intersect(colnames(counts), metadata$srr_id)
cat(sprintf("  Matched %d samples between counts and metadata\n", length(common_samples)))

# Filter both to common samples
counts <- counts[, common_samples, drop = FALSE]
metadata <- metadata[metadata$srr_id %in% common_samples, ]
metadata <- metadata[match(common_samples, metadata$srr_id), ]
rownames(metadata) <- metadata$srr_id

# Check we have both AD and Control
if (!"diagnosis" %in% colnames(metadata)) {
  cat("\nError: 'diagnosis' column not found in metadata\n")
  quit(status = 1)
}

n_ad <- sum(metadata$diagnosis == "AD", na.rm = TRUE)
n_control <- sum(metadata$diagnosis == "Control", na.rm = TRUE)

if (n_ad == 0 || n_control == 0) {
  cat("\nError: Need both AD and Control samples for comparison\n")
  cat(sprintf("  AD samples: %d\n", n_ad))
  cat(sprintf("  Control samples: %d\n", n_control))
  quit(status = 1)
}

cat(sprintf("\nFinal sample count: %d\n", ncol(counts)))
cat(sprintf("  AD samples: %d\n", n_ad))
cat(sprintf("  Control samples: %d\n", n_control))

# ============================================================================
# 2. IDENTIFY GENES vs TEs
# ============================================================================

cat("\nIdentifying genes vs TEs...\n")

# TEtranscripts convention: genes start with "ENS", TEs do not
is_gene <- grepl("^ENS", rownames(counts))
is_te <- !is_gene

cat(sprintf("  Genes: %d\n", sum(is_gene)))
cat(sprintf("  TEs: %d\n", sum(is_te)))

# ============================================================================
# 3. ANNOTATE TEs
# ============================================================================

cat("\nAnnotating TEs...\n")

# TE naming convention in TEtranscripts:
# Format: TEfamily or TEfamily:TEsubfamily
# Examples: L1HS, AluY, LINE1:L1HS, SINE:Alu:AluY

annotate_te <- function(te_names) {
  te_df <- data.frame(
    TE = te_names,
    TE_family = NA_character_,
    TE_class = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(te_names)) {
    te_name <- te_names[i]
    
    # Parse the TE name
    if (grepl(":", te_name)) {
      # Format: class:family:subfamily or family:subfamily
      parts <- strsplit(te_name, ":")[[1]]
      if (length(parts) == 3) {
        te_df$TE_class[i] <- parts[1]
        te_df$TE_family[i] <- parts[2]
      } else if (length(parts) == 2) {
        te_df$TE_family[i] <- parts[1]
      }
    } else {
      # Just the subfamily name
      # Try to infer family from name
      if (grepl("^L1", te_name)) {
        te_df$TE_family[i] <- "LINE-1"
        te_df$TE_class[i] <- "LINE"
      } else if (grepl("^L2", te_name)) {
        te_df$TE_family[i] <- "LINE-2"
        te_df$TE_class[i] <- "LINE"
      } else if (grepl("^Alu", te_name)) {
        te_df$TE_family[i] <- "Alu"
        te_df$TE_class[i] <- "SINE"
      } else if (grepl("^SVA", te_name)) {
        te_df$TE_family[i] <- "SVA"
        te_df$TE_class[i] <- "Retrotransposon"
      } else if (grepl("^HERV", te_name) || grepl("^LTR", te_name)) {
        te_df$TE_family[i] <- "LTR"
        te_df$TE_class[i] <- "LTR"
      } else {
        te_df$TE_family[i] <- "Unknown"
        te_df$TE_class[i] <- "Unknown"
      }
    }
  }
  
  return(te_df)
}

te_names <- rownames(counts)[is_te]
te_annotations <- annotate_te(te_names)

cat(sprintf("  Annotated %d TEs\n", nrow(te_annotations)))

# ============================================================================
# 4. RUN DESEQ2 ON TEs
# ============================================================================

cat("\n=================================================================\n")
cat("RUNNING DESEQ2 ON TEs\n")
cat("=================================================================\n\n")

# Extract TE counts
te_counts <- counts[is_te, , drop = FALSE]
cat(sprintf("Running DESeq2 on %d TEs...\n", nrow(te_counts)))

# Ensure counts are integers
te_counts <- round(te_counts)

# Create DESeq2 dataset
metadata$diagnosis <- factor(metadata$diagnosis, levels = c("Control", "AD"))

dds_te <- DESeqDataSetFromMatrix(
  countData = te_counts,
  colData = metadata,
  design = ~ diagnosis
)

cat("  Filtering low-count TEs (minimum 10 counts total)...\n")
dds_te <- dds_te[rowSums(counts(dds_te)) >= 10, ]
cat(sprintf("  Retained %d TEs for analysis\n", nrow(dds_te)))

cat("  Running DESeq2...\n")
dds_te <- DESeq(dds_te, quiet = TRUE)

cat("  Extracting results...\n")
res_te <- results(dds_te, contrast = c("diagnosis", "AD", "Control"))
res_te_df <- as.data.frame(res_te)
res_te_df$TE <- rownames(res_te_df)

# Add TE annotations
res_te_df <- merge(res_te_df, te_annotations, by = "TE", all.x = TRUE)

# Add significance flag
res_te_df$significance <- ifelse(
  res_te_df$padj < 0.05 & !is.na(res_te_df$padj),
  ifelse(res_te_df$log2FoldChange > 0, "Up in AD", "Down in AD"),
  "Not significant"
)

# Sort by adjusted p-value
res_te_df <- res_te_df[order(res_te_df$pvalue, na.last = TRUE), ]

# Save results
output_file <- file.path(output_dir, "differential_TE_results_AD_vs_Control.csv")
cat(sprintf("\nSaving results to: %s\n", output_file))
write.csv(res_te_df, output_file, row.names = FALSE, quote = FALSE)

# ============================================================================
# 5. SUMMARY STATISTICS
# ============================================================================

cat("\n=================================================================\n")
cat("RESULTS SUMMARY\n")
cat("=================================================================\n\n")

n_tested <- sum(!is.na(res_te_df$padj))
n_sig_fdr <- sum(res_te_df$padj < 0.05, na.rm = TRUE)
n_sig_nom <- sum(res_te_df$pvalue < 0.05, na.rm = TRUE)
n_up <- sum(res_te_df$padj < 0.05 & res_te_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_te_df$padj < 0.05 & res_te_df$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("TEs tested: %d\n", n_tested))
cat(sprintf("Significant at FDR < 0.05: %d\n", n_sig_fdr))
cat(sprintf("  Upregulated in AD: %d\n", n_up))
cat(sprintf("  Downregulated in AD: %d\n", n_down))
cat(sprintf("Nominally significant (p < 0.05): %d\n\n", n_sig_nom))

if (n_sig_fdr > 0) {
  cat("Top 20 significant TEs (by adjusted p-value):\n")
  top_sig <- head(res_te_df[res_te_df$padj < 0.05 & !is.na(res_te_df$padj), ], 20)
  print(top_sig[, c("TE", "TE_family", "TE_class", "log2FoldChange", "padj")])
} else {
  cat("No TEs passed FDR threshold. Top 20 by nominal p-value:\n")
  top_nom <- head(res_te_df[!is.na(res_te_df$pvalue), ], 20)
  print(top_nom[, c("TE", "TE_family", "TE_class", "log2FoldChange", "pvalue", "padj")])
}

# ============================================================================
# 6. VISUALIZATION
# ============================================================================

cat("\n=================================================================\n")
cat("GENERATING VISUALIZATIONS\n")
cat("=================================================================\n\n")

# Volcano plot
cat("Creating volcano plot...\n")
volcano_file <- file.path(output_dir, "volcano_plot_TE_AD_vs_Control.pdf")

res_te_plot <- res_te_df[!is.na(res_te_df$pvalue), ]

p_volcano <- ggplot(res_te_plot, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Up in AD" = "#E74C3C",
    "Down in AD" = "#3498DB",
    "Not significant" = "gray60"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "TE Differential Expression: AD vs Control",
    x = "Log2 Fold Change (AD / Control)",
    y = "-Log10(p-value)",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(volcano_file, p_volcano, width = 10, height = 7)
cat(sprintf("  Saved: %s\n", volcano_file))

# MA plot
cat("Creating MA plot...\n")
ma_file <- file.path(output_dir, "MA_plot_TE_AD_vs_Control.pdf")

p_ma <- ggplot(res_te_plot, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Up in AD" = "#E74C3C",
    "Down in AD" = "#3498DB",
    "Not significant" = "gray60"
  )) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(
    title = "MA Plot: TE Expression (AD vs Control)",
    x = "Log10(Mean Expression)",
    y = "Log2 Fold Change (AD / Control)",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(ma_file, p_ma, width = 10, height = 7)
cat(sprintf("  Saved: %s\n", ma_file))

# PCA plot
cat("Creating PCA plot...\n")
pca_file <- file.path(output_dir, "PCA_plot_TE_expression.pdf")

# Variance stabilizing transformation for PCA
vsd <- vst(dds_te, blind = TRUE)

# Compute PCA
pca_data <- plotPCA(vsd, intgroup = "diagnosis", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = diagnosis)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Control" = "#2ECC71", "AD" = "#E74C3C")) +
  labs(
    title = "PCA of TE Expression",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Diagnosis"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(pca_file, p_pca, width = 8, height = 6)
cat(sprintf("  Saved: %s\n", pca_file))

# ============================================================================
# 7. COMPLETE
# ============================================================================

cat("\n=================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=================================================================\n\n")

cat("Output files:\n")
cat(sprintf("  1. Results CSV: %s\n", output_file))
cat(sprintf("  2. Volcano plot: %s\n", volcano_file))
cat(sprintf("  3. MA plot: %s\n", ma_file))
cat(sprintf("  4. PCA plot: %s\n", pca_file))
cat("\n")

cat("Summary:\n")
cat(sprintf("  Tested TEs: %d\n", n_tested))
cat(sprintf("  FDR-significant (padj < 0.05): %d\n", n_sig_fdr))
cat(sprintf("  Nominally significant (p < 0.05): %d\n", n_sig_nom))
cat("\n")
