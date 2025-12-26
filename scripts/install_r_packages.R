#!/usr/bin/env Rscript

# Install required packages for TE differential expression analysis

cat("Checking and installing required R packages...\n\n")

# Function to install if not present
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    if (bioc) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "http://cran.us.r-project.org")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
    cat(sprintf("✓ %s installed\n\n", pkg))
  } else {
    cat(sprintf("✓ %s already installed\n", pkg))
  }
}

# Install packages
install_if_missing("ggplot2")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("RColorBrewer")
install_if_missing("pheatmap")
install_if_missing("DESeq2", bioc = TRUE)

cat("\n=============================================================\n")
cat("All required packages installed successfully!\n")
cat("=============================================================\n")
