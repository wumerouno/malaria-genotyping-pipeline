# R script for differential gene expression analysis of Plasmodium falciparum RNA-seq data
# Uses tximport to load Salmon abundances, runs DESeq2, and generates plots.

# Set output options to avoid encoding issues on Windows
options(encoding = "UTF-8")

# 1. Load Libraries
message("Loading required R packages...")
suppressPackageStartupMessages({
  library(yaml)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

# 2. Load Configuration and Metadata
message("Loading configurations and sample sheet...")
config <- yaml::read_yaml("config/config.yaml")

samples <- read.delim(config$samples_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample_id

# Verify quantification directories exist
quant_files <- file.path(config$quant_dir, samples$sample_id, "quant.sf")
names(quant_files) <- samples$sample_id

for (i in seq_along(quant_files)) {
  if (!file.exists(quant_files[i])) {
    stop(paste("Salmon quantification file not found at:", quant_files[i]))
  }
}

# 3. Dynamic tx2gene Mapping Creation
message("Generating transcript-to-gene mapping...")
first_quant <- read.delim(quant_files[1], header = TRUE, stringsAsFactors = FALSE)
tx_names <- first_quant$Name

# In PlasmoDB, transcripts are named like PF3D7_0100100.1 or PF3D7_0100100.1-p1
# The gene ID is the part before the first dot (e.g., PF3D7_0100100)
gene_ids <- sub("\\..*", "", tx_names)

tx2gene <- data.frame(
  TXNAME = tx_names,
  GENEID = gene_ids,
  stringsAsFactors = FALSE
)

# 4. Import Salmon Abundances
message("Importing Salmon abundances using tximport...")
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

# 5. Initialize DESeq2 Dataset
message("Running DESeq2 pipeline...")
# Convert design formula string to formula object
design_formula <- as.formula(config$design_formula)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)

# Filter low-count genes (retain genes with at least 10 reads in total across samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(paste("Filtered dataset contains", sum(keep), "genes out of", length(keep)))

# Set factor levels for contrast baseline
factor_col <- gsub("~", "", config$design_formula)
factor_col <- trimws(factor_col)
dds[[factor_col]] <- factor(dds[[factor_col]], levels = c(config$reference_level, config$comparison_level))

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c(factor_col, config$comparison_level, config$reference_level))

# Perform log fold change shrinkage (lfcShrink) to adjust for low-count dispersion
res_shrunk <- lfcShrink(dds, coef = paste0(factor_col, "_", config$comparison_level, "_vs_", config$reference_level), type = "normal")

# 6. Save Differential Expression Results
message("Saving results...")
dir.create(config$de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$plot_dir, recursive = TRUE, showWarnings = FALSE)

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
# Put gene_id as first column
res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]

write.csv(res_df, file = file.path(config$de_dir, "diff_exp_results.csv"), row.names = FALSE)

# Save significant DEGs
sig_degs <- subset(res_df, padj < config$padj_cutoff & abs(log2FoldChange) >= config$lfc_cutoff)
write.csv(sig_degs, file = file.path(config$de_dir, "sig_diff_exp_results.csv"), row.names = FALSE)
message(paste("Found", nrow(sig_degs), "statistically significant DEGs (padj <", config$padj_cutoff, ", |log2FC| >=", config$lfc_cutoff, ")"))

# 7. Generate Visualizations
message("Generating publication-quality plots...")

# 7a. PCA Plot
# Transform counts using Variance Stabilizing Transformation (vst)
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = factor_col, returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = .data[[factor_col]], label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title = "Principal Component Analysis (PCA)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(config$plot_dir, "pca_plot.png"), plot = pca_plot, width = 7, height = 6, dpi = 300)

# 7b. Volcano Plot
res_volcano <- res_df
res_volcano$sig <- "Not Significant"
res_volcano$sig[res_volcano$padj < config$padj_cutoff & res_volcano$log2FoldChange >= config$lfc_cutoff] <- "Upregulated"
res_volcano$sig[res_volcano$padj < config$padj_cutoff & res_volcano$log2FoldChange <= -config$lfc_cutoff] <- "Downregulated"
res_volcano$sig <- factor(res_volcano$sig, levels = c("Upregulated", "Downregulated", "Not Significant"))

# Label top 10 genes by p-value
res_volcano$label <- NA
top10 <- head(res_volcano$gene_id[res_volcano$sig != "Not Significant"], 10)
res_volcano$label[res_volcano$gene_id %in% top10] <- res_volcano$gene_id[res_volcano$gene_id %in% top10]

volcano_plot <- ggplot(res_volcano, aes(x = log2FoldChange, y = -log10(padj), color = sig, label = label)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not Significant" = "#999999")) +
  geom_vline(xintercept = c(-config$lfc_cutoff, config$lfc_cutoff), linetype = "dashed", color = "darkgray") +
  geom_hline(yintercept = -log10(config$padj_cutoff), linetype = "dashed", color = "darkgray") +
  geom_text_repel(size = 3.2, color = "black", max.overlaps = 15) +
  labs(
    title = paste("Volcano Plot:", config$comparison_level, "vs", config$reference_level),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(config$plot_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6, dpi = 300)

# 7c. Heatmap of Top 50 Genes
top50_genes <- head(res_df$gene_id, 50)
vsd_counts <- assay(vsd)[top50_genes, ]

# Mean center the rows (Z-score)
z_counts <- t(scale(t(vsd_counts)))

# Annotation details
annotation_df <- data.frame(
  Stage = samples[[factor_col]],
  row.names = samples$sample_id
)

# Render heatmap to file
png(file.path(config$plot_dir, "heatmap_top50.png"), width = 800, height = 1000, res = 120)
pheatmap(
  z_counts,
  annotation_col = annotation_df,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "none",
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(50),
  main = "Top 50 Differentially Expressed Genes (Z-score)"
)
dev.off()

message("Analysis completed successfully. All outputs generated in 'results/'")
