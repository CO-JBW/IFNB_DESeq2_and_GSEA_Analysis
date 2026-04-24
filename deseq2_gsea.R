#!/usr/bin/env Rscript

# Cell-type-specific DESeq2 and GSEA analysis for IFNB dataset
# Input: /data/ifnb_qc.rds (QC'd and annotated Seurat object)
# Outputs: Cell-type-specific DE results, GSEA results, and visualizations to /results

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(EnhancedVolcano)
})

# Create results directory
dir.create("/results", showWarnings = FALSE, recursive = TRUE)

# Load QC'd and annotated data
message("Loading annotated IFNB data...")
ifnb <- readRDS("/data/ifnb_qc/ifnb_qc.rds")

message(paste("Total cells:", ncol(ifnb)))
message(paste("CTRL cells:", sum(ifnb$stim == "CTRL")))
message(paste("STIM cells:", sum(ifnb$stim == "STIM")))

# Display cell type distribution
message("\nCell type distribution:")
print(table(ifnb$cell_type, ifnb$stim))

# ============================================
# CELL TYPE SELECTION
# ============================================
message("\n=== Selecting cell type for analysis ===")

# Identify cell types with sufficient cells in both conditions
cell_type_counts <- table(ifnb$cell_type, ifnb$stim)
min_cells_per_condition <- 50

cell_types_with_both <- rownames(cell_type_counts)[
  cell_type_counts[, "CTRL"] >= min_cells_per_condition & 
  cell_type_counts[, "STIM"] >= min_cells_per_condition
]

# Exclude "Unknown" cells from analysis (focus on well-annotated PBMC populations)
cell_types_with_both <- cell_types_with_both[cell_types_with_both != "Unknown"]

if (length(cell_types_with_both) == 0) {
  stop(paste("No well-annotated cell types have >=", min_cells_per_condition, "cells in both conditions."))
}

# Choose the most abundant suitable cell type
selected_cell_type <- names(which.max(rowSums(cell_type_counts[cell_types_with_both, , drop = FALSE])))

message(paste("\nSelected cell type for analysis:", selected_cell_type))
message(paste("  CTRL cells:", cell_type_counts[selected_cell_type, "CTRL"]))
message(paste("  STIM cells:", cell_type_counts[selected_cell_type, "STIM"]))

# Save cell type selection info
write.csv(
  data.frame(
    selected_cell_type = selected_cell_type,
    ctrl_cells = cell_type_counts[selected_cell_type, "CTRL"],
    stim_cells = cell_type_counts[selected_cell_type, "STIM"]
  ),
  "/results/selected_cell_type.csv",
  row.names = FALSE
)

# Subset to selected cell type
ifnb_subset <- subset(ifnb, subset = cell_type == selected_cell_type)
message(paste("\nSubset to", selected_cell_type, ":", ncol(ifnb_subset), "cells"))

# Visualize the selected cell type
p_umap_selected <- DimPlot(ifnb, cells.highlight = list(Selected = colnames(ifnb_subset)),
                           cols.highlight = "red", cols = "grey90") +
  ggtitle(paste("Selected Cell Type:", selected_cell_type))
ggsave("/results/umap_selected_celltype.png", p_umap_selected, width = 10, height = 8, dpi = 300)

# ============================================
# PSEUDO-BULK AGGREGATION (CELL-TYPE-SPECIFIC)
# ============================================
message("\n=== Performing cell-type-specific pseudo-bulk aggregation ===")

set.seed(42)

# Create pseudo-replicates by randomly assigning cells within each condition
cells_ctrl <- which(ifnb_subset$stim == "CTRL")
cells_stim <- which(ifnb_subset$stim == "STIM")

# Split into 3 pseudo-replicates per condition
n_reps <- 3
ifnb_subset$replicate <- NA
ifnb_subset$replicate[cells_ctrl] <- paste0("CTRL_rep", sample(1:n_reps, length(cells_ctrl), replace = TRUE))
ifnb_subset$replicate[cells_stim] <- paste0("STIM_rep", sample(1:n_reps, length(cells_stim), replace = TRUE))

message("Pseudo-replicates created:")
print(table(ifnb_subset$replicate))

# Aggregate counts by replicate
counts_matrix <- GetAssayData(ifnb_subset, layer = "counts", assay = "RNA")

pseudobulk_counts <- sapply(unique(ifnb_subset$replicate), function(rep) {
  cells_in_rep <- which(ifnb_subset$replicate == rep)
  Matrix::rowSums(counts_matrix[, cells_in_rep])
})

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = colnames(pseudobulk_counts),
  condition = ifelse(grepl("CTRL", colnames(pseudobulk_counts)), "CTRL", "STIM"),
  replicate = gsub(".*_rep", "rep", colnames(pseudobulk_counts)),
  cell_type = selected_cell_type,
  row.names = colnames(pseudobulk_counts)
)

write.csv(sample_metadata, "/results/pseudobulk_sample_metadata.csv", row.names = TRUE)

message("\nPseudo-bulk matrix dimensions:")
message(paste("  Genes:", nrow(pseudobulk_counts)))
message(paste("  Samples:", ncol(pseudobulk_counts)))
message(paste("  Cell type:", selected_cell_type))

# ============================================
# DESeq2 ANALYSIS
# ============================================
message("\n=== Running DESeq2 analysis ===")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = sample_metadata,
  design = ~ condition
)

# Pre-filtering: keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
message(paste("Genes retained after filtering:", nrow(dds)))

# Set reference level
dds$condition <- relevel(dds$condition, ref = "CTRL")

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "STIM", "CTRL"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$cell_type <- selected_cell_type

# Reorder columns
res_df <- res_df[, c("gene", "cell_type", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Save full results
write.csv(res_df, "/results/deseq2_results.csv", row.names = FALSE)

# Filter significant genes (padj < 0.05)
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
write.csv(sig_genes, "/results/de_genes_significant.csv", row.names = FALSE)

message(paste("\nSignificant genes (padj < 0.05) in", selected_cell_type, ":", nrow(sig_genes)))
message(paste("  Upregulated in STIM:", sum(sig_genes$log2FoldChange > 0)))
message(paste("  Downregulated in STIM:", sum(sig_genes$log2FoldChange < 0)))

# ============================================
# VISUALIZATIONS
# ============================================
message("\n=== Generating visualizations ===")

# Volcano plot
message("Creating volcano plot...")
p_volcano <- EnhancedVolcano(res_df,
  lab = res_df$gene,
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  title = paste('STIM vs CTRL in', selected_cell_type),
  subtitle = 'DESeq2 Differential Expression',
  pointSize = 2.0,
  labSize = 3.0,
  legendPosition = 'right',
  maxoverlapsConnectors = 20
)
ggsave("/results/volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 300)

# MA plot
message("Creating MA plot...")
png("/results/ma_plot.png", width = 10, height = 8, units = "in", res = 300)
plotMA(res, ylim = c(-5, 5), main = paste("MA Plot:", selected_cell_type, "STIM vs CTRL"))
dev.off()

# Heatmap of top 50 DE genes
message("Creating heatmap of top DE genes...")
if (nrow(sig_genes) > 0) {
  top_n_genes <- min(50, nrow(sig_genes))
  top_genes <- head(sig_genes$gene, top_n_genes)
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  top_genes_counts <- normalized_counts[top_genes, , drop = FALSE]
  
  # Log transform for visualization
  top_genes_counts_log <- log2(top_genes_counts + 1)
  
  # Create annotation
  annotation_col <- data.frame(
    Condition = sample_metadata$condition,
    row.names = rownames(sample_metadata)
  )
  
  png("/results/heatmap_top_genes.png", width = 10, height = 12, units = "in", res = 300)
  pheatmap(top_genes_counts_log,
           annotation_col = annotation_col,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           scale = "row",
           show_rownames = ifelse(top_n_genes <= 50, TRUE, FALSE),
           fontsize_row = 8,
           main = paste("Top", top_n_genes, "DE Genes in", selected_cell_type))
  dev.off()
} else {
  message("No significant genes found for heatmap.")
}

# ============================================
# GSEA ANALYSIS
# ============================================
message("\n=== Running GSEA analysis ===")

# Prepare ranked gene list
res_df_gsea <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$log2FoldChange), ]
res_df_gsea$rank_metric <- sign(res_df_gsea$log2FoldChange) * (-log10(res_df_gsea$pvalue))

# Remove non-finite values (Inf, -Inf, NaN)
res_df_gsea <- res_df_gsea[is.finite(res_df_gsea$rank_metric), ]

ranked_genes <- setNames(res_df_gsea$rank_metric, res_df_gsea$gene)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

message(paste("Ranked genes for GSEA:", length(ranked_genes)))

# Get MSigDB Hallmark gene sets (human)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

message(paste("Hallmark gene sets:", length(hallmark_list)))

# Run fgsea (using fgseaMultilevel, the recommended method)
set.seed(42)
fgsea_results <- fgsea(
  pathways = hallmark_list,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500
)

# Sort by adjusted p-value
fgsea_results <- fgsea_results[order(fgsea_results$padj), ]
fgsea_results$cell_type <- selected_cell_type

# Save results
fgsea_df <- as.data.frame(fgsea_results)
fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, function(x) paste(x, collapse = ";"))
write.csv(fgsea_df, "/results/gsea_results.csv", row.names = FALSE)

message(paste("\nSignificant pathways (padj < 0.05):", sum(fgsea_results$padj < 0.05, na.rm = TRUE)))

# Plot top enriched pathways
message("Creating GSEA enrichment plots...")
pdf("/results/gsea_enrichment_plots.pdf", width = 10, height = 6)

top_pathways <- head(fgsea_results$pathway[fgsea_results$padj < 0.05], 10)
if (length(top_pathways) > 0) {
  for (pathway in top_pathways) {
    p <- plotEnrichment(hallmark_list[[pathway]], ranked_genes) +
      ggtitle(paste(pathway, "-", selected_cell_type)) +
      theme_minimal()
    print(p)
  }
}

dev.off()

# Create summary barplot of top pathways
sig_pathways <- fgsea_results[fgsea_results$padj < 0.05, ]
if (nrow(sig_pathways) > 0) {
  top_n <- min(20, nrow(sig_pathways))
  plot_data <- head(sig_pathways, top_n)
  plot_data$pathway <- gsub("HALLMARK_", "", plot_data$pathway)
  
  p_bar <- ggplot(plot_data, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue") +
    labs(title = paste("Top GSEA Hallmark Pathways in", selected_cell_type),
         x = "Pathway",
         y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave("/results/gsea_barplot.png", p_bar, width = 10, height = 8, dpi = 300)
}

# ============================================
# SUMMARY
# ============================================
message("\n=== Analysis Complete ===")
message("\nCell-type-specific analysis performed on:", selected_cell_type)
message("\nOutputs saved to /results:")
message("  - selected_cell_type.csv: Cell type used for analysis")
message("  - umap_selected_celltype.png: UMAP highlighting selected cells")
message("  - pseudobulk_sample_metadata.csv: Pseudo-bulk sample information")
message("  - deseq2_results.csv: Full DESeq2 results")
message("  - de_genes_significant.csv: Significant DE genes (padj < 0.05)")
message("  - volcano_plot.png: Volcano plot visualization")
message("  - ma_plot.png: MA plot")
message("  - heatmap_top_genes.png: Heatmap of top DE genes")
message("  - gsea_results.csv: GSEA pathway enrichment results")
message("  - gsea_enrichment_plots.pdf: Enrichment plots for top pathways")
message("  - gsea_barplot.png: Bar plot of top pathways by NES")

message("\nTop 10 upregulated genes in", selected_cell_type, ":")
if (nrow(sig_genes) > 0) {
  top_up <- head(sig_genes[sig_genes$log2FoldChange > 0, c("gene", "log2FoldChange", "padj")], 10)
  print(top_up)
}

message("\nTop 10 enriched pathways:")
if (nrow(fgsea_results) > 0) {
  print(head(fgsea_results[, c("pathway", "NES", "padj")], 10))
}
