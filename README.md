# Seurat
#ScRNAseq
install.packages(c("Seurat", "dplyr", "ggplot2", "Matrix", "patchwork", 
                   "viridis", "cowplot", "ggrepel", "circlize", "plotly"))
BiocManager::install(c("ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))

# Load Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(patchwork)
library(viridis)
library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(enrichplot)
library(plotly)


setwd("C:/Users/ASUS/Desktop/Projects/Kahnemouei/1")


# Load and Process Data
# Read the data files
barcodes <- readLines("barcodes.tsv")
features <- read.table("features.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
matrix <- readMM("matrix.mtx")

# Set up the matrix
rownames(matrix) <- make.unique(features$V2)
colnames(matrix) <- barcodes

# Define CMO groups (excluding shpi4ka groups)
cmo_groups <- list(
  "shctrl_7" = "CMO307",
  "shctrl_37" = "CMO310",
  "shvps36_7" = "CMO309",
  "shvps36_37" = "CMO312"
)

# Process CMO Data
# Find CMO feature indices
cmo_pattern <- paste0("^", unlist(cmo_groups), "$")
cmo_features_idx <- which(features$V2 %in% unlist(cmo_groups))

if (length(cmo_features_idx) > 0) {
  # Extract CMO counts and RNA counts
  cmo_mat <- matrix[cmo_features_idx, ]
  rna_mat <- matrix[-cmo_features_idx, ]
  
  # Convert sparse matrix to dense for CMO assignment
  cmo_mat <- as.matrix(cmo_mat)
  
  # Make sure row names are set correctly for CMO matrix
  rownames(cmo_mat) <- features$V2[cmo_features_idx]
  
  # Normalize CMO counts
  cmo_sums <- colSums(cmo_mat)
  cmo_norm <- t(t(cmo_mat) / cmo_sums) * 100
  
  # Find maximum CMO for each cell
  max_cmos <- apply(cmo_norm, 2, function(x) {
    cmo_name <- rownames(cmo_norm)[which.max(x)]
    return(cmo_name)
  })
  
  # Map CMO IDs to conditions
  conditions <- sapply(max_cmos, function(x) {
    cond <- names(cmo_groups)[which(unlist(cmo_groups) == x)]
    if (length(cond) == 0) return(NA)
    return(cond)
  })
  
  # Filter cells
  valid_cells <- !is.na(conditions)
  message(sprintf("Retained %d cells after CMO filtering", sum(valid_cells)))
  
  if (sum(valid_cells) == 0) {
    stop("No cells matched the CMO patterns.")
  }
  
  # Subset the RNA matrix to keep only valid cells
  rna_mat <- rna_mat[, valid_cells]
  conditions <- conditions[valid_cells]
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = rna_mat,
    project = "scRNA_analysis"
  )
  seurat_obj$condition <- conditions
  
  # Print summary of cells per condition
  print("Cells per condition:")
  print(table(conditions))
  
} else {
  stop("No CMO features found.")
}

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# QC Plots
qc_plots <- VlnPlot(seurat_obj,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3,
                    pt.size = 0.1,
                    group.by = "condition") &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(qc_plots)

# Feature scatter plot
scatter_plot <- FeatureScatter(seurat_obj, 
                               feature1 = "nCount_RNA", 
                               feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
print(scatter_plot)

# Filter cells
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 2500 & 
                       percent.mt < 5)

# Normalization and Variable Feature Selection
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 3000)


# Continue from your current point
# Instead of scaling all genes, only scale variable features
var_features <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = var_features)

# Run PCA only on variable features
seurat_obj <- RunPCA(seurat_obj, features = var_features, npcs = 50)

# Check elbow plot to determine number of PCs to use
print(ElbowPlot(seurat_obj))

# Based on elbow plot, choose appropriate number of PCs (let's say 30 for now)
n_pcs <- 30

# Run UMAP and clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)

# Visualization
p1 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = TRUE) +
  ggtitle("UMAP plot - Clusters")

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "condition") +
  ggtitle("UMAP plot - Conditions")

print(p1 | p2)

# Find markers (using only top 100 most variable genes per cluster to save memory)
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25,
                          max.cells.per.ident = 500)  # Subsample cells to save memory

# Save markers
write.csv(markers, "cluster_markers.csv")

# Generate a more memory-efficient heatmap with top markers
top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Plot heatmap with subset of cells
cells_to_plot <- sample(colnames(seurat_obj), 
                        min(5000, ncol(seurat_obj)))  # Plot max 5000 cells

print(DoHeatmap(seurat_obj, 
                features = top10$gene,
                cells = cells_to_plot) + 
        NoLegend())

# Differential expression between conditions
# Use a more memory-efficient approach with subsampling
de_results <- list()

# Function to perform DE analysis with subsampling
run_de_analysis <- function(seurat_obj, ident.1, ident.2, max_cells = 1000) {
  # Subsample cells from each condition
  cells1 <- WhichCells(seurat_obj, expression = condition == ident.1)
  cells2 <- WhichCells(seurat_obj, expression = condition == ident.2)
  
  cells1 <- sample(cells1, min(length(cells1), max_cells))
  cells2 <- sample(cells2, min(length(cells2), max_cells))
  
  # Run DE analysis on subsampled cells
  de_genes <- FindMarkers(seurat_obj,
                          cells.1 = cells1,
                          cells.2 = cells2,
                          group.by = "condition")
  
  return(de_genes)
}



# Save a reduced version of the Seurat object (optional)
# This will save only necessary data to reduce file size
seurat_obj_minimal <- DietSeurat(seurat_obj, 
                                 counts = TRUE,
                                 data = TRUE,
                                 scale.data = FALSE,
                                 features = var_features,
                                 assays = "RNA",
                                 dimreducs = c("pca", "umap"))

saveRDS(seurat_obj_minimal, "seurat_object_minimal.rds")

# Print summary statistics
print("Analysis Summary:")
print(paste("Number of cells:", ncol(seurat_obj)))
print("Cells per condition:")
print(table(seurat_obj$condition))
print("Cells per cluster:")
print(table(seurat_obj$seurat_clusters))




# Advanced visualization with enhanced plots
library(ggplot2)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# 1. Enhanced UMAP Plot with Custom Theme
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Enhanced UMAP plots
p1 <- DimPlot(seurat_obj, 
              reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 4,
              repel = TRUE) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  ggtitle("Cell Clusters") +
  custom_theme

p2 <- DimPlot(seurat_obj,
              reduction = "umap",
              group.by = "condition",
              label = FALSE) +
  scale_color_brewer(palette = "Set2") +
  ggtitle("Experimental Conditions") +
  custom_theme

combined_umap <- p1 | p2
print(combined_umap)

# 2. Cluster Composition Barplot
cluster_comp <- table(seurat_obj$seurat_clusters, seurat_obj$condition)
cluster_comp_df <- as.data.frame(cluster_comp)
colnames(cluster_comp_df) <- c("Cluster", "Condition", "Count")

# Calculate percentages
cluster_comp_df <- cluster_comp_df %>%
  group_by(Cluster) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Create stacked barplot
composition_plot <- ggplot(cluster_comp_df, 
                           aes(x = Cluster, 
                               y = Percentage, 
                               fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(title = "Cluster Composition by Condition",
       x = "Cluster",
       y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

print(composition_plot)

# 3. Advanced Violin Plot for Key Markers
top_markers <- c("RPS8", "PTMA", "RPL6", "RPL28", "FTL")
markers_violin <- VlnPlot(seurat_obj,
                          features = top_markers,
                          stack = TRUE,
                          flip = TRUE,
                          sort = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Expression of Top Markers Across Clusters") +
  custom_theme

print(markers_violin)

# 4. Enhanced Feature Plot for Top Genes
feature_plot <- FeaturePlot(seurat_obj,
                            features = top_markers[1:4],
                            ncol = 2,
                            pt.size = 0.1,
                            min.cutoff = "q10",
                            max.cutoff = "q90") &
  scale_color_viridis(option = "magma") &
  custom_theme

print(feature_plot)

# 5. Enhanced Heatmap with Top Markers
top_genes_expr <- GetAssayData(seurat_obj, slot = "data")[top_markers, ]
cluster_avg <- aggregate(t(top_genes_expr), 
                         by = list(seurat_obj$seurat_clusters),
                         FUN = mean)
row.names(cluster_avg) <- cluster_avg$Group.1
cluster_avg$Group.1 <- NULL
cluster_avg <- t(cluster_avg)

# Create enhanced heatmap
col_fun <- colorRamp2(c(min(cluster_avg), mean(cluster_avg), max(cluster_avg)),
                      c("#313695", "#FFFFBF", "#A50026"))

hm <- Heatmap(cluster_avg,
              name = "Expression",
              col = col_fun,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 10),
              column_title = "Clusters",
              row_title = "Genes",
              heatmap_legend_param = list(title = "Expression"))

print(hm)

# Save all plots
pdf("advanced_visualization.pdf", width = 15, height = 20)
print(combined_umap)
print(composition_plot)
print(markers_violin)
print(feature_plot)
print(hm)
dev.off()







# Advanced DE and Complex Visualizations
library(ggplot2)
library(patchwork)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)

# 1. Enhanced Volcano Plot for DE genes between conditions
de_shvps36_7 <- FindMarkers(seurat_obj,
                            ident.1 = "shvps36_7",
                            ident.2 = "shctrl_7",
                            group.by = "condition",
                            test.use = "wilcox",
                            min.pct = 0.25)

# Create enhanced volcano plot
volcano_plot <- EnhancedVolcano(de_shvps36_7,
                                lab = rownames(de_shvps36_7),
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                title = 'Differential Expression: shvps36_7 vs shctrl_7',
                                pCutoff = 0.05,
                                FCcutoff = 0.5,
                                pointSize = 3.0,
                                labSize = 4.0,
                                col = c('gray', 'blue', 'red', 'green'),
                                colAlpha = 0.5,
                                legendPosition = 'right',
                                drawConnectors = TRUE,
                                widthConnectors = 0.2)

print(volcano_plot)

# 2. Time Course Expression Plot
# Compare expression changes over time points (7 vs 37)
timepoint_comparison <- function(gene) {
  expr_data <- FetchData(seurat_obj, 
                         vars = c(gene, "condition"))
  colnames(expr_data)[1] <- "expression"
  
  ggplot(expr_data, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.5) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    ggtitle(paste("Expression of", gene, "across conditions")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
}

# Create time course plots for top genes
top_genes <- c("RPS8", "PTMA", "RPL6", "RPL28", "FTL")
time_course_plots <- lapply(top_genes, timepoint_comparison)
combined_time_plots <- wrap_plots(time_course_plots, ncol = 2)
print(combined_time_plots)

# 3. Advanced Dot Plot with Multiple Features
dot_plot <- DotPlot(seurat_obj, 
                    features = top_genes,
                    group.by = "seurat_clusters") +
  scale_color_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

print(dot_plot)

# 4. Complex Heatmap with Multiple Annotations
# Prepare expression data
expr_matrix <- GetAssayData(seurat_obj, slot = "scale.data")[top_genes, ]
expr_matrix <- t(expr_matrix)


# Create complex heatmap
complex_heatmap <- pheatmap(expr_matrix,
                            annotation_col = annotation_col,
                            annotation_colors = ann_colors,
                            show_colnames = FALSE,
                            clustering_distance_rows = "correlation",
                            clustering_distance_cols = "correlation",
                            main = "Expression Patterns with Annotations")

print(complex_heatmap)

# 5. Cluster Relationship Network Plot
# Calculate cluster correlations
cluster_averages <- AverageExpression(seurat_obj, group.by = "seurat_clusters")
cluster_cor <- cor(cluster_averages$RNA)

# Create network data
library(igraph)
library(ggraph)

# Create graph from correlation matrix
graph <- graph_from_adjacency_matrix(
  cluster_cor,
  weighted = TRUE,
  mode = "undirected",
  diag = FALSE
)

# Create network plot
network_plot <- ggraph(graph, layout = 'fr') +
  geom_edge_link(aes(edge_alpha = abs(weight), edge_width = abs(weight)), 
                 edge_colour = "grey50") +
  geom_node_point(size = 10, color = viridis(vmax = length(V(graph)))) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() +
  ggtitle("Cluster Relationship Network")

print(network_plot)

# 6. Gene Expression Distribution Ridgeline Plot
library(ggridges)

# Prepare data for ridgeline plot
expr_data <- FetchData(seurat_obj, vars = c(top_genes, "seurat_clusters"))
expr_long <- tidyr::gather(expr_data, key = "gene", value = "expression", -seurat_clusters)

# Create ridgeline plot
ridge_plot <- ggplot(expr_data, aes(x = expression, y = gene, fill = gene)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Expression Distribution of Top Genes Across Clusters") +
  theme(legend.position = "none")

print(ridge_plot)

# 7. Circular Barplot for Cluster Proportions
library(tidyr)
library(circlize)

# Prepare data for circular plot
cluster_props <- prop.table(table(seurat_obj$seurat_clusters))
cluster_df <- data.frame(
  cluster = names(cluster_props),
  proportion = as.numeric(cluster_props)
)

# Create circular barplot
circos.clear()
circos.par(start.degree = 90, gap.degree = 2)
circos.initialize("a", xlim = c(0, nrow(cluster_df)))

# Plot circular bars
circos.track(
  ylim = c(0, max(cluster_df$proportion)),
  track.height = 0.3,
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.barplot(
      cluster_df$proportion,
      pos = seq_len(nrow(cluster_df)) - 0.5,
      col = viridis(nrow(cluster_df)),
      border = NA
    )
  }
)

# Add cluster labels
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      seq_len(nrow(cluster_df)) - 0.5,
      rep(max(cluster_df$proportion) + 0.02, nrow(cluster_df)),
      cluster_df$cluster,
      facing = "clockwise",
      adj = c(0, 0.5)
    )
  }
)

# Save all plots
pdf("advanced_de_visualization.pdf", width = 15, height = 25)
print(volcano_plot)
print(combined_time_plots)
print(dot_plot)
print(complex_heatmap)
print(network_plot)
print(ridge_plot)
dev.off()

# Print completion message
message("All advanced visualizations have been created and saved!")











#####  PI4KA knockdown cells (shpi4ka) cluster compared to controls #####
# Define PI4KA study CMO groups
cmo_groups <- list(
  "shctrl_7" = "CMO307",
  "shctrl_37" = "CMO310",
  "shpi4ka_7" = "CMO308",
  "shpi4ka_37" = "CMO311"
)

# Process CMO Data
cmo_features_idx <- which(features$V2 %in% unlist(cmo_groups))
if (length(cmo_features_idx) > 0) {
  # Extract CMO and RNA counts
  cmo_mat <- matrix[cmo_features_idx, ]
  rna_mat <- matrix[-cmo_features_idx, ]
  
  cmo_mat <- as.matrix(cmo_mat)
  rownames(cmo_mat) <- features$V2[cmo_features_idx]
  
  # Normalize CMO counts
  cmo_sums <- colSums(cmo_mat)
  cmo_norm <- t(t(cmo_mat) / cmo_sums) * 100
  
  # Assign cells to groups
  max_cmos <- apply(cmo_norm, 2, function(x) {
    cmo_name <- rownames(cmo_norm)[which.max(x)]
    return(cmo_name)
  })
  
  conditions <- sapply(max_cmos, function(x) {
    cond <- names(cmo_groups)[which(unlist(cmo_groups) == x)]
    if (length(cond) == 0) return(NA)
    return(cond)
  })
  
  valid_cells <- !is.na(conditions)
  
  # Create Seurat object with only PI4KA study cells
  rna_mat <- rna_mat[, valid_cells]
  conditions <- conditions[valid_cells]
  
  seurat_obj <- CreateSeuratObject(
    counts = rna_mat,
    project = "PI4KA_analysis"
  )
  seurat_obj$condition <- conditions
}

# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 2500 & 
                       percent.mt < 5)

# Normalize and Find Variable Features
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)

# Scale variable features only
var_features <- VariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = var_features)

# Dimensional reduction and clustering
seurat_obj <- RunPCA(seurat_obj, features = var_features)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Visualization
p1 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "seurat_clusters", 
              label = TRUE) +
  ggtitle("UMAP plot - Clusters")

p2 <- DimPlot(seurat_obj, 
              reduction = "umap", 
              group.by = "condition") +
  ggtitle("UMAP plot - Conditions")

print(p1 | p2)

# Calculate cluster composition
cluster_comp <- table(seurat_obj$seurat_clusters, seurat_obj$condition)
cluster_percentages <- prop.table(cluster_comp, margin = 1) * 100

# Find enriched clusters for PI4KA knockdown
enriched_clusters <- apply(cluster_percentages, 1, function(x) {
  pi4ka_percent <- sum(x[grep("shpi4ka", names(x))])
  return(pi4ka_percent > 60)  # Clusters with >60% PI4KA cells
})

enriched_cluster_ids <- names(enriched_clusters)[enriched_clusters]

# Find markers for enriched clusters
markers_list <- list()
for(cluster in enriched_cluster_ids) {
  markers <- FindMarkers(seurat_obj,
                         ident.1 = cluster,
                         group.by = "seurat_clusters",
                         min.pct = 0.25,
                         logfc.threshold = 0.25)
  markers_list[[cluster]] <- markers
}

# Save results
saveRDS(seurat_obj, "pi4ka_analysis.rds")
write.csv(cluster_percentages, "cluster_composition.csv")
for(cluster in names(markers_list)) {
  write.csv(markers_list[[cluster]], 
            file = paste0("markers_cluster_", cluster, ".csv"))
}

# Create enrichment visualization
cluster_comp_df <- as.data.frame(cluster_percentages)
colnames(cluster_comp_df) <- c("Cluster", "Condition", "Percentage")

ggplot(cluster_comp_df, aes(x = Cluster, y = Percentage, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cluster Composition",
       y = "Percentage",
       x = "Cluster")


# Advanced visualizations for PI4KA analysis
library(ggplot2)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(RColorBrewer)

# 1. Enhanced UMAP with custom theme
custom_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Create enhanced UMAP plots
p1 <- DimPlot(seurat_obj, 
              reduction = "umap",
              group.by = "seurat_clusters",
              label = TRUE,
              label.size = 4) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  ggtitle("Cell Clusters") +
  custom_theme

p2 <- DimPlot(seurat_obj,
              reduction = "umap",
              group.by = "condition") +
  scale_color_brewer(palette = "Set2") +
  ggtitle("Experimental Conditions") +
  custom_theme

combined_umap <- p1 | p2
print(combined_umap)

# 2. Circular Cluster Composition Plot
cluster_comp <- table(seurat_obj$seurat_clusters, seurat_obj$condition)
cluster_percentages <- prop.table(cluster_comp, margin = 1) * 100

# Prepare data for circular plot
melted_data <- as.data.frame(as.table(cluster_percentages))
colnames(melted_data) <- c("Cluster", "Condition", "Percentage")

# Create circular barplot
circos.clear()
circos.par(start.degree = 90, gap.degree = 5)

# Set up circular layout
cluster_order <- unique(melted_data$Cluster)
condition_colors <- brewer.pal(length(unique(melted_data$Condition)), "Set2")
names(condition_colors) <- unique(melted_data$Condition)

circos.initialize(factors = cluster_order, xlim = c(0, 100))

# Add tracks
circos.track(
  ylim = c(0, 1),
  factors = cluster_order,
  track.height = 0.3,
  bg.border = NA,
  panel.fun = function(x, y) {
    cluster <- CELL_META$sector.index
    cluster_data <- subset(melted_data, Cluster == cluster)
    cumsum_prev <- 0
    for(cond in unique(cluster_data$Condition)) {
      value <- subset(cluster_data, Condition == cond)$Percentage
      circos.rect(
        cumsum_prev, 0,
        cumsum_prev + value, 0.8,
        col = condition_colors[cond]
      )
      cumsum_prev <- cumsum_prev + value
    }
  }
)

# Add labels
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      x = 50,
      y = 0,
      labels = CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE
    )
  }
)

# 3. Advanced Heatmap with Multiple Annotations
# Get top markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)
top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# Prepare expression matrix
expr_matrix <- AverageExpression(seurat_obj, assays = "RNA", features = unique(top10$gene))$RNA
expr_matrix <- t(scale(t(expr_matrix)))

# Create annotation colors
cluster_colors <- viridis(length(unique(seurat_obj$seurat_clusters)))
names(cluster_colors) <- unique(seurat_obj$seurat_clusters)
condition_colors <- brewer.pal(length(unique(seurat_obj$condition)), "Set2")
names(condition_colors) <- unique(seurat_obj$condition)

# Create complex heatmap
heatmap <- Heatmap(expr_matrix,
                   name = "Expression",
                   col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8),
                   column_title = "Clusters",
                   row_title = "Genes")

print(heatmap)

# 4. Ridge Plot for Key Genes
key_genes <- unique(top10$gene)[1:10]
expr_data <- FetchData(seurat_obj, vars = c(key_genes, "seurat_clusters"))
expr_long <- tidyr::gather(expr_data, key = "gene", value = "expression", -seurat_clusters)

ridge_plot <- ggplot(expr_long, aes(x = expression, y = gene, fill = gene)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis_d() +
  custom_theme +
  labs(title = "Expression Distribution of Key Genes") +
  theme(legend.position = "none")

print(ridge_plot)

# 5. Bubble Plot
cluster_averages <- AverageExpression(seurat_obj, features = key_genes, group.by = "seurat_clusters")$RNA
percent_expressed <- AverageExpression(seurat_obj, features = key_genes, group.by = "seurat_clusters", slot = "data")$RNA > 0

bubble_data <- data.frame(
  Gene = rep(rownames(cluster_averages), ncol(cluster_averages)),
  Cluster = rep(colnames(cluster_averages), each = nrow(cluster_averages)),
  Expression = as.vector(cluster_averages),
  PercentExpressed = as.vector(percent_expressed) * 100
)

bubble_plot <- ggplot(bubble_data, 
                      aes(x = Cluster, y = Gene, size = PercentExpressed, color = Expression)) +
  geom_point() +
  scale_color_viridis() +
  scale_size_continuous(range = c(1, 10)) +
  custom_theme +
  labs(title = "Gene Expression Patterns Across Clusters",
       size = "% Expressed",
       color = "Expression Level")

print(bubble_plot)

# Save all plots
pdf("advanced_pi4ka_visualization.pdf", width = 15, height = 20)
print(combined_umap)
print(capture.output(circos.clear()))
print(heatmap)
print(ridge_plot)
print(bubble_plot)
dev.off()
