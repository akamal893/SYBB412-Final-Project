
```{r setup}
library(Seurat)
library(ggplot2)

#loading expr mat

data <- read.csv("/scratch/markov2/users/axm1810/sclc_project/adata.SCLC.010920.csv")

cell_names <- data$Cell
expr_data <- data[, -1]
expr_mat <- t(as.matrix(expr_data))
colnames(expr_mat) <- cell_names
#QC
seurat <- CreateSeuratObject(counts = expr_mat, min.features = 200, min.cells = 3)
str(data)

meta_data <- readRDS("/scratch/markov2/users/axm1810/seurat_normalize.rds")
head(rownames(meta_data))      # should be long barcodes like "RU1231A_120703424294126"
head(colnames(seurat))         # same here – these must match exactly

seurat$SCLC_subtype <- meta_data$SCLC_subtype
seurat$procedure <- meta_data$procedure

ncol(seurat)              
length(meta_data$SCLC_subtype)    



seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)


# Visualize QC metrics
pdf("QC_metrics.pdf")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.5)



seurat <- RunUMAP(seurat, dims = 1:15)
saveRDS(seurat, "sclc_clustered.rds")


seurat <- readRDS("/scratch/markov2/users/axm1810/sclc_project/sclc_clustered.rds")
pdf("sclc_umap_clusters.pdf", width = 8, height = 6)
p <- DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("UMAP of SCLC Clusters")
print(p)

## Marker Gene
marker_genes <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
top <- marker_genes %>%
    group_by(cluster) %>%
    top_n(n=1, wt = avg_log2FC)
for (i in 1:nrow(top)){
  gene <- top$gene[i]
  cluster <- top$cluster[i]
  if (gene %in% rownames(seurat)) {
    print(FeaturePlot(seurat, features = gene, reduction ="umap") + ggtitle(paste("Top Marker Gene for", cluster, ":", gene)))
  }
  else {
    cat("Skip", gene, "not in set\n")
  }
    }
dev.off()
##Subtype Distribution 
pdf("SCLC_Subtype_Distribution.pdf", width = 8, height = 6)
if ("SCLC_subtype" %in% colnames(seurat@meta.data)) {
  DimPlot(seurat, group.by = "SCLC_subtype", reduction = "umap") + ggtitle("SCLC Subtype (N vs A) Distributiion Plot")
} else {
  cat("not found in metadata.\n")
}

# end of PDF builder
dev.off()

## Subtype Counts and Percents
subtype_count <- table(seurat$seurat_clusters, seurat$SCLC_subtype)
subtype_df <- as.data.frame.matrix(subtype_count)
View(subtype_df)
subtype_pt <- prop.table(subtype_count, margin = 1) * 100
subtype_pt_df <- as.data.frame.matrix(round(subtype_pt, 2))
View(subtype_pt_df)
write.csv(subtype_percent_df, "subtype_counts_percents_by_cluster.csv")
#Plots for Subtype Counts/Percents
install.packages("pheatmap")
library(pheatmap)

## Enrichment Pathway
library(dplyr)
marker_genes <- readRDS("marker_genes_optimized.rds")
head(marker_genes)

gene_cluster1 <- marker_genes %>% filter(cluster == 1 & p_val_adj < 0.05) %>% pull(gene)

library(clusterProfiler)
library(stringr)

library(org.Hs.eg.db)

entrez_ids <- bitr(gene_cluster1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)
pdf("cluster1_enrichment_dotplot.pdf", width = 8, height = 6)
dotplot(ego, showCategory = 20) + ggtitle("GO Enrichment: Cluster 1") + theme(axis.text.y = element_text(size=10)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 40))
dev.off()






seurat <- readRDS("/scratch/markov2/users/axm1810/seurat_clusters_and_subtypes.rds")
head(seurat@meta.data)


meta_export <- seurat@meta.data

meta_export$gene_id <- rownames(seurat)
saveRDS(meta_export, file = "~/ondemand/seurat_metadata_with_clusters.rds")


```
