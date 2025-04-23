library(data.table)
library(Seurat)
library(ggplot2)

#loading expr mat
data <- fread("/scratch/markov2/users/axm1810/sclc_project/adata.SCLC.010920.csv", data.table = FALSE)
cell_names <- data$Cell
expr_data <- data[, -1]
expr_mat <- t(as.matrix(expr_data))
colnames(expr_matrix) <- cell_names
#QC
seurat_object <-CreateSeuratObject(counts = expr_mat, min.features = 200, min.cells = 3)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)

seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.5)

seurat <- RunUMAP(seurat, dims = 1:15)
saveRDS(seurat, "sclc_clustered.rds")