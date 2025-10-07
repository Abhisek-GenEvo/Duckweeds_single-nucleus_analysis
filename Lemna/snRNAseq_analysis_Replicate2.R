library(Seurat)
Lemna_seurat2 <- Read10X_h5("filtered_feature_bc_matrix_P2.h5")
seurat_obj2 <- CreateSeuratObject(counts = Lemna_seurat2, project = "seurat_obj2")
seurat_obj2 <- subset(seurat_obj2, subset = nCount_RNA > 500)
--------
#Removing the doublets
library(scDblFinder)
library(SingleCellExperiment)
set.seed(123)
sce_obj2 <- as.SingleCellExperiment(seurat_obj2)
sce_obj2 <- scDblFinder(sce_obj2)
seurat_obj2$scDblFinder.class <- sce_obj2$scDblFinder.class
seurat_obj2$scDblFinder.score <- sce_obj2$scDblFinder.score
seurat_obj2 <- subset(seurat_obj2, subset = scDblFinder.class == 'singlet')
--------
seurat_obj2[["percent.cp"]] <- PercentageFeatureSet(seurat_obj2, pattern = "^CP-")
seurat_obj2[["percent.mt"]] <- PercentageFeatureSet(seurat_obj2, pattern = "^MT-")
VlnPlot(seurat_obj2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4)
FeatureScatter(seurat_obj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(seurat_obj2, file = "Lemna_seurat2.rds")
