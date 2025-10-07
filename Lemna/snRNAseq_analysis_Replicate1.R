library(Seurat)
Lemna_seurat1 <- Read10X_h5("filtered_feature_bc_matrix_P1.h5")
seurat_obj1 <- CreateSeuratObject(counts = Lemna_seurat1, project = "seurat_obj1")
seurat_obj1 <- subset(seurat_obj1, subset = nCount_RNA > 500)

library(scDblFinder)
library(SingleCellExperiment)
set.seed(123)
sce_obj1 <- as.SingleCellExperiment(seurat_obj1)
sce_obj1 <- scDblFinder(sce_obj1)
seurat_obj1$scDblFinder.class <- sce_obj1$scDblFinder.class
seurat_obj1$scDblFinder.score <- sce_obj1$scDblFinder.score
seurat_obj1 <- subset(seurat_obj1, subset = scDblFinder.class == 'singlet')
seurat_obj1[["percent.cp"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^CP-")
seurat_obj1[["percent.mt"]] <- PercentageFeatureSet(seurat_obj1, pattern = "^MT-")
VlnPlot(seurat_obj1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4)
FeatureScatter(seurat_obj1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(seurat_obj1, file = "Lemna_seurat1.rds")
