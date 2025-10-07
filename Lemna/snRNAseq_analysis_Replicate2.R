library(Seurat)
P2_rna <- Read10X_h5("filtered_feature_bc_matrix_P2.h5")
seurat_obj2 <- CreateSeuratObject(counts = P2_rna, project = "seurat_obj2")
seurat_obj2 <- subset(seurat_obj2, subset = nCount_RNA > 500)

library(scDblFinder)
library(SingleCellExperiment)
set.seed(123)
sce_P2 <- as.SingleCellExperiment(P2_lemna)
sce_P2 <- scDblFinder(sce_P2)
P2_lemna$scDblFinder.class <- sce_P2$scDblFinder.class
P2_lemna$scDblFinder.score <- sce_P2$scDblFinder.score
P2_lemna <- subset(P2_lemna, subset = scDblFinder.class == 'singlet')
P2_lemna[["percent.cp"]] <- PercentageFeatureSet(P2_lemna, pattern = "^CP-")
P2_lemna[["percent.mt"]] <- PercentageFeatureSet(P2_lemna, pattern = "^MT-")
VlnPlot(P2_lemna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.cp"), ncol = 4)
FeatureScatter(P2_lemna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(P2_lemna, file = "/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/upto_scdblFinder_QC/P2_lemna_1.rds")
