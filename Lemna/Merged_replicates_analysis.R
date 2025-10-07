library(Seurat)
P1_lemna <- readRDS("/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/upto_scdblFinder_QC/P1_lemna_1.rds")
P2_lemna <- readRDS("/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/upto_scdblFinder_QC/P2_lemna_1.rds")
P1_lemna <- subset(P1_lemna, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 0.25 & percent.cp < 0.15)
print(P1_lemna)
P2_lemna <- subset(P2_lemna, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 0.25 & percent.cp < 0.15)
print(P2_lemna)
P_merged <- merge(P1_lemna,
              y = P2_lemna,
              add.cell.ids = c("P1", "P2"),
              project = "Merged_P")
P_merged <- NormalizeData(P_merged, normalization.method = "LogNormalize", scale.factor = 10000)
P_merged <- FindVariableFeatures(P_merged, selection.method = "vst", nfeatures = 3000)
all.genes.P_merged <- rownames(P_merged)
P_merged <- ScaleData(P_merged, features = all.genes.P_merged)
P_merged <- RunPCA(P_merged, features = VariableFeatures(object = P_merged))
P_merged <- FindNeighbors(P_merged, dims = 1:30, reduction = "pca")
P_merged <- FindClusters(P_merged, resolution = 0.5, cluster.name = "unintegrated_clusters")
P_merged <- RunUMAP(P_merged, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
P_merged <- IntegrateLayers(
  object = P_merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
P_merged[["RNA"]] <- JoinLayers(P_merged[["RNA"]])
P_merged <- FindNeighbors(P_merged, reduction = "integrated.cca", dims = 1:30)
P_merged <- FindClusters(P_merged, resolution = 0.5, cluster.name = "cca_clusters")
P_merged <- RunUMAP(P_merged, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
pdf(file = '/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/after_QC_merging/P1P2_umap_complete_cca1.pdf', width=10, height=10)
DimPlot(P_merged, reduction = "umap.cca", group.by = "orig.ident")
dev.off()
pdf(file = '/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/after_QC_merging/P1P2_umap_complete_cca2.pdf', width=10, height=10)
DimPlot(P_merged, reduction = "umap.cca", label = TRUE)
dev.off()
pdf(file = '/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/after_QC_merging/P1P2_umap_complete_cca3.pdf', width=10, height=10)
DimPlot(P_merged, reduction = "umap.cca", split.by = "orig.ident", label = TRUE)
dev.off()
saveRDS(P_merged, file = "/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/after_QC_merging/Before_markers_P1P2_CCA.rds")
markers_P1P2 <- FindAllMarkers(P_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_P1P2, "/lustre/project/m2_jgu-evoltroph/achakrab/LA_cellranger_analysis/Final_separate_seurat_analysis/after_QC_merging/P1P2_markers_louvian_CCA.csv", row.names = FALSE, quote = FALSE)
