library(Seurat)
seurat_obj1 <- readRDS("Lemna_seurat1.rds")
seurat_obj2 <- readRDS("Lemna_seurat2.rds")
seurat_obj1 <- subset(seurat_obj1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 0.25 & percent.cp < 0.15)
seurat_obj2 <- subset(seurat_obj2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 10000 & percent.mt < 0.25 & percent.cp < 0.15)
obj_merged <- merge(seurat_obj1,
              y = seurat_obj2,
              add.cell.ids = c("P1", "P2"),
              project = "Merged_P")
obj_merged <- NormalizeData(obj_merged, normalization.method = "LogNormalize", scale.factor = 10000)
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst", nfeatures = 3000)
all.genes.obj_merged <- rownames(obj_merged)
obj_merged <- ScaleData(obj_merged, features = all.genes.obj_merged)
obj_merged <- RunPCA(obj_merged, features = VariableFeatures(object = obj_merged))
obj_merged <- FindNeighbors(obj_merged, dims = 1:30, reduction = "pca")
obj_merged <- FindClusters(obj_merged, resolution = 0.5, cluster.name = "unintegrated_clusters")
obj_merged <- RunUMAP(obj_merged, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated")
obj_merged <- IntegrateLayers(
  object = obj_merged,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
obj_merged[["RNA"]] <- JoinLayers(obj_merged[["RNA"]])
obj_merged <- FindNeighbors(obj_merged, reduction = "integrated.cca", dims = 1:30)
obj_merged <- FindClusters(obj_merged, resolution = 0.5, cluster.name = "cca_clusters")
obj_merged <- RunUMAP(obj_merged, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(obj_merged, reduction = "umap.cca", group.by = "orig.ident")
DimPlot(P_merged, reduction = "umap.cca", label = TRUE)
DimPlot(P_merged, reduction = "umap.cca", split.by = "orig.ident", label = TRUE)
markers <- FindAllMarkers(obj_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "all_markers.csv", row.names = FALSE, quote = FALSE)
saveRDS(obj_merged, file = "Merged_seurat_object.rds")
