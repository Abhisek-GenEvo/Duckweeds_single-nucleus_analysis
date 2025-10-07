library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)

multiome_data <- Read10X_h5("filtered_feature_bc_matrix.h5")
-----------------------
#Create a GRanges object
library(GenomicRanges)
library(data.table)
gtf_file_multi <- "Wmic_Reference_Genome_Complete.gtf"
gtf_data_multi <- fread(gtf_file_multi, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(gtf_data_multi) <- c("seqnames", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
extract_attributes <- function(attribute_string, attribute_name) {
  sub(paste0(".*", attribute_name, " \"([^\"]+)\".*"), "\\1", attribute_string)
}
gtf_data_multi[, transcript_id := extract_attributes(attributes, "transcript_id")]
gtf_data_multi[, gene_name := extract_attributes(attributes, "gene_name")]
gtf_data_multi[, gene_id := extract_attributes(attributes, "gene_id")]
gtf_data_multi[, gene_biotype := extract_attributes(attributes, "gene_biotype")]
gr_multi <- GRanges(seqnames = gtf_data_multi$seqnames,
               ranges = IRanges(start = gtf_data_multi$start, end = gtf_data_multi$end),
               strand = gtf_data_multi$strand,
               tx_id = gtf_data_multi$transcript_id,
               gene_name = gtf_data_multi$gene_name,
               gene_id = gtf_data_multi$gene_id,
               gene_biotype = gtf_data_multi$gene_biotype,
               type = gtf_data_multi$type)

-----------------------
seurat_obj <- CreateSeuratObject(counts = multiome_data$`Gene Expression`, assay = "RNA", project = "wolffia_multiome")
seurat_obj[['ATAC']] <- CreateChromatinAssay(counts = multiome_data$`Peaks`, annotation = gr_multi, fragments = "atac_fragments.tsv.gz", sep = c(":", "-"))
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^CP-", col.name = "percent.cp", assay = "RNA")
seurat_obj <- NucleosomeSignal(seurat_obj, assay = "ATAC")
seurat_obj <- TSSEnrichment(seurat_obj, assay = "ATAC")
VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt", "percent.cp", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"), ncol = 6, pt.size = 0)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 1 & percent.cp < 0.5 & nFeature_ATAC > 100 & nFeature_ATAC < 18000 & TSS.enrichment > 1 & nucleosome_signal < 1)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

seurat_obj <- FindNeighbors(seurat_obj,
                        reduction = "umap_rna",
                        dims = 1:ncol(Embeddings(seurat_obj,"umap_rna"))) %>%
                        FindClusters(resolution = 0.5)

DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 50)
seurat_obj <- RunTFIDF(seurat_obj, method = 1)
seurat_obj <- RunSVD(seurat_obj, n = 50)
p1 <- ElbowPlot(seurat_obj, ndims = 30, reduction="lsi")
p2 <- DepthCor(seurat_obj, n = 30)
p1 | p2
seurat_obj <- RunUMAP(seurat_obj,
                  reduction = "lsi",
                  dims = 2:30,
                  reduction.name = "umap_atac",
                  reduction.key = "UMAPATAC_")
seurat_obj <- FindMultiModalNeighbors(seurat_obj,
                                  reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:ncol(Embeddings(seurat_obj,"pca")),
                                                   1:ncol(Embeddings(seurat_obj,"lsi"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)
seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", assay = "RNA")
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", resolution = 0.5)
UMAPPlot(seurat_obj, group.by = "orig.ident") & NoAxes()
UMAPPlot(seurat_obj, group.by = "wsnn_res.0.5", label=T) 

#Marker genes
DefaultAssay(seurat_obj) <- "RNA"
markers_seurat_obj <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




