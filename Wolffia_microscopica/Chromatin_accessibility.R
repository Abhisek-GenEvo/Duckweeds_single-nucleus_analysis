#Linking marker genes to nearby peaks
library(Biostrings)
library(GenomicRanges)
library(data.table)
custom_genome_fasta <- "Wmic_Reference_Genome_Complete.fasta"
custom_genome <- readDNAStringSet(custom_genome_fasta, format = "fasta")
seurat_obj <- RegionStats(seurat_obj, genome = custom_genome)
seurat_obj <- LinkPeaks(
  object = seurat_obj,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = top_marker_genes$feature
)
-------------------------
#Chromatin accessibility profiles
new_cluster_ids <- c("Cluster 0", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6",  "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Cluster 14", "Cluster 15")
names(new_cluster_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)
seurat_obj$celltype <- Idents(seurat_obj)
print(head(seurat_obj$celltype))
CoveragePlot(seurat_obj,
                   region = "Wolffia04G002500",
                   features = "Wolffia04G002500",
                   group.by = "celltype",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
--------------------------
#TF binding motif enrichment analysis
library(TFBSTools)
library(JASPAR2020)
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'plants', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

library(Biostrings)
custom_genome <- readDNAStringSet("Wmic_Reference_Genome_Complete.fasta")
DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- AddMotifs(seurat_obj, genome = custom_genome, pfm = pfm)

open_peaks <- AccessiblePeaks(seurat_obj)
peaks_matched <- MatchRegionStats(meta.feature = seurat_obj[['ATAC']]@meta.features[open_peaks, ],
                                  query.feature = seurat_obj[['ATAC']]@meta.features[top_peaks$feature, ], n = 38542)

motif_enrichment_endo <- FindMotifs(seurat_obj,
                                    features = top_peaks$feature[top_peaks$group == "15"],
                                    background = peaks_matched) %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol), df_pfm$id)[motif]) %>%
  mutate(padj = p.adjust(pvalue, method="BH"))
enriched_motif_endo <- motif_enrichment_endo %>%
  dplyr::filter(padj < 0.01 & fold.enrichment > 0.25) %>%
  top_n(-4, wt = padj)
print(head(enriched_motif_endo))
p1 <- MotifPlot(seurat_obj, motifs = enriched_motif_endo$motif[1:2], ncol=2)
p1
