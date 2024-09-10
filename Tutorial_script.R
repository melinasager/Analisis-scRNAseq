install.packages("Seurat")

library(Seurat)


counts <- Read10X(data.dir = "C:/Users/Usuario/Desktop/meli/r/Tutorial_github_scrnaseq/DS1/")
seurat <- CreateSeuratObject(counts, project="DS1")

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)


#Opcion 1: Normalizacion log
seurat <- NormalizeData(seurat)

seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

seurat <- ScaleData(seurat)
## si quiero sacar fuentes de variación ####seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))



#Opcion 2: usar SCTransform
seurat <- SCTransform(seurat, variable.features.n = 3000)

## si quiero sacar fuentes de variación ## seurat <- SCTransform(seurat,
#                      vars.to.regress = c("nFeature_RNA", "percent.mt"),
#                     variable.features.n = 3000)


seurat <- RunPCA(seurat, npcs = 50)

ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

PCHeatmap(seurat, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)


seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

plot1 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "umap")
plot1 / plot2


seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 1)

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2


ct_markers <- c("MKI67","NES","DCX","FOXG1", # G2M, NPC, neuron, telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon related
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related
DoHeatmap(seurat, features = ct_markers) + NoLegend()


cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##### Si tarda mucho el FindAllMarkers se puede usar Presto.
###install.packages("presto")
###install.packages("DESeq2")

###library(BiocManager)
###library(DESeq2)
###library(presto)

######cl_markers_presto <- wilcoxauc(seurat)
#cl_markers_presto %>%
#  filter(logFC > log(1.2) & pct_in > 20 & padj < 0.05) %>%
#  group_by(group) %>%
#  arrange(desc(logFC), .by_group=T) %>%
#  top_n(n = 2, wt = logFC) %>%
#  print(n = 40, width = Inf)

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()

plot1 <- FeaturePlot(seurat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(seurat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1 | plot2 + plot_layout(widths = c(1,1))
combined_plot <- wrap_plots(plot1, plot2, widths = c(0.75, 2))
combined_plot

plot3 <- FeaturePlot(seurat, c("CLU"), ncol = 1)
plot4 <- VlnPlot(seurat, features = c("CLU"), pt.size = 0)
plot3 | plot4




