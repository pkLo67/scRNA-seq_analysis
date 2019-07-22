library(reticulate)
library(Seurat)
library(dplyr)


#Set working directory
setwd("C:/Users/Pang-Kuo/Desktop/RD_HFD_Seurat_analysis/HFD_243")


#Read10X() only reads zipped files from CellRanger3.0.
HFD.data <- Read10X(data.dir = "../HFD_243/", gene.column = 2)

#Create Seurat object
HFD <- CreateSeuratObject(counts = HFD.data, project = "HFD", min.cells = 3, min.features = 200)
HFD

#Calculate the percentage of detected mitochondrial genes
HFD[["percent.mt"]] <- PercentageFeatureSet(HFD, pattern = "^MT-")
VlnPlot(HFD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(HFD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HFD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


#Set the optimal cell subset for the following analysis
HFD <- subset(HFD, subset = nFeature_RNA > 200 & nFeature_RNA < 3300 & percent.mt < 5)
VlnPlot(HFD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Perform normalization and identify variable expressed genes
HFD <- NormalizeData(HFD, normalization.method = "LogNormalize", scale.factor = 10000)
HFD <- FindVariableFeatures(HFD, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(HFD), 10)
top10
plot1 <- VariableFeaturePlot(HFD)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1
plot2


#Linear scaling and run PCA analysis
HFD <- ScaleData(HFD)
HFD <- RunPCA(HFD, features = VariableFeatures(object = HFD))
VizDimLoadings(HFD, dims = 1:2, reduction = "pca")
DimPlot(HFD, reduction = "pca")
DimHeatmap(HFD, dims = 1:6, cells = 500, balanced = TRUE)
ElbowPlot(HFD)


#Perform neighboring and clustering analysis
HFD <- FindNeighbors(HFD, dims = 1:15)
HFD <- FindClusters(HFD, resolution = 0.5)
head(Idents(HFD), 5)

#Through reticulate, import python function to run UMAP
use_virtualenv("base")
use_python("/users/Pang-Kuo/Anaconda3/python.exe")
HFD <- RunUMAP(HFD, dims = 1:15)
DimPlot(HFD, reduction = "umap")


#Save Seurat object for the future analysis
saveRDS(HFD, file = "../HFD_243/HFD_243.rds")


#Identify marker genes for clusters
cluster1.markers <- FindMarkers(HFD, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(HFD, ident.1 = 5, ident.2 = c(0, 2), min.pct = 0.25)
head(cluster5.markers, n = 5)


#Identify marker genes for every cell cluster
HFD.markers <- FindAllMarkers(HFD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HFD.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(HFD, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(HFD, features = c("Vim", "Rgcc", "Jun", "Hspa1a", "Igfbp5", "Cd14", "Gzmb", "Cd79a", "Cd3e"), slot = "counts", log = TRUE)
FeaturePlot(HFD, features = c("Vim", "Rgcc", "Jun", "Hspa1a", "Igfbp5", "Cd14", "Gzmb", "Cd79a", "Cd3e"))
top10 <- HFD.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(HFD, features = top10$gene) + NoLegend()


#Replace cluster numbers with cell type names
new.cluster.ids <- c("Stromal_cluster1", "Stromal_cluster2", "Stromal_cluster3", "Cd14+ Mono", "NK", "B", "Cd3+ T", "C7", "C8")
names(new.cluster.ids) <- levels(HFD)
HFD <- RenameIdents(HFD, new.cluster.ids)
DimPlot(HFD, reduction = "umap", label = TRUE, pt.size = 1.0) + NoLegend()


#Make an object including gene names and average expression values for two cell clusters
Stromal.cluster1 <- subset(HFD, idents = "Stromal_cluster1")
avg.s.c1 <- log1p(AverageExpression(Stromal.cluster1, verbose = FALSE)$RNA)
head(avg.s.c1, 5)
SF.c1_c2 <- avg.s.c1

Stromal.cluster2 <- subset(HFD, idents = "Stromal_cluster2")
avg.s.c2 <- log1p(AverageExpression(Stromal.cluster2, verbose = FALSE)$RNA)
head(avg.s.c2, 5)
SF.c1_c2$Stromal_cluster2 <- avg.s.c2[, "Stromal_cluster2"]
head(SF.c1_c2, 5)



#ggplot of Stromal_cluster1 vs Stromal_cluster2
library(ggplot2)
ggplot(SF.c1_c2, aes(Stromal_cluster1, Stromal_cluster2)) + geom_point() + xlab("Stromal_cluster1") + ylab("Stromal_cluster2") + ggtitle("Stromal cluster1 vs cluster2") + geom_jitter() + geom_smooth(method = 'lm')


#Identify differentially expressed genes between Stromal_cluster1 and Stromal_cluster2
cluster1.markers <- FindMarkers(HFD, ident.1 = "Stromal_cluster1", ident.2 = "Stromal_cluster2", min.pct = 0.25, logfc.threshold = 0.25)
cluster1.markers.sorted <- cluster1.markers[order(cluster1.markers$avg_logFC),]
head(cluster1.markers.sorted, 8)
tail(cluster1.markers.sorted, 8)


#Label the top 6 upregulated and top 6 downregulated genes with gene names in the ggplot
genes.to.label = c("Rgcc", "Il1rl1", "Postn", "Mmp3", "Srxn1", "Sdc1", "Tpm1", "Atf3", "Fos", "Jun", "Klf2", "Txnip")
plot1 <- ggplot(SF.c1_c2, aes(x=Stromal_cluster1, y=Stromal_cluster2)) + geom_point() + geom_jitter() + geom_smooth(method = 'lm') + labs(x = "Stromal_cluster1", y = "Stromal_cluster2") + ggtitle("Stromal cluster1 vs cluster2")
plot1 <- LabelPoints(plot = plot1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
plot1


#Fetch the expression and statistical information of a specific gene
View(cluster1.markers["Rgcc",])
Rgcc <- cluster1.markers["Rgcc",]
Rgcc


HFD.markers[c("Fgf7", "Vim"),]





