library(reticulate)
library(Seurat)
library(cowplot)
library(ggplot2)


#Load 10X datasetscreate their respective Seurat objects
RD.data <- Read10X(data.dir = "/Users/Pang-Kuo/Desktop/RD_HFD_Seurat_analysis/RD_236/", gene.column = 2)
HFD.data <- Read10X(data.dir = "/Users/Pang-Kuo/Desktop/RD_HFD_Seurat_analysis/HFD_243/", gene.column = 2)


#Create their respective Seurat objects
RD <- CreateSeuratObject(counts = RD.data, project = "RD", min.cells = 5)
RD$HFD <- "RD"
RD <- subset(RD, subset = nFeature_RNA > 500)
RD <- NormalizeData(RD, verbose = FALSE)
RD <- FindVariableFeatures(RD, selection.method = "vst", nfeatures = 2000)


HFD <- CreateSeuratObject(counts = HFD.data, project = "HFD", min.cells = 5)
HFD$HFD <- "HFD"
HFD <- subset(HFD, subset = nFeature_RNA > 500)
HFD <- NormalizeData(HFD, verbose = FALSE)
HFD <- FindVariableFeatures(HFD, selection.method = "vst", nfeatures = 2000)


#Perform integration
RD_HFD.anchors <- FindIntegrationAnchors(object.list = list(RD, HFD), dims = 1:15)
RD_HFD.combined <- IntegrateData(anchorset = RD_HFD.anchors, dims = 1:15)


#Perform an integrated analysis
DefaultAssay(RD_HFD.combined) <- "integrated"
RD_HFD.combined <- ScaleData(RD_HFD.combined, verbose = FALSE)
RD_HFD.combined <- RunPCA(RD_HFD.combined, npcs = 30, verbose = FALSE)


#Import python to use UMAP
use_virtualenv("base")
use_python("/users/Pang-Kuo/Anaconda3/python.exe")


#UMAP and Clustering
RD_HFD.combined <- FindNeighbors(RD_HFD.combined, reduction = "pca", dims = 1:15)
RD_HFD.combined <- FindClusters(RD_HFD.combined, resolution = 0.5)
RD_HFD.combined <- RunUMAP(RD_HFD.combined, reduction = "pca", dims = 1:15)


#Save integrated Seurat object as a RDS file for the future use
setwd("C:/Users/Pang-Kuo/Desktop/RD_HFD_Seurat_analysis/RD-HFD_integration")
saveRDS(RD_HFD.combined, file = "../RD-HFD_integration/RD-HFD_combined.rds")


#Visualization
p1 <- DimPlot(RD_HFD.combined, reduction = "umap", group.by = "HFD")
p2 <- DimPlot(RD_HFD.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


#Split visualization
DimPlot(RD_HFD.combined, reduction = "umap", split.by = "HFD")


#Identify conserved cell type markers
DefaultAssay(RD_HFD.combined) <- "RNA"
c0.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 0, grouping.var = "HFD", verbose = FALSE)
head(c0.markers)

c1.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 1, grouping.var = "HFD", verbose = FALSE)
head(c1.markers)

c2.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 2, grouping.var = "HFD", verbose = FALSE)
head(c2.markers)

c3.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 3, grouping.var = "HFD", verbose = FALSE)
head(c3.markers)

c4.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 4, grouping.var = "HFD", verbose = FALSE)
head(c4.markers)

c5.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 5, grouping.var = "HFD", verbose = FALSE)
head(c5.markers)

c6.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 6, grouping.var = "HFD", verbose = FALSE)
head(c6.markers)

c7.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 7, grouping.var = "HFD", verbose = FALSE)
head(c7.markers)

c8.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 8, grouping.var = "HFD", verbose = FALSE)
head(c8.markers)

c9.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 9, grouping.var = "HFD", verbose = FALSE)
head(c9.markers)

c10.markers <- FindConservedMarkers(RD_HFD.combined, ident.1 = 10, grouping.var = "HFD", verbose = FALSE)
head(c10.markers)



#Visualization of cell type marker genes
FeaturePlot(RD_HFD.combined, features = c("Vim", "Sema3c", "Hspa1a", "Col4a1", "Cd14", "Rgcc", "Cxcl1", 
                                          "Nkg7", "Cd3e", "Cd79a"), min.cutoff = "q9")


#Classify cell types in UMAP plots
RD_HFD.combined <- RenameIdents(RD_HFD.combined, `0` = "Stroma-C1", `1` = "Stroma-C2", `2` = "Stroma-C3", 
                                `3` = "Cd14+ Mono", `4` = "Stroma-C4", `5` = "Stroma-C5", `6` = "NK", `7` = "T cells", `8` = "B cells", `9` = "C9", 
                                `10` = "C10")

DimPlot(RD_HFD.combined, split.by = "HFD", label = TRUE)


#Visualization of the relationship between cell type marker genes and cell types
Idents(RD_HFD.combined) <- factor(Idents(RD_HFD.combined), levels = c("Stroma-C1", "Stroma-C2", "Stroma-C3", "Cd14+ Mono", "Stroma-C4", "Stroma-C5", "NK", "T cells", "B cells", "C9", "C10")) 
                                                                      
markers.to.plot <- c("Sema3c", "Hspa1a", "Col4a1", "Cd14", "Rgcc", "Cxcl1", "Nkg7", "Cd3e", "Cd79a", "Cldn4", "Ccr9")
                     
DotPlot(RD_HFD.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "HFD") + RotatedAxis()
       

#Identify differential expressed genes across conditions
##Create scatter plots for expressed genes accross conditions
Stroma.C1 <- subset(RD_HFD.combined, idents = "Stroma-C1")
Idents(Stroma.C1) <- "HFD"
avg.Stroma.C1 <- log1p(AverageExpression(Stroma.C1, verbose = FALSE)$RNA)
avg.Stroma.C1$gene <- rownames(avg.Stroma.C1)

Stroma.C2 <- subset(RD_HFD.combined, idents = "Stroma-C2")
Idents(Stroma.C2) <- "HFD"
avg.Stroma.C2 <- log1p(AverageExpression(Stroma.C2, verbose = FALSE)$RNA)
avg.Stroma.C2$gene <- rownames(avg.Stroma.C2)

Stroma.C3 <- subset(RD_HFD.combined, idents = "Stroma-C3")
Idents(Stroma.C3) <- "HFD"
avg.Stroma.C3 <- log1p(AverageExpression(Stroma.C3, verbose = FALSE)$RNA)
avg.Stroma.C3$gene <- rownames(avg.Stroma.C3)

Stroma.C4 <- subset(RD_HFD.combined, idents = "Stroma-C4")
Idents(Stroma.C4) <- "HFD"
avg.Stroma.C4 <- log1p(AverageExpression(Stroma.C4, verbose = FALSE)$RNA)
avg.Stroma.C4$gene <- rownames(avg.Stroma.C4)

Stroma.C5 <- subset(RD_HFD.combined, idents = "Stroma-C5")
Idents(Stroma.C5) <- "HFD"
avg.Stroma.C5 <- log1p(AverageExpression(Stroma.C5, verbose = FALSE)$RNA)
avg.Stroma.C5$gene <- rownames(avg.Stroma.C5)


p1 <- ggplot(avg.Stroma.C1, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 1") + geom_jitter() + geom_smooth(method = 'lm')
p2 <- ggplot(avg.Stroma.C2, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 2") + geom_jitter() + geom_smooth(method = 'lm')
p3 <- ggplot(avg.Stroma.C3, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 3") + geom_jitter() + geom_smooth(method = 'lm')
p4 <- ggplot(avg.Stroma.C4, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 4") + geom_jitter() + geom_smooth(method = 'lm')
p5 <- ggplot(avg.Stroma.C5, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 5") + geom_jitter() + geom_smooth(method = 'lm')
plot_grid(p1, p2, p3, p4, p5)


##Identify differentially expressed genes in a specific cell type across conditions
RD_HFD.combined$celltype.HFD <- paste(Idents(RD_HFD.combined), RD_HFD.combined$HFD, sep = "_")
RD_HFD.combined$celltype <- Idents(RD_HFD.combined)
Idents(RD_HFD.combined) <- "celltype.HFD"
Stroma__C5.HFD <- FindMarkers(RD_HFD.combined, ident.1 = "Stroma-C5_HFD", ident.2 = "Stroma-C5_RD", verbose = FALSE)
Stroma__C5.HFD.sorted <- Stroma__C5.HFD[order(Stroma__C5.HFD$avg_logFC, decreasing = TRUE),]
head(Stroma__C5.HFD.sorted, 8)
tail(Stroma__C5.HFD.sorted, 8)

##Visualization of top differentially expressed genes in scatter plots
genes.to.label = c("Col3a1", "Igfbp6", "Sparc", "Eln", "Meg3", "Inmt", "Iigp1", "Ccl2", "Cxcl13", "Lif", "Ier3", "Ccl7")
plot5 <- ggplot(avg.Stroma.C5, aes(RD, HFD)) + geom_point() + ggtitle("Stromal cell cluster 5") + geom_jitter() + geom_smooth(method = 'lm')
plot5 <- LabelPoints(plot = plot5, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
plot5


##Visualization of top differentially expressed genes in UMAP plots
FeaturePlot(RD_HFD.combined, features = c("Col3a1", "Igfbp6"), split.by = "HFD", max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(RD_HFD.combined, features = c("Iigp1", "Ccl2"), split.by = "HFD", max.cutoff = 3, cols = c("grey", "red"))         


##Visualization of top differentially expressed genes in Vlnplots
plots <- VlnPlot(RD_HFD.combined, features = c("Col3a1", "Igfbp6", "Iigp1", "Ccl2"), split.by = "HFD", group.by = "celltype", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)






