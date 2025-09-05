## single cell data analysis by Seurat V3
## Author: Tamrin Chowdhury
## Date: 3/25/2020

library(dplyr)
library(Seurat)

#import single cell data counts
scflu. <- read.delim("E:/R projects/190508/scflu..txt", row.names=1)

#set unique rownames if shows duplicate rownames error
nams = df$Gene
rownames(df) = make.names(nams, unique=TRUE)
df$Gene <- NULL

#Create seuratobject
sct <- CreateSeuratObject(counts = scflu., project = "sct1k")
sct

# Visualize QC metrics as a violin plot
VlnPlot(sct, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#view scatter plot for data filtering
FeatureScatter(sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") #not run

#subset accordingly if necessary
sct <- subset(sct, subset = nCount_RNA > 200 & nFeature_RNA < 2500)

#Normalize data
sct <- NormalizeData(sct)

#Find variables
sct <- FindVariableFeatures(sct, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sct), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sct)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#scale data
all.genes <- rownames(sct)
sct <- ScaleData(sct, features = all.genes)

#Dimensional reduction
sct <- RunPCA(sct, features = VariableFeatures(object = sct))
print(sct[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sct, dims = 1:2, reduction = "pca")
DimPlot(sct, reduction = "pca")
DimHeatmap(sct, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sct, dims = 1:15, cells = 500, balanced = TRUE)

#Determine demensionality
sct <- JackStraw(sct, num.replicate = 100)
sct <- ScoreJackStraw(sct, dims = 1:20)
JackStrawPlot(sct, dims = 1:15)

ElbowPlot(sct)

#cluster cells
sct <- FindNeighbors(sct, dims = 1:10)
sct <- FindClusters(sct, resolution = 0.5)

head(Idents(sct), 5)

#Non-linear dimesional reduction
sct <- RunUMAP(sct, dims = 1:10)
DimPlot(sct, reduction = "umap")

sct <- RunTSNE(sct, dims = 1:25, method = "FIt-SNE")
DimPlot(sct, reduction = "tsne", label = TRUE) + NoLegend()

#extracting coordinates of tsne, umap, pca
tplotco <- Embeddings(object = sct[["tsne"]])
uplotco <- Embeddings(object = sct[["umap"]])
pplotco <- Embeddings(object = sct[["pca"]])


#Find top20 markers in every cluster
sct.markers <- FindAllMarkers(sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- sct.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

#Number of cells in each cluster
cell.num <- table(Idents(sct))

#Write the tables for top 20 markers and cell numbers 
write.table(top20, file = "top20_genepercluster.txt", sep = "\t")
write.table(cell.num, file = "cellnumberpercluster.txt", sep = "\t")

#save the object for later use
saveRDS(sct, file = "seurat_test.rds")


