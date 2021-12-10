library(dplyr)
install.packages("Seurat")
library(Seurat)
install.packages("googledrive")
library("googledrive")
install.packages("gprofiler2")
library(gprofiler2)
library(patchwork)
library(ggplot2)
library(ggfortify)

#set working directory to folder containing cardiomyocyte data
setwd("~/Desktop/Capstone")

# Load all cardiomyocyte data set
DOES.data.t <-read.csv("~/Desktop/Capstone/GSM3728971_101894_DOES_merged.genematrix.csv")
d2.meso.t <-read.csv("~/Desktop/Capstone/GSM3728972_101888_d2mesoderm_merged.genematrix.csv")
d5.progenitor.t <-read.csv("~/Desktop/Capstone/GSM3728973_100306_d5cpc_merged.genematrix.csv")
d9.induction.t <-read.csv("~/Desktop/Capstone/GSM3728974_100286_h_icm_merged.genematrix.csv")
d14.cardio.t <-read.csv("~/Desktop/Capstone/GSM3728975_100352_d14CM_merged.genematrix.csv")
d60.mature.t <-read.csv("~/Desktop/Capstone/GSM3728976_105125d60_merged.genematrix.csv")
# Initialize the Seurat object with the raw (non-normalized data).

#transpose data to get genes as rows and cell id tags as columns
DOES.data <- t(DOES.data.t)
d2.mesoderm <-t(d2.meso.t)
d5.progenitor <-t(d5.progenitor.t)
d9.induction <-t(d9.induction.t)
d14.cardiomyocyte <-t(d14.cardio.t)
d60.maturecardio <-t(d60.mature.t)

#CONVERT EMBRYONIC DATA TO HGNC NOTATION
rows <- rownames(DOES.data)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames <-conversion[,4]

#CONVERT MESODERM DATA TO HGNC NOTATION
rows.meso <- rownames(d2.mesoderm)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows.meso, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames.meso <-conversion[,4]

#CONVERT PROGENITOR DATA TO HGNC NOTATION
rows.progenitor <- rownames(d5.progenitor)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows.progenitor, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames.progenitor <-conversion[,4]

#CONVERT INDUCTION DATA TO HGNC NOTATION
rows.induction <- rownames(d9.induction)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows.induction, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames.induction <-conversion[,4]

#CONVERT CARDIO DATA TO HGNC NOTATION
rows.cardio <- rownames(d14.cardiomyocyte)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows.cardio, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames.cardio <-conversion[,4]

#CONVERT MATURE CARDIO DATA TO HGNC NOTATION
rows.mature <- rownames(d60.maturecardio)
#convert genes from ensemble to HGNC notation
conversion <- gconvert(query = rows.mature, organism = "hsapiens", target = "HGNC", mthreshold = 1, filter_na = FALSE)
geneNames.mature <-conversion[,4]

# Create Seurat Object (a representation of single-cell expression data)
does <- CreateSeuratObject(counts = DOES.data, project = "DOES", min.cells = 3, min.features = 200)
mesoderm <- CreateSeuratObject(counts = d2.mesoderm, project = "mesoderm", min.cells = 3, min.features = 200)
progenitor <- CreateSeuratObject(counts = d5.progenitor, project = "progenitor", min.cells = 3, min.features = 200)
induction <- CreateSeuratObject(counts = d9.induction, project = "induction", min.cells = 3, min.features = 200)
cardio <- CreateSeuratObject(counts = d14.cardiomyocyte, project = "cardiomyocyte", min.cells = 3, min.features = 200)
maturecardio <- CreateSeuratObject(counts = d60.maturecardio, project = "mature cardiomyocyte", min.cells = 3, min.features = 200)

# Visualize QC metrics as a violin plot
VlnPlot(does, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(mesoderm, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(progenitor, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(induction, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(cardio, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(maturecardio, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


cardiomyocyte.combined <- merge(does, y = c(mesoderm, progenitor, induction, cardio, maturecardio), add.cell.ids = c("Day 0", "Day 2", "Day 5", "Day 9", "Day 14", "Day 60"), 
                                project = "Cardiomyocyte Transition")
cardiomyocyte.combined

cardio.middle <- merge(mesoderm, y = c(progenitor, induction, cardio), add.cell.ids = c("Day 2", "Day 5", "Day 9", "Day 14"), project = "Cardiomyocyte Middle Transition")

cardio.middle

head(colnames(cardiomyocyte.combined))
## [1] "Day 0_AACCAACGCAA" "Day 0_AACCAACGGTT" "Day 0_AACCAACTAGA" "Day 0_AACCAACTGAC" "Day 0_AACCAAGACGA" "Day 0_AACCAAGATTC"

tail(colnames(cardiomyocyte.combined))
## [1] "Day 60_TTGGTCAAGCT" "Day 60_TTGGTCATAGG" "Day 60_TTGGTCCTACT" "Day 60_TTGGTCGCTAG" "Day 60_TTGGTCTGACG"
## [6] "Day 60_TTGGTTAGATG"

unique(sapply(X = strsplit(colnames(cardiomyocyte.combined), split = "_"), FUN = "[", 1))
## [1] "Day 0"  "Day 2"  "Day 5"  "Day 9"  "Day 14" "Day 60"

table(cardiomyocyte.combined$orig.ident)

does <- NormalizeData(does, normalization.method = "LogNormalize")
mesoderm <-NormalizeData(mesoderm, normalization.method = "LogNormalize")
progenitor <-NormalizeData(progenitor, normalization.method = "LogNormalize")
induction <-NormalizeData(induction, normalization.method = "LogNormalize")
cardio <-NormalizeData(cardio, normalization.method = "LogNormalize")
maturecardio <-NormalizeData(maturecardio, normalization.method = "LogNormalize")

cardiomyocyte.combined <- merge(does, y = c(mesoderm, progenitor, induction, cardio, maturecardio), add.cell.ids = c("Day 0", "Day 2", "Day 5", "Day 9", "Day 14", "Day 60"), project = "Cardiomyocyte",
                         merge.data = TRUE)
GetAssayData(cardiomyocyte.combined)[1:10, 1:15]

cardiomyocyte.combined <- FindVariableFeatures(cardiomyocyte.combined, selection.method = "vst", nfeatures = 2000)

cardio.middle <- merge(mesoderm, y = c(progenitor, induction, cardio), add.cell.ids = c("Day 2", "Day 5", "Day 9", "Day 14"), project = "Cardiomyocyte Middle Transition",
                       merge.data = TRUE)
GetAssayData(cardio.middle)[1:10, 1:15]

cardio.middle <- FindVariableFeatures(cardio.middle, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cardiomyocyte.combined), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cardiomyocyte.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1

plot2

all.genes <- rownames(cardiomyocyte.combined)
some.genes <- rownames(cardio.middle)
cardiomyocyte.combined <- ScaleData(cardiomyocyte.combined, features = all.genes)
cardio.middle <- ScaleData(cardio.middle, features = some.genes)

## Centering and scaling data matrix

cardiomyocyte.combined <- RunPCA(cardiomyocyte.combined, features = VariableFeatures(object = cardiomyocyte.combined), verbose = FALSE)
mat <- Seurat::GetAssayData(cardiomyocyte.combined, assay = "RNA", slot = "scale.data")
pca <- cardiomyocyte.combined[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

#cardiomyocyte.combined <- RunPCA(cardiomyocyte.combined, features = VariableFeatures(object = cardiomyocyte.combined))

# Examine and visualize PCA results a few different ways
print(cardiomyocyte.combined[["pca"]], dims = 1:2, nfeatures = 5)

VizDimLoadings(cardiomyocyte.combined, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(cardiomyocyte.combined, reduction = "pca") + plot_annotation(xlab = eigValues[1], ylab = eigValues[2])

cardio.middle <- RunPCA(cardio.middle, features = VariableFeatures(object = cardio.middle), verbose = FALSE)
mat2 <- Seurat::GetAssayData(cardio.middle, assay = "RNA", slot = "scale.data")
pca2 <- cardio.middle[["pca"]]

# Get the total variance:
total_variance2 <- sum(matrixStats::rowVars(mat2))

eigValues2 = (pca2@stdev)^2  ## EigenValues
varExplained2 = eigValues2 / total_variance2

print(cardio.middle[["pca"]], dims = 1:2, nfeatures = 5)

VizDimLoadings(cardio.middle, dims = 1:2, nfeatures = 15, reduction = "pca")

DimPlot(cardio.middle, reduction = "pca") + plot_annotation(xlab = eigValues2[1], ylab = eigValues2[2])

#pca2<-prcomp(cardio.middle, scale=TRUE)
#plot(pca2$x[,1],pca2$x[,2])

cardiomyocyte.combined <- RunTSNE(object = cardiomyocyte.combined, dims.use = 1:10, do.fast = TRUE, check_duplicates = FALSE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = cardiomyocyte.combined, reduction = "tsne")

VizDimLoadings(cardiomyocyte.combined, dims = 1:2, nfeatures = 15, reduction = "tsne")

DimPlot(cardiomyocyte.combined, reduction = "tsne")

DimHeatmap(cardiomyocyte.combined, dims = 1, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
## APPROX 2 minutes PROCESSING TIME
#cardiomyocyte.combined <- JackStraw(cardiomyocyte.combined, num.replicate = 100)
#cardiomyocyte.combined <- ScoreJackStraw(does, dims = 1:20)

#JackStrawPlot(cardiomyocyte.combined, dims = 1:15)

#ElbowPlot(cardiomyocyte.combined)

cardiomyocyte.combined <- FindNeighbors(cardiomyocyte.combined, dims = 1:10)
## Computing nearest neighbor graph

## Computing SNN
cardiomyocyte.combined <- FindClusters(cardiomyocyte.combined, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(cardiomyocyte.combined), 5)

cardiomyocyte.combined <- RunUMAP(cardiomyocyte.combined, dims = 1:10)

DimPlot(cardiomyocyte.combined, reduction = "umap")

# find all markers of cluster 1
# APPROX RUNNING TIME: 7 MINS
cluster1.markers <- FindMarkers(cardiomyocyte.combined, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

##                        p_val avg_log2FC pct.1 pct.2    p_val_adj
## ENSG00000138768 1.044122e-74 -0.3192285 0.077 0.722 2.904852e-70
## ENSG00000153310 2.039642e-72 -0.4117794 0.072 0.699 5.674487e-68
## ENSG00000259020 1.718568e-71 -0.4291131 0.137 0.835 4.781227e-67
## ENSG00000178913 9.748938e-71 -0.3890494 0.127 0.793 2.712252e-66
## ENSG00000136628 1.120056e-69 -0.3645899 0.147 0.839 3.116109e-65

VlnPlot(cardiomyocyte.combined, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

# find all markers of cluster 2
# APPROX RUNNING TIME: 8 MINS
cluster2.markers <- FindMarkers(cardiomyocyte.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

##                        p_val avg_log2FC pct.1 pct.2    p_val_adj
## ENSG00000250568 1.037847e-43  0.3423989 0.830 0.285 2.887393e-39
## ENSG00000198408 1.089872e-38  0.6377184 0.873 0.369 3.032133e-34
## ENSG00000164754 2.894465e-38  0.6356005 0.939 0.448 8.052692e-34
## ENSG00000213866 3.558106e-38  0.3481164 0.842 0.352 9.899007e-34
## ENSG00000153561 7.972235e-38  0.2608199 0.709 0.210 2.217955e-33

VlnPlot(cardiomyocyte.combined, features = c(row.names(cluster2.markers)[1], row.names(cluster2.markers)[2]))

# find all markers distinguishing cluster 5 from clusters 0 and 3
## cluster5.markers <- FindMarkers(does, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
## head(cluster5.markers, n = 5)
## VlnPlot(does, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2]))

# find markers for every cluster compared to all remaining cells, report only the positive ones
cardio.markers <- FindAllMarkers(cardiomyocyte.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

x <- cardio.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
FeaturePlot(cardiomyocyte.combined, features = x$gene[1:4])
FeaturePlot(cardiomyocyte.combined, features = x$gene[5:8])
FeaturePlot(cardiomyocyte.combined, features = x$gene[9:12])
FeaturePlot(cardiomyocyte.combined, features = x$gene[13:16])
FeaturePlot(cardiomyocyte.combined, features = x$gene[17:18])

p <- FeaturePlot(cardiomyocyte.combined, features = c("ENSG00000197061", "ENSG00000233839", "ENSG00000281181"), combine = FALSE)

p <- lapply(X = p, FUN = function(x) x + 
              theme(plot.title = element_text(size = 8)) +
              theme(axis.title.y = element_text(size = 5)) +
              theme(axis.title.x = element_text(size = 5)) +
              theme(axis.text.y = element_text(size = 5)) +
              theme(axis.text.x = element_text(size = 5)) +
              theme(legend.position = "none")  )

CombinePlots(plots = p)

top10 <- cardio.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

top10

## [1] "ENSG00000186767" "ENSG00000183036" "ENSG00000180354" "ENSG00000254427" "ENSG00000255458" "ENSG00000270164"
## [7] "ENSG00000259658" "ENSG00000206120" "ENSG00000121350" "ENSG00000143786"

p2 <- DoHeatmap(cardiomyocyte.combined, features = top10$gene, group.bar.height = 0.01,size=3,combine = FALSE) 

p2 <- lapply(X = p2, FUN = function(x) x + 
               theme(plot.title = element_text(size = 8)) +
               theme(axis.title.y = element_text(size = 5)) +
               theme(axis.title.x = element_text(size = 5)) +
               theme(axis.text.y = element_text(size = 3)) +
               theme(legend.position = "left")  )

CombinePlots(plots = p2)

new.cluster.ids <- c("0", "1", "2")
names(new.cluster.ids) <- levels(cardiomyocyte.combined)

does <- RenameIdents(cardiomyocyte.combined, new.cluster.ids)
DimPlot(cardiomyocyte.combined, reduction = "pca", label = TRUE, pt.size = 0.5)

cardiomyocyte.combined

DimPlot(cardiomyocyte.combined, reduction = "umap", label = TRUE, pt.size = 0.5)

mergedData <- as.matrix(Seurat::GetAssayData(cardiomyocyte.combined, slot = "data"))
write.table(mergedData,row.names = TRUE,col.names = TRUE, sep = "\t", file = file.path("merged_cardio_data","combined.tsv"),quote = FALSE)

mergedmiddle <- as.matrix(Seurat::GetAssayData(cardio.middle, slot = "data"))
write.table(mergedmiddle,row.names = TRUE,col.names = TRUE, sep = "\t", file = file.path("merged_cardio_data", "middle.tsv"),quote = FALSE)

