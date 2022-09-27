
################################################################################
# SPINAL CORD DATASET EXPLORATION
################################################################################

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(gplots)
library(ggplot2)
library(Polychrome)

# inside system
setwd('D:/Presh/h5d_dataset/')
# inside drive
setwd('Y:/PD/spinal_cord/')
# on the outside drive
setwd('U:/eng_research_economo/PD/spinal_cord')
# Load in custom functions:
source('levine_new/scripts/functions.R')




sc_object<-readRDS('levine_new/final_meta_dataset.rds')

# this gives us the first integration plot in the paper
DefaultAssay(sc_object)<-'integrated'
UMAPPlot(sc_object) 


split.list<-SplitObject(sc_object, split.by = "dataset")

# normalize and identify variable features for each dataset independently
split.list <- lapply(X = split.list, FUN = function(x) {
  DefaultAssay(x)<-'raw'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
})

features <- SelectIntegrationFeatures(object.list = split.list, nfeatures = 4000)

split.list <- lapply(X = split.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = T)
  x <- RunPCA(x, features = features, verbose = T)
})


rpca.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features, reduction ="rpca")
cca.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features, reduction ="cca")
saveRDS(rpca.anchors,'levine_new/rpca_anchors.rds')
saveRDS(cca.anchors,'levine_new/cca_anchors.rds')


combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"


combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)

ElbowPlot(combined, ndims = 50)

combined <- RunUMAP(combined, reduction = "pca", dims = 1:40)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:40)
combined <- FindClusters(combined, resolution = 5)

p1 <- DimPlot(combined, reduction = "umap", repel = TRUE,shuffle = TRUE, label = TRUE) 
p1