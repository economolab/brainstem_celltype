#########################################################################
## INTEGRATING BRAINSTEM DATA WITH THE LEVINE DATASET
#########################################################################


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


#### performing integration ####

# loading in the saved file
sc_object <- readRDS('levine_new/levine_dataset/integrate2.rds')
# set Default assay as RNA
DefaultAssay(sc_object)<-'integrated'
DimPlot(sc_object, reduction = "umap", repel = TRUE,shuffle = TRUE, label = TRUE, group.by = 'final_cluster_assignment')  + NoLegend()


bs_object <- readRDS('raw_bs.rds')

bs_object[['percent.mt']] <- PercentageFeatureSet(bs_object, pattern = '^mt-') # You don't need to do this for nuclei
sc_object[['percent.mt']] <- PercentageFeatureSet(sc_object, pattern = '^mt-')


plot1<-VlnPlot(bs_object, features = c('nFeature_RNA', 'nCount_RNA'), ncol =2)

plot2<-VlnPlot(sc_object, features = c('nFeature_RNA', 'nCount_RNA'), ncol =2)

plot1
plot2

plot1 <- FeatureScatter(bs_object, 
                        feature1 = 'nCount_RNA', 
                        feature2 = 'nFeature_RNA')


plot2 <- FeatureScatter(sc_object, 
                        feature1 = 'nCount_RNA', 
                        feature2 = 'nFeature_RNA')
plot1+plot2


bs_object <- subset(bs_object,
                    subset = nCount_RNA > 200
                    #& nFeature_RNA < 5000
                    & percent.mt < 5)

sc_object <- subset(sc_object,
                    subset = nCount_RNA > 200
                    #& nFeature_RNA < 5000
                    & percent.mt < 5)

### For some reason, this ain't working 

# obj.list <- list()
# obj.list[["brainstem"]] <- bs_object
# obj.list[["spinal_cord"]] <- sc_object
# 
# for (i in obj.list) {
#   obj.list[[i]] <- NormalizeData(i, verbose = T)
#   obj.list[[i]] <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000, verbose = T)
# }


bs_object<-NormalizeData(bs_object)
bs_object<-FindVariableFeatures(bs_object, selection.method = "vst", nfeatures = 2000)
# 
# sc_object<-NormalizeData(sc_object)
# sc_object<-FindVariableFeatures(sc_object, selection.method = "vst", nfeatures = 2000)


# plotting the most variable genes 
# Identify the 20 most highly variable genes
top20.bs <- head(x = VariableFeatures(object = bs_object), 
              n =20)

top20.sc <- head(x = VariableFeatures(object = sc_object), 
              n =20)

# Plot variable features with labels
plot1 <- VariableFeaturePlot(object = bs_object, assay = 'RNA')

p1<- LabelPoints(plot = plot1, 
            points = top20.bs, 
            repel = TRUE)


plot2 <- VariableFeaturePlot(object = sc_object, assay = 'RNA')

p2<- LabelPoints(plot = plot2, 
            points = top20.sc, 
            repel = TRUE)


p1+p2


obj.list <- list()
obj.list[["brainstem"]] <- bs_object
obj.list[["spinal_cord"]] <- sc_object

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)


obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction ="rpca")

combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"

# check if the sum of all the NA values in the 2 columns is the same number of cells in the brainstem dataset 


combined@meta.data$final_coarse_types[is.na(combined@meta.data$final_coarse_types)]<-"Brainstem Cells"
combined@meta.data$final_cluster_assignment[is.na(combined@meta.data$final_cluster_assignment)]<-"Brainstem Cells"

metadata.df<-as.data.frame(combined@meta.data)

metadata.df<-metadata.df %>%
  mutate(old.ident = case_when(
    endsWith(orig.ident, "sc_brainstem") ~ "Brainstem Cells",
  ))
metadata.df$old.ident[is.na(metadata.df$old.ident)]<-"Spinalcord Cells"

combined@meta.data<-metadata.df


saveRDS(combined, file = 'levine_new/levine_integrate1.rds')

############# Run the standard workflow for visualization and clustering #################

combined<- readRDS('levine_new/levine_integrate1.rds')


combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, verbose = FALSE)

ElbowPlot(combined, ndims = 50)

combined <- RunUMAP(combined, reduction = "pca", dims = 1:40)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:40)
combined <- FindClusters(combined, resolution = 5)

# Visualization

p1 <- DimPlot(combined, reduction = "umap", repel = TRUE,shuffle = TRUE, label = TRUE) 
p1



saveRDS(combined, file = 'levine_new/levine_integrate1_clustered.rds')

#######################################

combined<-readRDS('levine_new/levine_integrate1_clustered.rds')
# DefaultAssay(combined) <- "RNA"

options(ggrepel.max.overlaps = 10)
colors.needed<-CustomColor(combined,cluster.by='seurat_clusters', seed = T)
p1<-DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'seurat_clusters', cols = colors.needed) + NoLegend()
colors.needed<-CustomColor(combined,cluster.by='final_cluster_assignment', seed = T)
p2<-DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment', cols = colors.needed) + NoLegend()
p1+p2


# filter out dataset below 10% mix
# stacked bar chart to understand the proportions of cells in each cluster
cell_type <- combined@meta.data$old.ident
cluster <- combined@meta.data$seurat_clusters
data <- data.frame(cell_type, cluster)
data <-data %>% group_by_all %>% count

stacked_plot<-ggplot(data, aes(fill=cell_type, y=n, x=cluster)) +
  geom_bar(position="fill", stat="identity")+ geom_hline(yintercept=0.1)+ geom_hline(yintercept=0.9)
stacked_plot


stacked_plot.data<-stacked_plot$data %>%
  group_by(cluster) %>% mutate(Total = sum(n))


stacked_plot.data$percent <- stacked_plot.data$n/stacked_plot.data$Total
stacked_plot.data$percent10<-stacked_plot.data$percent > 0.1
unchosen.clusters<-unique(stacked_plot.data$cluster[stacked_plot.data$percent10==FALSE | stacked_plot.data$percent==1])
unchosen.clusters<-sort(unchosen.clusters)
unchosen.clusters
length(unchosen.clusters)

refined<-subset(combined, idents = unchosen.clusters, invert = TRUE)

saveRDS(refined, file = 'levine_new/levine_integrate1_clustered_thresh10.rds')

refined<-readRDS('levine_new/levine_integrate1_clustered_thresh10.rds')
colors.needed<-CustomColor(refined,cluster.by='seurat_clusters', seed = F)
p1<-DimPlot(refined, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'seurat_clusters', cols = colors.needed) + NoLegend()
colors.needed<-CustomColor(refined,cluster.by='final_cluster_assignment', seed = F)
p2<-DimPlot(refined, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment', cols = colors.needed) + NoLegend()
p1+p2


















