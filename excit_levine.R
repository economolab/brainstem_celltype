#########################################################################
## EXPLORING THE EXCIT DATASET 
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
source('levine_new/scripts/brainstem_celltype/functions.R')



#### performing integration ####

# loading in the saved file

sc_object <- readRDS('levine_new/excit_sc.rds')
DefaultAssay(sc_object)<-'raw'
bs_object <- readRDS('levine_new/excit_bs.rds')



plot1<-VlnPlot(bs_object, features = c('nFeature_RNA', 'nCount_RNA'), ncol =2)
plot1
plot2<-VlnPlot(sc_object, features = c('nFeature_RNA', 'nCount_RNA'), ncol =2)
plot2
plot1+plot2


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


bs_object<-NormalizeData(bs_object)
bs_object<-FindVariableFeatures(bs_object, selection.method = "vst", nfeatures = 2000)


obj.list <- list()
obj.list[["brainstem"]] <- bs_object
obj.list[["spinal_cord"]] <- sc_object


# sc_object<-NormalizeData(sc_object)
# sc_object<-FindVariableFeatures(sc_object, selection.method = "vst", nfeatures = 2000)



features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction ="rpca")


combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"



combined@meta.data$final_coarse_types[is.na(combined@meta.data$final_coarse_types)]<-"Brainstem Cells"
combined@meta.data$final_cluster_assignment[is.na(combined@meta.data$final_cluster_assignment)]<-"Brainstem Cells"

metadata.df<-as.data.frame(combined@meta.data)

metadata.df<-metadata.df %>%
  mutate(old.ident = case_when(
    endsWith(orig.ident, "sc_brainstem") ~ "Brainstem Cells",
  ))
metadata.df$old.ident[is.na(metadata.df$old.ident)]<-"Spinalcord Cells"

combined@meta.data<-metadata.df


saveRDS(combined, file = 'levine_new/excit_integrate1_raw.rds')


combined_excit<-combined



#### Run the standard workflow for visualization and clustering ####
combined_excit<- readRDS('levine_new/excit_integrate1_raw.rds')

DefaultAssay(combined_excit) <- "integrated"

combined_excit <- ScaleData(combined_excit, verbose = FALSE)
combined_excit <- RunPCA(combined_excit, verbose = FALSE)

ElbowPlot(combined_excit, ndims = 50)

combined_excit <- RunUMAP(combined_excit, reduction = "pca", dims = 1:40)
combined_excit <- FindNeighbors(combined_excit, reduction = "pca", dims = 1:40)
combined_excit <- FindClusters(combined_excit, resolution = 5)

# Visualization

p1 <- DimPlot(combined_excit, reduction = "umap", repel = TRUE,shuffle = TRUE, label = TRUE) 
p1






# combining the metadata fields of predicted.id and Final.clusters so that one can use it 

select.df<-combined_excit@meta.data


for(i in 1:length(select.df$final_cluster_assignment)) {
  if(select.df$final_cluster_assignment[i]=='Brainstem Cells') {
    select.df$final_cluster_assignment[i]<-select.df$predicted.id[i]
  }

}


# put select.df back into object
combined_excit@meta.data<-select.df


DimPlot(combined_excit, reduction = "umap", repel = TRUE,shuffle = TRUE, label = TRUE, group.by = 'final_cluster_assignment') 


saveRDS(combined_excit, file = 'levine_new/excit_integrate1_clustered_raw.rds')

combined_excit<-readRDS('levine_new/excit_integrate1_clustered_raw.rds')
combined_excit1<-readRDS('levine_new/excit_integrate1_clustered.rds')

