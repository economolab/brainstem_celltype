#########################################################################
## MAPPING BRAINSTEM DATA TO THE LEVINE DATASET
#########################################################################

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(gplots)
library(ggplot2)
library(Polychrome)


setwd('Y:/PD/spinal_cord/')

# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')

# load in query dataset 
bs_object<-readRDS('raw_bs.rds')

bs_object[['percent.mt']] <- PercentageFeatureSet(bs_object, pattern = '^mt-')

bs_object <- subset(bs_object,
                    subset = nCount_RNA > 200
                    #& nFeature_RNA < 5000
                    & percent.mt < 5)


bs_object<-NormalizeData(bs_object)
bs_object<-FindVariableFeatures(bs_object, selection.method = "vst", nfeatures = 2000)


# load in reference dataset 
reference.data<-readRDS('levine_new/final_meta_dataset.rds')
DefaultAssay(reference.data)<-'integrated'


ElbowPlot(reference.data, ndims = 50)

# find integration anchors between the 2 datasets
# anchors <- FindTransferAnchors(reference = reference.data, query = bs_object, dims = 1:25, reference.reduction = "pca")

# One can save these anchors for future use
# saveRDS(anchors,'levine_new/anchors_bs_pca.rds')


# If already saved, one can read them in 
anchors<-readRDS('levine_new/anchors_bs_pca.rds')


predictions <- TransferData(anchorset = anchors, refdata = reference.data@meta.data$final_cluster_assignment,
                            dims = 1:25)

# subset out the final predictions and scores 
predictions.subset <- subset(predictions, select = c('predicted.id','prediction.score.max'))
write.csv(predictions.subset,"levine_new/mapping.csv", row.names = TRUE)


# input the predicted cell labels in the bs_object metadata field
bs_object[['predicted.id']]<-predictions$predicted.id
bs_object[['predicted.score']]<- predictions$prediction.score.max

bs.meta<-bs_object@meta.data
# adding your own labels to these cells 
my.label<-c()
for (i in 1:nrow(bs.meta)){
  if(startsWith(bs.meta[i,'predicted.id'], 'Excit')){
    my.label<-c(my.label,'Excit')
  }
  else if(startsWith(bs.meta[i,'predicted.id'], 'Inhib')){
    my.label<-c(my.label,'Inhib')
  }
  else if(startsWith(bs.meta[i,'predicted.id'], 'MN')){
    my.label<-c(my.label,'Motor')
  }
  else{
    my.label<-c(my.label, 'Non_neuron')
  }
}


bs_object[['broad_type']]<-my.label



# saving the brainstem labeled data 
saveRDS(bs_object, 'levine_new/bs_levine_label_raw.rds')

############################### UNIMODAL UMAP PROJECTION ###########################

# running UMAP on levine dataset
ref <- RunUMAP(reference.data, dims = 1:25, reduction = "pca", return.model = TRUE)

anchors<-readRDS('levine_new/anchors_bs_pca.rds')

combined<- MapQuery(anchorset = anchors, reference = ref, query = bs_object,
                      refdata = ref@meta.data$Final.clusters, reference.reduction = "pca", reduction.model = "umap")


options(ggrepel.max.overlaps = Inf)
p1 <- DimPlot(ref, reduction = "umap", group.by = "final_coarse_types", label = TRUE, label.size = 3.5,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(combined, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3.5, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2


bs_labeled<-readRDS('levine_new/bs_levine_label_raw.rds')
# visualizing just the projected umap of the brainstem data 
colors.needed<-CustomColor(bs_labeled,cluster.by='predicted.id', seed = T)

p3<-DimPlot(bs_labeled, reduction = "ref.umap", label = TRUE, pt.size = 0.5, group.by = 'predicted.id', cols = colors.needed, shuffle = T, repel = T)
p3


################################################################################

# Applying these cell labels to the brainstem cells in the whole integrated levine and brainstem dataset #####
bs_labeled<-readRDS('levine_new/bs_levine_label_raw.rds')
# load in the whole integrated dataset which is clustered
integrated<-readRDS('levine_new/levine_integrate1_clustered.rds')
# extracting the meta data of both datasets 
integrated.meta<-integrated@meta.data
bs.meta<-bs_labeled@meta.data
# assigning final cluster ids to the integrated dataset 

for (i in 1:nrow(bs.meta)){
  cell<-rownames(bs.meta)[i]
  integrated.meta[cell,'final_cluster_assignment']<-bs.meta[cell,'predicted.id']
  integrated.meta[cell,'prediction_score']<-bs.meta[cell,'predicted.score']
}
# putting it back into integrated dataset metadata 
integrated@meta.data<-integrated.meta


saveRDS(integrated, 'levine_new/integrated_levine_label.rds')




################################################################################
# apply these cell labels to the brainstem cells in the 5% threshold refined integrated levine and brainstem dataset #####
integrated.whole<-readRDS('levine_new/integrated_levine_label.rds')
# load in the thresholded integrated dataset which is clustered
integrated.thresh<-readRDS('levine_new/levine_integrate1_clustered_thresh5.rds')
# extract the cells from this dataset into a list 
integrated<-integrated.whole[,colnames(integrated.thresh)]

integrated.meta<-integrated@meta.data

# adding your own labels to these cells 
my.label<-c()
for (i in 1:nrow(integrated.meta)){
  if(startsWith(integrated.meta[i,'final_cluster_assignment'], 'Excit')){
    my.label<-c(my.label,'Excit')
  }
  else if(startsWith(integrated.meta[i,'final_cluster_assignment'], 'Inhib')){
    my.label<-c(my.label,'Inhib')
  }
  else if(startsWith(integrated.meta[i,'final_cluster_assignment'], 'MN')){
    my.label<-c(my.label,'Motor')
  }
  else{
    my.label<-c(my.label, 'Non_neuron')
  }
}

integrated[['broad_type']]<-my.label


# adding in lamina and lineage info 
lam_fam <- read.csv('Y:/PD/spinal_cord/levine_new/lamina_info_modified.csv', header = T)
rownames(lam_fam) <-  lam_fam$Cluster

laminae.meta.fam<-c()
lineage.meta.fam<-c()
for (i in integrated@meta.data$final_cluster_assignment){
  
  type<- i
  if(type %in% lam_fam$Cluster){
    lamina.info<- lam_fam[type,]$Family_lamina
    laminae.meta.fam<-c(laminae.meta.fam,lamina.info)
    lineage.info<- lam_fam[type,]$Family_Lineage
    lineage.meta.fam<-c(lineage.meta.fam,lineage.info)
    
  }
  else{
    laminae.meta.fam<-c(laminae.meta.fam,'useless')
    lineage.meta.fam<-c(lineage.meta.fam,'useless')
  }
}

integrated[['Laminae_fam']]<-laminae.meta.fam
integrated[['Lineage_fam']]<-lineage.meta.fam

integrated<-subset(integrated, subset = Laminae_fam == 'useless', invert = TRUE)

integrated@meta.data$Laminae_fam<-factor(integrated@meta.data$Laminae_fam, levels = c('1/2o','1/2o/2i','2/3','3','3/4','4','4/5','5','6','7','8','9','10'))





saveRDS(integrated, 'levine_new/integrated_levine_label_thresh5.rds')


################################################################################


############# SPLITTING THE DATASET INTO EXCIT AND INHIB DATASETS ##############

# Separating the cells into E and I clusters
integrated<- readRDS('levine_new/integrated_levine_label_thresh5.rds')
integrated.meta<- integrated@meta.data

excit.clust.df<- integrated.meta %>% filter(integrated.meta$broad_type=='Excit')
inhib.clust.df<- integrated.meta %>% filter(integrated.meta$broad_type=='Inhib')


E.cells.interested<-rownames(excit.clust.df)
I.cells.interested<-rownames(inhib.clust.df)


sc_object <- readRDS('levine_new/final_meta_dataset.rds')
bs_object <- readRDS('levine_new/bs_levine_label_raw.rds')

setwd('Y:/PD/spinal_cord/levine_new/')


excit_sc <- sc_object[,E.cells.interested]
saveRDS(excit_sc, file = "excit_sc.rds")

inhib_sc <- sc_object[,I.cells.interested]
saveRDS(inhib_sc, file = "inhib_sc.rds")

excit_bs <- bs_object[,E.cells.interested]
saveRDS(excit_bs, file = "excit_bs.rds")

inhib_bs <- bs_object[,I.cells.interested]
saveRDS(inhib_bs, file = "inhib_bs.rds")

################################################################################

