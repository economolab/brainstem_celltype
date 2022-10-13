################################################################################
# FIG 2
################################################################################

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(gplots)
library(ggplot2)
library(Polychrome)
library(viridis)

# inside system
setwd('D:/Presh/h5d_dataset/')
# inside drive
setwd('Y:/PD/spinal_cord/')
# on the outside drive
setwd('U:/eng_research_economo/PD/spinal_cord')
# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')


combined<-readRDS('levine_new/levine_integrate1_clustered.rds')

################################################################################
# FIG a - integrated dataset
################################################################################



p1<-DimPlot(combined, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE, group.by = 'old.ident', cols = c('#87CEFA','#F08080')) 
p1


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


unchosen.clusters<-unique(stacked_plot.data$cluster[stacked_plot.data$percent10==FALSE | stacked_plot.data$percent == 1])
unchosen.clusters<-sort(unchosen.clusters)
unchosen.clusters
length(unchosen.clusters)



# relabel clusters as above or below 10
# creating an ident for labeling the clusters as greater than or less than 10% mix
cluster.mix = c()
seurat_clusters<-  combined[['seurat_clusters']]
for (i in 1:length(seurat_clusters[[1]])){
  if (seurat_clusters[[1]][i] %in% unchosen.clusters){
    cluster.mix <- c(cluster.mix, 'Homogeneous')
  }
  else{
    cluster.mix <- c(cluster.mix, 'Heterogeneous')
  }
}

# adding it to the combined object 
combined[['cluster_comp']]<-cluster.mix


p4<-DimPlot(combined, reduction = "umap", label = FALSE, repel = TRUE,shuffle = TRUE, group.by = 'cluster_comp', cols = c('#8250C4','#DDA0DD')) 
p4




################################################################################
# FIG b - excit and inhib dataset 
################################################################################



combined_excit<-readRDS('levine_new/excit_integrate1_clustered.rds')
combined_inhib<-readRDS('levine_new/inhib_integrate1_clustered.rds')

options(ggrepel.max.overlaps = 50)

colors.needed.excit_type<-CustomColor(combined_excit,cluster.by='final_cluster_assignment', seed = T)
p2<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment', cols = colors.needed.excit_type) + NoLegend() 
p2


colors.needed.inhib_type<-CustomColor(combined_inhib,cluster.by='final_cluster_assignment', seed = T)
p2<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment', cols = colors.needed.inhib_type) + NoLegend() 
p2



# Adding laminae family for excit and inhib dataset 


lam_fam<-read.csv('Y:/PD/spinal_cord/levine_new/lamina_info_modified.csv', header = T)


rownames(lam_fam) <- lam_fam$Cluster

### for combined_excit
laminae.meta.fam<-c()
for (i in combined_excit@meta.data$final_cluster_assignment){
  
  type<- i
  if(type %in% lam_fam$Cluster){
    lamina.info<- lam_fam[type,]$Family_lamina
    laminae.meta.fam<-c(laminae.meta.fam,lamina.info)
    
  }
  else{
    laminae.meta.fam<-c(laminae.meta.fam,'useless')
    
  }
}

combined_excit[['Laminae_fam']]<-laminae.meta.fam


### for combined_inhib
laminae.meta.fam<-c()
for (i in combined_inhib@meta.data$final_cluster_assignment){
  
  type<- i
  if(type %in% lam_fam$Cluster){
    lamina.info<- lam_fam[type,]$Family_lamina
    laminae.meta.fam<-c(laminae.meta.fam,lamina.info)
    
  }
  else{
    laminae.meta.fam<-c(laminae.meta.fam,'useless')
    
  }
}

combined_inhib[['Laminae_fam']]<-laminae.meta.fam



p2<-DimPlot(combined_excit, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Laminae_fam") 
p2


p3<-DimPlot(combined_inhib, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Laminae_fam") 
p3



################################################################################
# Mapping brainstem to spinal cord dataset 
################################################################################

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

# reference.data<-NormalizeData(reference.data)
# reference.data<-FindVariableFeatures(reference.data, selection.method = "vst", nfeatures = 2000)
# reference.data <- ScaleData(reference.data, verbose = T)
# reference.data <- RunPCA(reference.data, verbose = T)

ElbowPlot(reference.data, ndims = 50)

# find integration anchors between the 2 datasets
anchors <- FindTransferAnchors(reference = reference.data, query = bs_object,
                               dims = 1:25, reference.reduction = "pca")

saveRDS(anchors,'levine_new/anchors_bs_pca.rds')


anchors<-readRDS('levine_new/anchors_bs_pca.rds')

predictions <- TransferData(anchorset = anchors, refdata = reference.data@meta.data$final_cluster_assignment,
                            dims = 1:25)

predicted.id<-predictions$predicted.id

bs_object[['predicted.id']]<-predicted.id



# UNIMODAL UMAP PROJECTION #####

# running UMAP on levine dataset
ref <- RunUMAP(reference.data, dims = 1:25, reduction = "pca", return.model = TRUE)

anchors<-readRDS('levine_new/anchors_bs_pca.rds')

combined<- MapQuery(anchorset = anchors, reference = ref, query = bs_object,
                    refdata = ref@meta.data$Final.clusters, reference.reduction = "pca", reduction.model = "umap")

################# create coarse labels for the brainstem cells final_cluster_assignment
# create a table for the final cluster types and coarse types from the reference or spinal cord data
reference.data.meta<-ref@meta.data
coarse.table<-as.data.frame(reference.data.meta[,c('final_cluster_assignment','final_coarse_types')] %>% distinct(final_cluster_assignment, final_coarse_types, .keep_all = TRUE))


# make a metadata field in the combined brainstem object 
combined.meta<-combined@meta.data
coarse_types<-c()
for (i in 1:length(combined.meta$predicted.id)){
  value<-combined.meta$predicted.id[i]
  assign<-as.character(subset(coarse.table, final_cluster_assignment == value, select = c('final_coarse_types'))[[1]])
  coarse_types<-c(coarse_types,assign)
}

combined[['final_coarse_types']]<-coarse_types
View(combined@meta.data)


################## create simplified labels for the datasets 
# use gray for non-neuronal and another color for neuronal (light blue)
simplified<-c()
for(i in 1:length(reference.data.meta$final_coarse_types)){
  value<-reference.data.meta$final_coarse_types[i]
  if(value != 'Neurons'){
    simplified<-c(simplified, 'Non-neurons')
  }
  else{
    simplified<-c(simplified, 'Neurons')
  }
}

ref[['simplify_label']]<-simplified



simplified<-c()
for(i in 1:length(combined@meta.data$final_coarse_types)){
  value<-combined@meta.data$final_coarse_types[i]
  if(value != 'Neurons'){
    simplified<-c(simplified, 'Non-neurons')
  }
  else{
    simplified<-c(simplified, 'Neurons')
  }
}

combined[['simplify_label']]<-simplified



# options(ggrepel.max.overlaps = Inf)

p1 <- DimPlot(ref, reduction = "umap", group.by = "simplify_label", label = TRUE, label.size = 3.5,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(combined, reduction = "ref.umap", group.by = "simplify_label", label = TRUE,
              label.size = 3.5, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 

p2






