################################################################################
# FIG 2 - heatmap and cluster tree
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
library(tidyverse)

# inside system
setwd('D:/Presh/h5d_dataset/')
# inside drive
setwd('Y:/PD/spinal_cord/')
# on the outside drive
setwd('U:/eng_research_economo/PD/spinal_cord')
# Load in custom functions:
source('levine/scripts/functions.R')


combined_excit<-readRDS('levine/method1/excit_integrate1_clustered.rds')
combined_inhib<-readRDS('levine/method1/inhib_integrate1_clustered.rds')
combined<- readRDS('levine/integrated_levine_label_thresh10.rds')

##### heatmaps and cluster trees #####
#### Heatmap for the different sets of gene markers that distinguish the different cell types


excit.markers<-c('Slc17a6', 'Slc17a7')
inhib.markers<-c('Gad1', 'Gad2','Slc32a1', 'Slc6a1','Slc6a5')
glial_markers<-c('Aqp4','Mbp','Trem2')
# cholin_markers <- c('Slc5a7', 'Chat')
int_genes<-c(excit.markers, inhib.markers, glial_markers)


file<-readRDS('Y:/PD/final_meta_dataset.rds')





# excit.data<-NormalizeMarkers(combined_excit, int_genes, method = 'max',cluster_by = 'seurat_clusters')
# rownames(excit.data)<-paste("E",excit.data$cluster,sep='')
# rownames(excit.data)<-factor(rownames(excit.data), levels = unique(as.character(rownames(excit.data))))
# 
# 
# inhib.data<-NormalizeMarkers(combined_inhib, int_genes, method = 'max',cluster_by = 'seurat_clusters')
# rownames(inhib.data)<-paste("I",inhib.data$cluster,sep='')
# rownames(inhib.data)<-factor(rownames(inhib.data), levels = unique(as.character(rownames(inhib.data))))
# 
# 
# 
# 
# excit.data.2 <- excit.data[1:ncol(excit.data)-1] %>%
#   rownames_to_column() %>%
#   gather(colname, value, -rowname)
# 
# 
# excit.data.2$rowname<-factor(excit.data.2$rowname, levels = unique(as.character(excit.data.2$rowname)))
# excit.data.2$colname<-factor(excit.data.2$colname, levels = unique(as.character(excit.data.2$colname)))
# 
# 
# inhib.data.2 <- inhib.data[1:ncol(inhib.data)-1] %>%
#   rownames_to_column() %>%
#   gather(colname, value, -rowname)
# 
# inhib.data.2$rowname<-factor(inhib.data.2$rowname, levels = unique(as.character(inhib.data.2$rowname)))
# inhib.data.2$colname<-factor(inhib.data.2$colname, levels = unique(as.character(inhib.data.2$colname)))
# 
# 
# 
# ggplot(excit.data.2, aes(x = rowname, y = colname, fill = value)) + scale_fill_gradient(low = 'blue', high = 'red') + geom_tile()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggplot(inhib.data.2, aes(x = rowname, y = colname, fill = value)) + 
#   geom_tile()  + scale_fill_gradient(low = 'blue', high = 'red') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# 
# 





####### trees clustering ######

# Reading in the excit clustered dataset 
combined_excit<-readRDS('levine/method1/excit_integrate1_clustered.rds')

combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clusters']])
######### Hierarchial clustering #########
combined_excit<-BuildClusterTree(combined_excit,assay="integrated", dims = 1:40)
PlotClusterTree(combined_excit, direction = "downwards")






#### EXCIT FAM #####
combined_excit<-readRDS('levine/method1/excit_integrate1_clustered.rds')
combined_excit <- RenameIdents(combined_excit, "19"="Reln","37"="Reln","8"="Reln","10"="Reln","23"="Reln","5"="Sox5_Adarb2","14"="Sox5_Adarb2","16"="Sox5_Adarb2","4"="Sox5_Adarb2","46"="Sox5","0"="Sox5","28"="Sox5","7"="Sox5","50"="Sox5","22"="Sox5","38"="Sox5","32"="Rgs6","42"="Rgs6","21"="Rgs6","31"="Rgs6","35"="Rgs6","44"="Cpne4_Rgs6","36"="Cpne4_Rgs6","43"="Cpne4_Rgs6","34"="Nts_Adarb2","18"="Nts_Adarb2","24"="Nts_Adarb2","27"="Cck_Adarb2","48"="Cck_Adarb2","52"="Cck_Adarb2","17"="Cck_Adarb2","54"="Cck_Adarb2","49"="Pnoc","51"="Meis2","25"="Meis2","39"="Meis2","9"="Meis2","41"="Meis2","45"="Meis2","11"="Meis2","13"="Meis2","3"="Adarb2","15"="Adarb2","20"="Adarb2","33"="Adarb2","12"="Adarb2","1"="Adarb2","30"="Adarb2","40"="Adarb2","2"="Adarb2","53"="Adarb2","29"="Zhfx3","47"="Zhfx3","6"="Tcf4","26"="Tcf4")
combined_excit <- StashIdent(combined_excit, save.name = 'Fam.label')

colors.needed<-CustomColor(combined_excit,seed = F)
options(ggrepel.max.overlaps = 50)
p1<-DimPlot(combined_excit, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE, cols = colors.needed) 
p1



p2<-FeaturePlot(combined_excit, features = c('Zfhx3'), cols = c('lightgreen', 'red'), label = F, min.cutoff = 0, max.cutoff = 5)
p2




#### INHIB FAM #####

combined_inhib<-readRDS('levine/method1/inhib_integrate1_clustered.rds')
combined_inhib <- RenameIdents(combined_inhib,"53"="Adarb2","51"="Adarb2_Chat","25"="Adarb2","47"="Rorb","18"="Rorb","48"="Rorb","1"="Chrm3","41"="Chrm3","22"="Rorb","0"="Rorb","7"="Rorb","23"="Adarb2","44"="Adarb2","46"="Adarb2","26"="Adarb2","36"="Adarb2","29"="Adarb2","16"="Adarb2","28"="Adarb2","39"="Adarb2","21"="Adarb2","5"="Slit2","43"="Adarb2","33"="Slit2","40"="Slit2","11"="Slit2","32"="Slit2","10"="Rorb","31"="Slit2","45"="Slit2","3"="Chrm3","14"="Chrm3","2"="Chrm3","8"="Chrm3","9"="Chrm3","50"="Adarb2","6"="Adarb2","35"="Adarb2","38"="Slit2","42"="Slit2","27"="Gal","12"="Npy","52"="Adarb2","15"="Npy","34"="Npy","37"="Npy","17"="Rorb_Chrm3","30"="Rorb_Chrm3","13"="Rorb_Chrm3","19"="Rorb_Chrm3","20"="Adarb2_Gal","24"="Rorb_Adarb2_Gal","4"="Rorb_Adarb2_Gal","49"="Adarb2")
combined_inhib <- StashIdent(combined_inhib, save.name = 'Fam.label')

colors.needed<-CustomColor(combined_inhib,seed = T)
p2<-DimPlot(combined_inhib, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE, cols = colors.needed) 
p2



p2<-FeaturePlot(combined_inhib, features = c('Slit2'), cols = c('lightgreen', 'red'), label = F, min.cutoff = 0, max.cutoff = 5)
p2
















