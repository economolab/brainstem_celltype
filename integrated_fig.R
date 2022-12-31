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
# Load Packages
library(tidyverse)
library(patchwork)
library(scales)


setwd('Y:/PD/spinal_cord/')

# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')


combined<-readRDS('levine_new/levine_integrate1_clustered.rds')

################################################################################
# FIG - integrated dataset
################################################################################

cell_type <- combined@meta.data$old.ident
cluster <- combined@meta.data$seurat_clusters
data <- data.frame(cell_type, cluster)
data <-data %>% group_by_all %>% count
data_bs<- data[data$cell_type== 'Brainstem Cells',]
data_sc<- data[data$cell_type== 'Spinalcord Cells',]

# merge the 2 data frames
total <- merge(data_bs,data_sc,by="cluster", all.x = TRUE)

total$n.y[is.na(total$n.y)] <- 0
total$ratio <- total$n.x/(total$n.y + total$n.x)
total.sort <- total[order(total$ratio),]

#total.sort <- factor(total.sort, levels = total.sort$cluster)
data$cluster<- ordered(data$cluster, levels = total.sort$cluster)

stacked_plot<-ggplot(data, aes(fill=cell_type, y=n, x=cluster)) +
  geom_bar(position="fill", stat="identity")+ geom_hline(yintercept=0.05)+ geom_hline(yintercept=0.95) + coord_flip()
stacked_plot


saveRDS(combined,'levine_new/levine_integrate1_clustered.rds')

################################################################################
# FIG - excit and inhib dataset 
################################################################################

combined_excit<-readRDS('levine_new/excit_integrate1_clustered_raw.rds')
combined_inhib<-readRDS('levine_new/inhib_integrate1_clustered_raw.rds')


# reorder the factor levels in the excit dataset 
combined_excit@meta.data$final_cluster_assignment<-factor(combined_excit@meta.data$final_cluster_assignment, levels = c("Excit-1","Excit-2","Excit-3","Excit-4","Excit-5","Excit-6","Excit-7","Excit-8","Excit-9","Excit-10","Excit-11","Excit-12","Excit-13","Excit-14","Excit-15","Excit-16","Excit-17","Excit-18","Excit-19","Excit-20","Excit-21","Excit-22","Excit-23","Excit-24","Excit-25","Excit-26","Excit-27","Excit-28","Excit-29","Excit-30","Excit-31","Excit-32","Excit-33","Excit-34","Excit-35","Excit-36","Excit-37","Excit-38"))

combined_inhib@meta.data$final_cluster_assignment<-factor(combined_inhib@meta.data$final_cluster_assignment, levels = c("Inhib-1","Inhib-2","Inhib-3","Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9","Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14","Inhib-15","Inhib-16","Inhib-17","Inhib-18","Inhib-19","Inhib-20","Inhib-21","Inhib-22","Inhib-23","Inhib-24","Inhib-25","Inhib-26","Inhib-27"))



options(ggrepel.max.overlaps = 100)

#colors.needed.excit_type<-CustomColor(combined_excit,cluster.by='final_cluster_assignment', seed = T)
p2<-DimPlot(combined_excit, reduction = "umap", label = FALSE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment',cols= viridis(38, option = "plasma")) + NoLegend()
p2

#colors.needed.inhib_type<-CustomColor(combined_inhib,cluster.by='final_cluster_assignment', seed = T)
p2<-DimPlot(combined_inhib, reduction = "umap", label = FALSE, repel = TRUE,shuffle = TRUE, group.by = 'final_cluster_assignment', cols = viridis(27, option = "plasma")) + NoLegend()
p2



# creating the legend block for the plots
colfunc <- rev(viridis(27, option = 'plasma'))
legend_image <- as.raster(matrix(colfunc, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 0, 0, 1,1)


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

combined_excit@meta.data$Laminae_fam<-factor(combined_excit@meta.data$Laminae_fam, levels = c('1/2o','1/2o/2i','2/3','3','3/4','4','4/5','5','6','7','10'))

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



saveRDS(combined_excit, 'levine_new/excit_integrate1_clustered_raw.rds')
saveRDS(combined_inhib, 'levine_new/inhib_integrate1_clustered_raw.rds')

