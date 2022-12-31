################################################################################
# FIGURES FOR SPINAL CORD
################################################################################

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



################################################################################
# NEUROTRANSMITTER STATUS OF SPINAL CORD
################################################################################


sc_neuron<-readRDS('levine_new/levine_dataset/integrate2.rds')

sc_neuron <- FindClusters(sc_neuron, resolution = 12) # 1.5 # Resolution above 1 lets you break the data up into more communities


DimPlot(sc_neuron, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 1.2, label = TRUE)


### Method 1 of labeling ####
# splitting the clusters into excit, inhib, acetylcholine and non-neuronal 

E_markers <- c('Slc17a7', 'Slc17a6')
I_markers_gaba <- c('Gad1', 'Gad2','Slc32a1', 'Slc6a1')
I_markers_glyc<-c('Slc6a5')
astro_markers <- c('Gfap', 'Aqp4')
sero_markers <- c('Slc6a4')
dopa_markers <- c('Th', 'Slc6a2')
microglia_markers <- c('Trem2')
endo_markers <- c('Flt1', 'Ly6c1')
oligo_markers <- c('Mbp', 'Mobp')
cholin_markers <- c('Slc5a7', 'Chat')


int_genes<-c(E_markers, I_markers_gaba, I_markers_glyc, astro_markers, sero_markers, dopa_markers, microglia_markers, endo_markers, oligo_markers, cholin_markers)


object_means<- NormalizeMarkers(sc_neuron, int_genes, method = 'min_max', cluster_by = 'seurat_clusters', means = TRUE)

object_means$E<-pmax(object_means$Slc17a6,object_means$Slc17a7)
object_means$I<-pmax(object_means$Gad1,object_means$Gad2,object_means$Slc6a1,object_means$Slc32a1, object_means$Slc6a5)
object_means$astro<-pmax(object_means$Gfap , object_means$Aqp4)
object_means$microglia <- object_means$Trem2
object_means$endo<- pmax(object_means$Flt1 , object_means$Ly6c1)
object_means$oligo<-pmax(object_means$Mbp , object_means$Mobp)
object_means$cholin<-pmax(object_means$Chat , object_means$Slc5a7)


interested <- object_means[c('E','I','cholin')]

interested$max_type<-colnames(interested)[apply(interested,1,which.max)]
interested$max_value<- apply(interested[,1:(ncol(interested)-1)], 1, max)

interested$second_type <- apply(interested[,1:(ncol(interested)-2)], 1, FUN = function(x) sort(x, decreasing = TRUE)[2])

interested$cluster<- object_means$cluster


interested %>%
  group_by(max_type) %>%
  summarize("count" = n())

# write out the interested df to a csv file, and then create the idents label using Excel 
write.csv(interested,"Y:/PD/spinal_cord/levine_new/figures/fig1/sc_labels_fig1a.csv", row.names = FALSE)

# In Excel, put quotes using CHAR(34)&<column>CHAR(34) and merge 2 columns using cellA&"="&cellB
# or use this directly -> =(CHAR(34)&G2&CHAR(34))&"="&(CHAR(34)&D2&CHAR(34))


# Need to rename seurat_clusters by their types 
sc_neuron <- SetIdent(sc_neuron, value = sc_neuron[['seurat_clusters']])
sc_neuron <- RenameIdents(sc_neuron,"0"="I","1"="I","10"="I","11"="E","12"="I","13"="E","14"="I","15"="I","16"="E","17"="I","18"="I","19"="I","2"="I","20"="E","21"="I","22"="E","23"="I","24"="E","25"="E","26"="E","27"="I","28"="I","29"="I","3"="E","30"="E","31"="I","32"="E","33"="E","34"="I","35"="I","36"="I","37"="E","38"="E","39"="E","4"="I","40"="E","41"="I","42"="I","43"="E","44"="I","45"="I","46"="I","47"="I","48"="cholin","49"="I","5"="I","50"="E","51"="E","52"="E","53"="I","54"="E","55"="E","56"="E","57"="I","58"="E","59"="E","6"="E","60"="E","61"="E","62"="E","63"="I","64"="I","65"="E","66"="cholin","67"="E","68"="cholin","69"="E","7"="E","70"="E","71"="I","72"="I","73"="I","74"="E","75"="E","76"="I","77"="E","78"="E","79"="cholin","8"="E","80"="I","81"="I","82"="I","83"="E","84"="E","85"="E","86"="I","87"="cholin","88"="I","89"="I","9"="E","90"="E","91"="E","92"="cholin","93"="E","94"="I","95"="I","96"="E")
sc_neuron <- StashIdent(sc_neuron, save.name = 'type.label')


# options(ggrepel.max.overlaps = 50)


p1<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "type.label", cols = c('#FFD700','#FF4500','#9370DB')) 
p1



#### Method 2 of labeling #### - This is what I've used 


sc_neuron.meta<-sc_neuron@meta.data


# adding your own labels to these cells 
my.label<-c()
for (i in 1:nrow(sc_neuron.meta)){
  if(startsWith(as.character(sc_neuron.meta[i,'final_cluster_assignment']), 'Excit')){
    my.label<-c(my.label,'Exc')
  }
  else if(startsWith(as.character(sc_neuron.meta[i,'final_cluster_assignment']), 'Inhib')){
    my.label<-c(my.label,'Inh')
  }
  else{
    my.label<-c(my.label, 'Ach')
  }
}
sc_neuron.meta$my_labels<-my.label 
sc_neuron[['my.label']]<-my.label


p1<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "my.label", cols = c('#9370DB', '#FF4500','#FFD700'))
p1


# 9370DB - purple
# FFD700 - yellow 
# FF4500 - red



################################################################################
# FIG creating histograms for spinal cord data based on nFeature and nCounts
################################################################################


### UMIS
umi = bs_object@meta.data$nCount_RNA
par(mar = c(4,2,4,2))
hist(umi , breaks=200 , col='#B0E0E6' , border=F , main="", xlab="# UMIs per cell", xlim=c(0,25000))

### Genes
genes = bs_object@meta.data$nFeature_RNA
par(mar = c(4,2,4,2))
hist(genes , breaks=200 , col='#3CB371' , border=F , main="", xlab="# Unique genes per cell", ylim = c(0,1000))


################################################################################
# FIG Grouping the laminae and lineage into families 
################################################################################


lam_fam <- read.csv('Y:/PD/spinal_cord/levine_new/lamina_info_modified.csv', header = T)
rownames(lam_fam) <-  lam_fam$Cluster

laminae.meta.fam<-c()
lineage.meta.fam<-c()
for (i in sc_neuron@meta.data$final_cluster_assignment){
  
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

sc_neuron[['Laminae_fam']]<-laminae.meta.fam
sc_neuron[['Lineage_fam']]<-lineage.meta.fam



sc_neuron@meta.data$Laminae_fam<-factor(sc_neuron@meta.data$Laminae_fam, levels = c('1/2o','1/2o/2i','2/3','3','3/4','4','4/5','5','6','7','8','9','10'))

p2<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Laminae_fam") 
p2

p2<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Lineage_fam") 
p2


saveRDS(sc_neuron, 'levine_new/levine_dataset/integrate2.rds')

