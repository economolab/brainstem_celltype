################################################################################
# FIG 1 
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

################################################################################
# FIG b NEUROTRANSMITTER STATUS OF BRAINSTEM
################################################################################

bs_object <- readRDS('raw_bs.rds')

bs_object[['percent.mt']] <- PercentageFeatureSet(bs_object, pattern = '^mt-') # You don't need to do this for nuclei
# bs_object[['identity']] <- bs_ids[,'Cluster']

VlnPlot(bs_object, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 0.0)
plot1 <- FeatureScatter(bs_object, 
                        feature1 = 'nCount_RNA', 
                        feature2 = 'nFeature_RNA')
plot1


# Subset to include only genes that high detected levels, and cells that had low levels of mito genes.

bs_object <- subset(bs_object,
                    subset = nCount_RNA > 200
                    #& nFeature_RNA < 5000
                    & percent.mt < 5)

# Normalization:

bs_object <- NormalizeData(bs_object,
                           normalization.method = 'LogNormalize',
                           scale.factor = 10000)

# Find highest variances

bs_object <- FindVariableFeatures(bs_object,
                                  selection.method = 'vst',
                                  nfeatures = 5000) # Default is usually 2000

top20 <- head(VariableFeatures(bs_object), 20)
plot1 <- VariableFeaturePlot(bs_object)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2
# Scaling data:

all.genes <- rownames(bs_object)

bs_object <- ScaleData(bs_object, features = all.genes)

# Linear Dimensionality reduction!

bs_object <- RunPCA(bs_object, features = VariableFeatures(object = bs_object), npcs = 100)


ElbowPlot(bs_object, ndims = 100) #Shows std devs for each PC

# Choose the number of dims that will actually show separation of cell types in the clustering analyses

bs_object <- FindNeighbors(bs_object, dims = 1:100) # Dimensions here specify the number of PCs to consider in the clustering
bs_object <- FindClusters(bs_object, resolution = 5) # 1.5 # Resolution above 1 lets you break the data up into more communities

# Generating umap/tsne

bs_object <- RunUMAP(bs_object, dims = 1:100)#dims = 1:20) # Maybe make this number match the total number of clusters that come out of clustering


DimPlot(bs_object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 1.2, label = TRUE) + NoLegend()
FeaturePlot(bs_object, features = c('Tac2'))

# splitting the clusters into excit, inhib, acetylcholine and non-neuronal #####

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


object_means<- NormalizeMarkers(bs_object, int_genes, method = 'min_max', cluster_by = 'seurat_clusters', means = TRUE)

object_means$E<-pmax(object_means$Slc17a6,object_means$Slc17a7)
object_means$I<-pmax(object_means$Gad1,object_means$Gad2,object_means$Slc6a1,object_means$Slc32a1, object_means$Slc6a5)
object_means$astro<-pmax(object_means$Gfap , object_means$Aqp4)
object_means$microglia <- object_means$Trem2
object_means$endo<- pmax(object_means$Flt1 , object_means$Ly6c1)
object_means$oligo<-pmax(object_means$Mbp , object_means$Mobp)
object_means$cholin<-pmax(object_means$Chat , object_means$Slc5a7)


interested <- object_means[c('E','I','astro','microglia','endo','oligo','cholin')]

interested$max_type<-colnames(interested)[apply(interested,1,which.max)]
interested$max_value<- apply(interested[,1:(ncol(interested)-1)], 1, max)

interested$second_type <- apply(interested[,1:(ncol(interested)-2)], 1, FUN = function(x) sort(x, decreasing = TRUE)[2])

interested$cluster<- object_means$cluster


interested %>%
  group_by(max_type) %>%
  summarize("count" = n())

# write out the interested df to a csv file, and then create the idents label using Excel 
write.csv(interested,"Y:/PD/spinal_cord/levine_new/figures/fig1/bs_labels_fig1a.csv", row.names = FALSE)

# In Excel, put quotes using CHAR(34) and merge 2 columns using cellA&"="&cellB


# Need to rename seurat_clusters by their types 
bs_object <- SetIdent(bs_object, value = bs_object[['seurat_clusters']])
bs_object <- RenameIdents(bs_object,"0"="Exc","1"="Exc","10"="Exc","100"="Inh","101"="Exc","102"="Inh","103"="Exc","104"="Exc","105"="Exc","106"="Exc","107"="Inh","108"="Exc","109"="Inh","11"="Non-neuronal","110"="Exc","111"="Inh","112"="Inh","113"="Exc","114"="Exc","115"="Non-neuronal","116"="Inh","117"="Non-neuronal","118"="Inh","12"="Inh","13"="Exc","14"="Inh","15"="Inh","16"="Inh","17"="Exc","18"="Inh","19"="Inh","2"="Inh","20"="Inh","21"="Exc","22"="Inh","23"="Inh","24"="Exc","25"="Inh","26"="Exc","27"="Inh","28"="Exc","29"="Inh","3"="Exc","30"="Exc","31"="Exc","32"="Exc","33"="Inh","34"="Exc","35"="Exc","36"="Exc","37"="Exc","38"="Exc","39"="Exc","4"="Exc","40"="Exc","41"="Exc","42"="Inh","43"="Inh","44"="Inh","45"="Exc","46"="Exc","47"="Exc","48"="Exc","49"="Inh","5"="Inh","50"="Non-neuronal","51"="Inh","52"="Inh","53"="Inh","54"="Inh","55"="Exc","56"="Exc","57"="Inh","58"="Inh","59"="Inh","6"="Inh","60"="Exc","61"="Exc","62"="Inh","63"="Exc","64"="Exc","65"="Inh","66"="Exc","67"="Inh","68"="Ach","69"="Exc","7"="Inh","70"="Inh","71"="Inh","72"="Non-neuronal","73"="Inh","74"="Exc","75"="Inh","76"="Inh","77"="Exc","78"="Exc","79"="Exc","8"="Exc","80"="Exc","81"="Inh","82"="Inh","83"="Inh","84"="Inh","85"="Exc","86"="Exc","87"="Exc","88"="Exc","89"="Exc","9"="Inh","90"="Exc","91"="Inh","92"="Exc","93"="Exc","94"="Inh","95"="Exc","96"="Exc","97"="Exc","98"="Exc","99"="Exc")
bs_object <- StashIdent(bs_object, save.name = 'type.label')

colors.needed<-CustomColor(bs_object,seed = F, cluster.by = "type.label")
# options(ggrepel.max.overlaps = 50)



p1<-DimPlot(bs_object, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "type.label", cols = c('#FF4500','#FFD700','#9370DB')) 
p1

################################################################################
# FIG a creating histograms for brainstem data based on nFeature and nCounts
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
# FIG b NEUROTRANSMITTER STATUS OF SPINAL CORD
################################################################################


sc_object<-readRDS('levine_new/final_meta_dataset.rds')



# processing the labels
# sc_object@meta.data$Final.clusters<-gsub('-0','.',sc_object@meta.data$Final.clusters)
# sc_object@meta.data$Final.clusters<-gsub('-','.',sc_object@meta.data$Final.clusters)



# only interested in neurons and motor neurons, so subset that data from the sc_object, look at the final_coarse_clusters

sc_neuron<-subset(x = sc_object, subset = final_coarse_types %in% c("Neurons"))


sc_neuron<-NormalizeData(sc_neuron)
sc_neuron<-FindVariableFeatures(sc_neuron, selection.method = "vst", nfeatures = 4000)
sc_neuron <- ScaleData(sc_neuron,  verbose = T)
sc_neuron <- RunPCA(sc_neuron,verbose = T)

sc_neuron <- FindNeighbors(sc_neuron, dims = 1:40) # Dimensions here specify the number of PCs to consider in the clustering
sc_neuron <- FindClusters(sc_neuron, resolution = 3) # 1.5 # Resolution above 1 lets you break the data up into more communities

# Generating umap/tsne

sc_neuron <- RunUMAP(sc_neuron, dims = 1:2)#dims = 1:20) # Maybe make this number match the total number of clusters that come out of clustering
# bs_object <- RunTSNE(bs_object, dims = 1:100, verbose = TRUE) #dims = 1:20)

DimPlot(sc_neuron, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 1.2, label = TRUE) 


saveRDS(sc_neuron, 'sc_neuron.rds')
sc_neuron<-readRDS('sc_neuron.rds')


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

# In Excel, put quotes using CHAR(34) and merge 2 columns using cellA&"="&cellB


# Need to rename seurat_clusters by their types 
sc_neuron <- SetIdent(sc_neuron, value = sc_neuron[['seurat_clusters']])
sc_neuron <- RenameIdents(sc_neuron,"0"="Inh","1"="Inh","2"="Inh","3"="Inh","4"="Exc","5"="Inh","6"="Inh","7"="Inh","8"="Inh","9"="Exc","10"="Inh","11"="Inh","12"="Inh","13"="Exc","14"="Ach","15"="Inh","16"="Ach","17"="Inh","18"="Inh","19"="Exc","20"="Exc","21"="Inh","22"="Inh","23"="Exc","24"="Inh","25"="Exc","26"="Inh","27"="Exc","28"="Inh","29"="Inh","30"="Exc","31"="Inh","32"="Exc","33"="Inh","34"="Ach","35"="Exc","36"="Ach","37"="Inh","38"="Inh")
sc_neuron <- StashIdent(sc_neuron, save.name = 'type.label')

colors.needed<-CustomColor(sc_neuron,seed = F, cluster.by = "type.label")
# options(ggrepel.max.overlaps = 50)



p1<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "type.label", cols = c('#FFD700','#FF4500','#9370DB')) 
p1







####### need to label the seurat_clusters based on the composition of the cell types. - like Exc, Inh, Ach #####
# let us label the different types as Exc or Inh

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

sc_neuron@meta.data <-sc_neuron.meta
# # getting the percentage of each cell type in each seurat cluster
# 
# df.percent<-sc_neuron.meta %>% group_by(seurat_clusters) %>% count(sc_neuron.meta$my_labels)
# df.percent<-df.percent %>% group_by(seurat_clusters)%>% mutate(percent = n/sum(n))
# 
# # getting the max cell type for each seurat_cluster
# df.max<-df.percent%>% group_by(seurat_clusters) %>% slice(which.max(percent))
# 
# # getting the number of clusters for each labeled type
# table(df.max$`sc_neuron.meta$my_labels`)
# 
# # getting the number of cells labeled as a type and adding them
# df.max.count<-df.max%>% group_by(`sc_neuron.meta$my_labels`) %>% summarise(sum_type = sum(n))


# write out the df.max to a csv file, and then create the idents label using Excel
write.csv(df.max,"Y:/PD/spinal_cord/levine_new/figures/fig1/sc_labels_maxpercent.csv", row.names = FALSE)


# Need to rename seurat_clusters by their types
sc_neuron <- SetIdent(sc_neuron, value = sc_neuron[['seurat_clusters']])
sc_neuron <- RenameIdents(sc_neuron,"0"="Inh","1"="Exc","2"="Exc","3"="Exc","4"="Exc","5"="Inh","6"="Exc","7"="Exc","8"="Inh","9"="Exc","10"="Inh","11"="Inh","12"="Exc","13"="Exc","14"="Exc","15"="Inh","16"="Ach","17"="Exc","18"="Inh","19"="Exc","20"="Exc","21"="Exc","22"="Exc","23"="Exc","24"="Exc","25"="Exc","26"="Inh","27"="Exc","28"="Inh","29"="Exc","30"="Exc","31"="Ach","32"="Exc","33"="Inh","34"="Exc","35"="Exc","36"="Ach","37"="Exc","38"="Ach")
sc_neuron <- StashIdent(sc_neuron, save.name = 'type.label')


p1<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "my_labels", cols = c('#9370DB','#FF4500','#FFD700'))
p1

# 9370DB - purple
# FFD700 - yellow 
# FF4500 - red

################################################################################
# FIG c Grouping the laminae and lineage into families 
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

# sc_neuron_lam<-subset(sc_neuron, subset = Laminae_fam == 'useless', invert = TRUE)
# 
# sc_neuron_lin<-subset(sc_neuron, subset = Lineage_fam == 'useless', invert = TRUE)



p2<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Laminae_fam") 
p2

p2<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Lineage_fam") 
p2


################################################################################
# FIG a creating histograms for spinal cord data based on nFeature and nCounts
################################################################################

### UMIS
umi = sc_neuron@meta.data$nCount_RNA
par(mar = c(4,2,4,2))
hist(umi , breaks=2000 , col='#B0E0E6' , border=F , main="", xlab="# UMIs per cell", xlim=c(0,25000))

### Genes
genes = sc_neuron@meta.data$nFeature_RNA
par(mar = c(4,2,4,2))
hist(genes , breaks=100 , col='#B0E0E6' , border=F , main="", xlab="# Unique genes per cell", xlim = c(0,10000))








################################################################################
int.assay.genes<- sc_object@assays$integrated@var.features


sc_neuron@assays$raw@var.features<-int.assay.genes













