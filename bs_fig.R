################################################################################
# FIGURES FOR BRAINSTEM 
################################################################################

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(gplots)
library(ggplot2)
library(Polychrome)
library(scales)

setwd('Y:/PD/spinal_cord/')

# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')

################################################################################
# NEUROTRANSMITTER STATUS OF BRAINSTEM
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

bs_object <- FindNeighbors(bs_object, dims = 1:75) # Dimensions here specify the number of PCs to consider in the clustering
bs_object <- FindClusters(bs_object, resolution = 4) # 1.5 # Resolution above 1 lets you break the data up into more communities

# Generating umap/tsne

bs_object <- RunUMAP(bs_object, dims = 1:75)#dims = 1:20) # Maybe make this number match the total number of clusters that come out of clustering


DimPlot(bs_object, reduction = 'umap', group.by = 'seurat_clusters', pt.size = 1.2, label = TRUE)


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

# In Excel, put quotes using CHAR(34)&<column>CHAR(34) and merge 2 columns using cellA&"="&cellB
# or use this directly -> =(CHAR(34)&K2&CHAR(34))&"="&(CHAR(34)&H2&CHAR(34))

# Need to rename seurat_clusters by their types 
bs_object <- SetIdent(bs_object, value = bs_object[['seurat_clusters']])
bs_object <- RenameIdents(bs_object,"0"="E","1"="E","10"="I","100"="E","101"="E","102"="astro","103"="I","11"="E","12"="I","13"="I","14"="E","15"="I","16"="I","17"="I","18"="I","19"="E","2"="E","20"="E","21"="I","22"="I","23"="I","24"="I","25"="E","26"="E","27"="I","28"="E","29"="I","3"="E","30"="I","31"="E","32"="E","33"="I","34"="E","35"="E","36"="I","37"="I","38"="I","39"="E","4"="I","40"="E","41"="I","42"="E","43"="E","44"="I","45"="E","46"="I","47"="I","48"="E","49"="I","5"="I","50"="cholin","51"="I","52"="E","53"="E","54"="I","55"="E","56"="I","57"="E","58"="I","59"="E","6"="E","60"="I","61"="E","62"="I","63"="E","64"="I","65"="E","66"="I","67"="E","68"="cholin","69"="E","7"="I","70"="I","71"="I","72"="I","73"="I","74"="E","75"="E","76"="E","77"="E","78"="E","79"="E","8"="I","80"="E","81"="E","82"="E","83"="I","84"="I","85"="E","86"="E","87"="E","88"="E","89"="E","9"="E","90"="I","91"="E","92"="E","93"="E","94"="E","95"="E","96"="I","97"="E","98"="microglia","99"="I")
bs_object <- StashIdent(bs_object, save.name = 'type.label')


# options(ggrepel.max.overlaps = 50)


p1<-DimPlot(bs_object, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "type.label", cols = c('#FF4500','#FFD700','#9370DB','green','lightblue')) 
p1

# save this object 
saveRDS(bs_object, "levine_new/bs_clustered.rds")




################################################################################
# FIG OF UMAP OF BRAINSTEM DATA WITH ALL LABELS PROJECTED ON IT 
################################################################################

# labels are stored in this object 
bs_labeled<- readRDS('levine_new/bs_levine_label_raw.rds')
label.df<-bs_labeled@meta.data[,c('predicted.id','broad_type','predicted.score')]

# clusterings are stored in this object 
bs_object<-readRDS("levine_new/bs_clustered.rds")
cluster.df<-bs_object@meta.data


data_frame_merge <- merge(cluster.df, label.df,
                          by = 'row.names', all = TRUE)


data_frame_merge2 <- data_frame_merge[,-1]
rownames(data_frame_merge2) <- data_frame_merge[,1]

# putting this metadata back into bs_object 
bs_object@meta.data<-data_frame_merge2

# need to reorder the labels
bs_object@meta.data$predicted.id<-factor(bs_object@meta.data$predicted.id, levels = c("Excit-1","Excit-2","Excit-3","Excit-4","Excit-5","Excit-6","Excit-7","Excit-8","Excit-9","Excit-10","Excit-11","Excit-12","Excit-13","Excit-14","Excit-15","Excit-16","Excit-17","Excit-18","Excit-19","Excit-20","Excit-21","Excit-22","Excit-23","Excit-24","Excit-25","Excit-26","Excit-27","Excit-28","Excit-29","Excit-30","Excit-31","Excit-32","Excit-33","Excit-34","Excit-35","Excit-36","Excit-37","Excit-38","Inhib-1","Inhib-2","Inhib-3","Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9","Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14","Inhib-15","Inhib-16","Inhib-17","Inhib-18","Inhib-19","Inhib-20","Inhib-21","Inhib-22","Inhib-23","Inhib-24","Inhib-25","Inhib-26","Inhib-27","MN-alpha","PGC","Astrocytes-1","Endothelial","Ependymal"," Microglia","OPC","Oligo Progen-2","Oligos-1"))


# choosing suitable colors for plotting 

excit.colors <- colorRampPalette(c('#ab2830', '#fff900'))(27) # should be 38 total, but not all types mapped to bs cells
show_col(excit.colors)

inhib.colors <- colorRampPalette(c('#095778', '#cfe2f3'))(21) # should be 27
show_col(inhib.colors)

discrete.colors <- c('#fb9def','#c07430','#a22ba3','#ff7fa5','#3a5450','#bdcdd5','#d70254','#fa5c44')
show_col(discrete.colors)

all.colors<- c(excit.colors, inhib.colors, discrete.colors)


DimPlot(bs_object, reduction = "umap", cols = all.colors, label = T, repel = TRUE,shuffle = TRUE,  group.by = "predicted.id")


########### adding in metadata about cells which are homo or heterogeneous ############ 


bs_meta <-bs_object@meta.data
# want to color only those cells that are "heterogeneous"
integrated<- readRDS('levine_new/levine_integrate1_clustered_thresh5.rds')


bs_thresh <- subset(integrated, subset = old.ident == 'Brainstem Cells')
cells.int <- colnames(bs_thresh)

# create new meta data fields in bs_clustered_labeled 
cluster.type<-c()
mislabel <- c()
# iterate through all cells 
for (i in 1:length(rownames(bs_meta))){
  if (rownames(bs_meta)[i] %in% cells.int){
    mislabel<- c(mislabel, as.character(bs_meta$predicted.id[i]))
    cluster.type<-c(cluster.type, "Homo")
  }
  else{
    mislabel<- c(mislabel, 'wrong')
    cluster.type<-c(cluster.type, "Hetero")
  }
}

bs_meta$corrected <- mislabel
bs_meta$cluster.type<- cluster.type
bs_object@meta.data<- bs_meta

bs_object@meta.data$corrected<- factor(bs_object@meta.data$corrected, levels = c("Excit-1","Excit-2","Excit-3","Excit-4","Excit-5","Excit-6","Excit-7","Excit-8","Excit-9","Excit-10","Excit-11","Excit-12","Excit-13","Excit-14","Excit-15","Excit-16","Excit-17","Excit-18","Excit-19","Excit-20","Excit-21","Excit-22","Excit-23","Excit-24","Excit-25","Excit-26","Excit-27","Excit-28","Excit-29","Excit-30","Excit-31","Excit-32","Excit-33","Excit-34","Excit-35","Excit-36","Excit-37","Excit-38","Inhib-1","Inhib-2","Inhib-3","Inhib-4","Inhib-5","Inhib-6","Inhib-7","Inhib-8","Inhib-9","Inhib-10","Inhib-11","Inhib-12","Inhib-13","Inhib-14","Inhib-15","Inhib-16","Inhib-17","Inhib-18","Inhib-19","Inhib-20","Inhib-21","Inhib-22","Inhib-23","Inhib-24","Inhib-25","Inhib-26","Inhib-27","MN-alpha","PGC","Astrocytes-1","Endothelial","Ependymal"," Microglia","OPC","Oligo Progen-2","Oligos-1","wrong"))


all.colors<- c(all.colors,"gray80")
p1<-DimPlot(bs_object, reduction = 'umap', label=F, repel = TRUE,shuffle = TRUE,  group.by = "corrected", cols = all.colors) + NoLegend()
p1

# creating the legend block for the plots

legend_image <- as.raster(matrix(excit.colors, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 0, 0, 1,1)

legend_image <- as.raster(matrix(inhib.colors, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 0, 0, 1,1)



############################ Adding the lamina and lineage info ############################

lam_fam <- read.csv('Y:/PD/spinal_cord/levine_new/lamina_info_modified.csv', header = T)
rownames(lam_fam) <-  lam_fam$Cluster

laminae.meta.fam<-c()
lineage.meta.fam<-c()
for (i in bs_object@meta.data$corrected){
  
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

bs_object[['Laminae_fam']]<-laminae.meta.fam
bs_object[['Lineage_fam']]<-lineage.meta.fam


bs_object@meta.data$Laminae_fam<-factor(bs_object@meta.data$Laminae_fam, levels = c('1/2o','1/2o/2i','2/3','3','3/4','4','4/5','5','6','7','8','9','10'))

p2<-DimPlot(bs_object, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Laminae_fam", cols = c(hue_pal()(12), 'gray80')) 
p2

p2<-DimPlot(sc_neuron, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE,  group.by = "Lineage_fam") 
p2



saveRDS(bs_object, "levine_new/bs_clustered_labeled.rds")
saveRDS(bs_object, "D:/Presh/temp_storage/levine_new/bs_clustered_labeled.rds")


################################################################################
# ANALYZING THE MAPPING OF LABELS IN HOMO AND HETERO CLUSTERS
################################################################################



bs_object<- readRDS("levine_new/bs_clustered_labeled.rds")



score<- bs_object[['predicted.score']]
bs_object@meta.data$score<- factor(as.integer(unlist(score*100)))

colors <- colorRampPalette(c('#fff900','#ab2830','#000000'))(length(unique(bs_object@meta.data$score)))
show_col(colors)


# get all cells that are homo 
cells.homo <- colnames(subset(bs_object, subset = cluster.type == "Homo"))

p1<-DimPlot(bs_object, reduction = "umap", label = F, repel = TRUE, shuffle = TRUE, cells = cells.homo, group.by = "score", cols = colors) + NoLegend()
p1
# get all cell that are hetero 

cells.hetero <- colnames(subset(bs_object, subset = cluster.type == "Hetero"))

p2<-DimPlot(bs_object, reduction = "umap", label = F, repel = TRUE, shuffle = TRUE, cells = cells.hetero, group.by = "score", cols = colors) + NoLegend()
p2

p1+p2


# creating a palette for the plot
legend_image <- as.raster(matrix(colors, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
rasterImage(legend_image, 0, 0, 1,1)





