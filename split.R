#########################################################################
## SPLITTING INTEGRATED BRAINSTEM AND LEVINE DATASET
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
source('scripts/functions.R')


# Method 3 
integrated<-readRDS('levine/integrated_levine_label_thresh10.rds')


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


object_means<- NormalizeMarkers(integrated, int_genes, method = 'min_max', cluster_by = 'seurat_clusters', means = TRUE)


object_means$E<-pmax(object_means$Slc17a6,object_means$Slc17a7)
object_means$I<-pmax(object_means$Gad1,object_means$Gad2,object_means$Slc6a1,object_means$Slc32a1, object_means$Slc6a5)
object_means$astro<-pmax(object_means$Gfap , object_means$Aqp4)
object_means$sero<-object_means$Slc6a4
object_means$dopa<-pmax(object_means$Th , object_means$Slc6a2)
object_means$microglia <- object_means$Trem2
object_means$endo<- pmax(object_means$Flt1 , object_means$Ly6c1)
object_means$oligo<-pmax(object_means$Mbp , object_means$Mobp)
object_means$cholin<-pmax(object_means$Chat , object_means$Slc5a7)



interested <- object_means[c('E','I','astro','sero','dopa','microglia','endo','oligo','cholin')]
row.names(interested)<- object_means$cluster


interested$max_type<-colnames(interested)[apply(interested,1,which.max)]
interested$max_value<- apply(interested[,1:(ncol(interested)-1)], 1, max)

interested$second_type <- apply(interested[,1:(ncol(interested)-2)], 1, FUN = function(x) sort(x, decreasing = TRUE)[2])

interested$cluster<- object_means$cluster


interested %>%
  group_by(max_type) %>%
  summarize("count" = n())



E.clusters<-interested[interested$max_type=="E",]$cluster
I.clusters<-interested[interested$max_type=="I",]$cluster


# total.clusters<-c(E.clusters,I.clusters)

E.cells.interested<-WhichCells(integrated, idents = E.clusters)
I.cells.interested<-WhichCells(integrated, idents = I.clusters)

sc_object <- readRDS('levine/modified_levine.rds')
bs_object <- readRDS('levine/bs_levine_label_raw.rds')

setwd('Y:/PD/spinal_cord/levine/method3')

excit_sc <- sc_object[,E.cells.interested]
saveRDS(excit_sc, file = "excit_sc3.rds")

inhib_sc <- sc_object[,I.cells.interested]
saveRDS(inhib_sc, file = "inhib_sc3.rds")

excit_bs <- bs_object[,E.cells.interested]
saveRDS(excit_bs, file = "excit_bs3.rds")

inhib_bs <- bs_object[,I.cells.interested]
saveRDS(inhib_bs, file = "inhib_bs3.rds")

# the 3 files are thresholded files that have been classified using the max method


# Method 1
# Take all cells labeled Excit and all Labeled Inhib and put take their cells out 
setwd('Y:/PD/spinal_cord/')
integrated<-readRDS('levine/integrated_levine_label_thresh10.rds')



integrated.meta<-integrated@meta.data

my.label<-c()
for (i in 1:nrow(integrated.meta)){
  if(startsWith(integrated.meta[i,'Final.clusters'], 'Excit')){
    my.label<-c(my.label,'Excit')
  }
  else if(startsWith(integrated.meta[i,'Final.clusters'], 'Inhib')){
    my.label<-c(my.label,'Inhib')
  }
  else{
    my.label<-c(my.label, integrated.meta[i,'Final.clusters'])
  }
}
integrated.meta$my_labels<-my.label 

excit.clust.df<-  integrated.meta%>% filter(integrated.meta$my_labels=='Excit')
inhib.clust.df<-  integrated.meta%>% filter(integrated.meta$my_labels=='Inhib')



E.cells.interested<-rownames(excit.clust.df)
I.cells.interested<-rownames(inhib.clust.df)


sc_object <- readRDS('levine/modified_levine.rds')
bs_object <- readRDS('levine/bs_levine_label_raw.rds')

# inside drive
setwd('Y:/PD/spinal_cord/levine/method1')
# outside drive
setwd('U:/eng_research_economo/PD/spinal_cord/levine')

excit_sc <- sc_object[,E.cells.interested]
saveRDS(excit_sc, file = "excit_sc1.rds")

inhib_sc <- sc_object[,I.cells.interested]
saveRDS(inhib_sc, file = "inhib_sc1.rds")

excit_bs <- bs_object[,E.cells.interested]
saveRDS(excit_bs, file = "excit_bs1.rds")

inhib_bs <- bs_object[,I.cells.interested]
saveRDS(inhib_bs, file = "inhib_bs1.rds")



































