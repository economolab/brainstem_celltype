#########################################################################
## MARKER ANALYSIS FOR THE EXCIT DATASET 
#########################################################################


library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(future)
library(gplots)
library(ggplot2)
library(Polychrome)
library(formattable)


setwd('Y:/PD/spinal_cord/')


# load in functions file
source('Y:/PD/spinal_cord/levine_new/scripts/brainstem_celltype/functions.R')

# source the tree seurat functions
source('Y:/PD/spinal_cord/levine_new/scripts/brainstem_celltype/tree_functions_seurat.R')


# function to get markers from a node
get.markers<- function(node){
  leaves<-GetLeavesOnly(combined_excit, tree, node)
  leaves<-subset(tree.df, node %in% leaves, select = c('cluster'))$cluster
  node.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.4,
                            test.use = 'wilcox',
                            min.diff.pct = 0.3, 
                            only.pos = TRUE)
  
  
  node.markers.top10table <- node.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
  node.markers.top10<-rownames(node.markers.top10table)
  return(node.markers.top10table)
  
}


get.leaves.plot<-function(node){
  leaves<-GetLeavesOnly(combined_excit, tree, node)
  leaves<-subset(tree.df, node %in% leaves, select = c('cluster'))$cluster
  # Highlight the cells of the cluster you're looking at 
  cells.int<-WhichCells(combined_excit, idents = leaves)
  return(cells.int)
  
}


#########################################################################

# Reading in the excit clustered dataset 
combined_excit<-readRDS('levine_new/excit_integrate1_clustered_raw.rds')

combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clusters']])

# get the majority Excit type present in each seurat cluster 
df.meta<-subset(combined_excit@meta.data,select = c('final_cluster_assignment','seurat_clusters'))
df.pct<-df.meta %>% group_by(seurat_clusters, final_cluster_assignment) %>%
  summarise(cnt = n()) %>%
  mutate(freq = formattable::percent(cnt / sum(cnt)))



max2 <- df.pct %>%                                     
  arrange(desc(freq)) %>% 
  group_by(seurat_clusters) %>%
  slice(1:2) %>% data.frame()


new_label<-c()
for (i in 0:((nrow(max2)-1)/2)){
  subset.df<- subset(max2, seurat_clusters==i)
  val<- as.numeric(subset.df$freq) > 0.8
  if (unique(val) == FALSE){
    types<-c(sub(".*-", "", subset.df$final_cluster_assignment[1]),  sub(".*-", "", subset.df$final_cluster_assignment[2]))
    types<-types[order(as.numeric(types))]
    types<-paste('Excit-',types, sep = '')
    label<-paste(types[1], types[2], sep = '/')
  }
  else{
    label<- subset.df$final_cluster_assignment[1]
  }
  new_label<-c(new_label, label)
}

label.table<-df.pct %>% group_by(seurat_clusters) %>% slice(which.max(freq)) %>% data.frame()
label.table$final.labels<-new_label

# renaming the seurat_clusters by the predominant type found 
# create an empty vector to hold the values
seurat_clust_type<-c()
for (i in combined_excit@meta.data$seurat_clusters){
  value<- as.character(subset(label.table, seurat_clusters == i, select = c('final.labels')))
  seurat_clust_type<-c(seurat_clust_type,value)
}

# input the metdata field into the object
combined_excit[['seurat_clust_type']]<- seurat_clust_type

#set it as an active ident in the object
combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clust_type']])


######### Hierarchical clustering #########
combined_excit<-BuildClusterTree(combined_excit,assay="integrated", dims = 1:40)
PlotClusterTree(combined_excit, direction = "downwards")

# saving the file
saveRDS(combined_excit,'levine_new/excit_integrate1_clustered.rds')
# reading in the file if it is already created 
combined_excit<- readRDS('levine_new/excit_integrate1_clustered.rds')

tree <- Tool(combined_excit, slot = 'BuildClusterTree') # This extracts the phylo tree from our object

tree.df<-data.frame(node=c(1:length(tree[['tip.label']])), cluster=tree[["tip.label"]])




# Giving an example of how to utilize the functions 
# node 64
top10markers.64<-get.markers(64)

cells.int<-get.leaves.plot(64)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Pde11a'), cols = c('lightgreen', 'red'), label = T)
p1+p2





















