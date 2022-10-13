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

# inside system
setwd('D:/Presh/h5d_dataset/')
# inside drive
setwd('Y:/PD/spinal_cord/')
# on the outside drive
setwd('U:/eng_research_economo/PD/spinal_cord/levine/')

# load in functions file
# inside drive
source('Y:/PD/spinal_cord/levine_new/scripts/brainstem_celltype/functions.R')



# source the tree seurat functions
# inside drive
source('Y:/PD/spinal_cord/levine_new/scripts/brainstem_celltype/tree_functions_seurat.R')


get.markers<- function(node){
  leaves<-GetLeavesOnly(combined_excit, tree, node)
  leaves<-subset(tree.df, node %in% leaves, select = c('cluster'))$cluster
  node.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.6,
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
combined_excit<-readRDS('levine_new/excit_integrate1_clustered.rds')

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
for (i in 0:((nrow(max2)-2)/2)){
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

# combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clusters']])

######### Hierarchial clustering #########
combined_excit<-BuildClusterTree(combined_excit,assay="integrated", dims = 1:40)
PlotClusterTree(combined_excit, direction = "downwards")

tree <- Tool(combined_excit, slot = 'BuildClusterTree') # This extracts the phylo tree from our object

tree.df<-data.frame(node=c(1:length(tree[['tip.label']])), cluster=tree[["tip.label"]])


######### Group 1 -  node 45  #########

top10markers.45<-get.markers(45)

cells.int<-get.leaves.plot(45)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()


p2<-FeaturePlot(combined_excit, features = c('Reln'), cols = c('lightgreen', 'red'), label = T, keep.scale = 'all')
p1+p2




######### Group 2 -  node 47  #########
top10markers.47<-get.markers(47)

cells.int<-get.leaves.plot(47)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()


p2<-FeaturePlot(combined_excit, features = c('Cpne8','Syt10','Trhr'), cols = c('lightgreen', 'red'), label = T, keep.scale = 'all')
p1+p2




######### Group 3 -  node 53  #########
top10markers.53<-get.markers(53)

cells.int<-get.leaves.plot(53)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Nts'), cols = c('lightgreen', 'red'), label = T, keep.scale = 'all')
p1+p2



######### Group 4 -  node 57  #########
top10markers.57<-get.markers(57)

cells.int<-get.leaves.plot(57)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Trhde','Otof'), cols = c('lightgreen', 'red'), label = T)
p1+p2


######### Group 5 -  node 58  #########
top10markers.58<-get.markers(58)

cells.int<-get.leaves.plot(58)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Adarb2','Maf','Tox','Ryr3'), cols = c('lightgreen', 'red'), label = T)
p1+p2



######### Group 6 -  node 51  #########
top10markers.51<-get.markers(51)

cells.int<-get.leaves.plot(51)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Fbn2','Sox5','Cadps2'), cols = c('lightgreen', 'red'), label = T)
p1+p2


######### Group 7 -  node 55  #########
top10markers.55<-get.markers(55)

cells.int<-get.leaves.plot(55)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Zeb2','Slit2','Man1a', 'Rgs6','Gpc6'), cols = c('lightgreen', 'red'), label = F)
p1+p2



######### Group 8 -  node 56  #########
top10markers.56<-get.markers(56)

cells.int<-get.leaves.plot(56)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Htr2c','Meis2','Slit2'), cols = c('lightgreen', 'red'), label = T)
p1+p2



# try to split up node 56 
# node 63
top10markers.63<-get.markers(63)

cells.int<-get.leaves.plot(63)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'seurat_clust_type') + NoLegend()

p2<-FeaturePlot(combined_excit, features = c('Cntn5','Slit2'), cols = c('lightgreen', 'red'), label = T)
p1+p2


################################################################################
################################################################################
# Final annotations of the cluster families 

combined_excit <- RenameIdents(combined_excit, "19"="Reln","37"="Reln","8"="Reln","10"="Reln","23"="Reln","5"="Sox5_Adarb2","14"="Sox5_Adarb2","16"="Sox5_Adarb2","4"="Sox5_Adarb2","46"="Sox5","0"="Sox5","28"="Sox5","7"="Sox5","50"="Sox5","22"="Sox5","38"="Sox5","32"="Rgs6","42"="Rgs6","21"="Rgs6","31"="Rgs6","35"="Rgs6","44"="Cpne4_Rgs6","36"="Cpne4_Rgs6","43"="Cpne4_Rgs6","34"="Nts_Adarb2","18"="Nts_Adarb2","24"="Nts_Adarb2","27"="Cck_Adarb2","48"="Cck_Adarb2","52"="Cck_Adarb2","17"="Cck_Adarb2","54"="Cck_Adarb2","49"="Pnoc","51"="Meis2","25"="Meis2","39"="Meis2","9"="Meis2","41"="Meis2","45"="Meis2","11"="Meis2","13"="Meis2","3"="Adarb2","15"="Adarb2","20"="Adarb2","33"="Adarb2","12"="Adarb2","1"="Adarb2","30"="Adarb2","40"="Adarb2","2"="Adarb2","53"="Adarb2","29"="Zhfx3","47"="Zhfx3","6"="Tcf4","26"="Tcf4")
combined_excit <- StashIdent(combined_excit, save.name = 'Fam.label')

colors.needed<-CustomColor(combined_excit,seed = F)
options(ggrepel.max.overlaps = 50)
p1<-DimPlot(combined_excit, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE, cols = colors.needed) 
p1



################################################################################
# Get the cell composition in a seurat cluster
combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clusters']])
cells.int<-WhichCells(combined_excit, idents = c(49))
data.int<-combined_excit[,cells.int]
table(data.int[['Final.clusters']])






# To highlight the cells belonging to a particular ident, and seeing its distribution across seurat_clusters
combined_excit<-SetIdent(combined_excit, value = combined_excit[['Final.clusters']])
cells.int<-WhichCells(combined_excit, idents = c('Excit.16'))
combined_excit<-SetIdent(combined_excit, value = combined_excit[['Final.clusters']])
p3<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int) + NoLegend()
p3





pt<-as.data.frame(table(combined_excit[['Final.clusters']]))
pt$Var1 <- as.character(pt$Var1)
library(RColorBrewer)
ggplot(pt, aes(x = Var1, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

