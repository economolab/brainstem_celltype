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
library(pheatmap)

setwd('Y:/PD/spinal_cord/')

# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')

# take inspiration from : https://www.nature.com/articles/s41586-018-0654-5/figures/10


# load in the hetero dataset or thresh_5 dataset 
combined<- readRDS('levine_new/integrated_levine_label_thresh5.rds')


combined_bs <- subset(combined, subset = old.ident == 'Brainstem Cells')
combined_sc <- subset(combined, subset = old.ident == 'Spinalcord Cells')


# trying to get names of the seurat clusters based on the majority of cell types in a certain lamina (both inhib and excit) in each cluster
combined_bs<-SetIdent(combined_bs, value = combined_bs[['seurat_clusters']])

# get the majority of cell types present in each seurat cluster 
df.meta<-subset(combined_bs@meta.data,select = c('Laminae_fam','seurat_clusters'))
df.pct<-df.meta %>% group_by(seurat_clusters, Laminae_fam) %>%
  summarise(cnt = n()) %>%
  mutate(freq = formattable::percent(cnt / sum(cnt)))



max1 <- df.pct %>%                                     
  arrange(desc(freq)) %>% 
  group_by(seurat_clusters) %>%
  slice(1) %>% data.frame()



label.table<-df.pct %>% group_by(seurat_clusters) %>% slice(which.max(freq)) %>% data.frame()
label.table$seurat_clusters<-sub("^","RNA.",label.table$seurat_clusters)
row.names(label.table) <- label.table$seurat_clusters
label.table<- label.table[c(2)]
label.table$Laminae_fam <- factor(label.table$Laminae_fam, levels = c("1/2o","1/2o/2i","2/3","3","3/4","4/5","5","6","7","8","9"))


############## HEATMAP 
marker_genes<-c('Slc32a1','Slc17a6','Slc17a7','Rorb','Npy','Chrm2','Penk','Fxyd6','Pvalb','Tac1','Reln','Sncg','Cck','Tac2','Nts','Sst', 'Calb1')

bs.avg<-as.data.frame(AverageExpression(combined_bs, features = marker_genes, group.by = 'seurat_clusters', slot = 'data', assays = 'RNA'))

bs.colmax<-apply(bs.avg,1,max)
bs.norm<-sweep(bs.avg,1,bs.colmax,FUN="/")

# adding annotation for laminae 
my_colour = list(Laminae_fam = c("1/2o" = "red","1/2o/2i" = "red","2/3" = "red","3" = "red","3/4" = "red","4/5" = "blue","5" = "blue","6" = "green","7"="green","8"= "green","9" = "green"))

# plotting the heatmap
pheatmap(bs.norm, cluster_rows = FALSE, annotation_col = label.table, annotation_colors = my_colour)

pheatmap(bs.norm,  border_color = NA, color=colorRampPalette(c("black", "red"))(100), cluster_rows = 1)




# highlight some cells
int.cells<- WhichCells(combined, ident = c(3,50))
p1<- DimPlot(combined, reduction = "umap", cells.highlight = int.cells, label = T) + NoLegend()

# plot final cluster assignments
colors.needed<- CustomColor(combined, cluster.by = 'final_cluster_assignment', seed = T)
p2<- DimPlot(combined, group.by = 'final_cluster_assignment',label = T, cols = colors.needed, label.size = 3.5, repel = TRUE,shuffle = TRUE)


p1+p2






