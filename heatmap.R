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
library(ComplexHeatmap)
library(circlize)
library(seriation)

setwd('Y:/PD/spinal_cord/')

setwd('/Volumes/My\ Passport/spinal_cord/')
# Load in custom functions:
source('levine_new/scripts/brainstem_celltype/functions.R')

# take inspiration from : https://www.nature.com/articles/s41586-018-0654-5/figures/10

marker_genes<-c('Slc17a6','Slc17a7','Slc32a1','Cck','Tac2','Nts','Sst', 'Calb1','Rorb','Npy','Chrm2','Penk','Fxyd6','Pvalb','Tac1','Reln','Sncg')

# load in the hetero dataset or thresh_5 dataset 
combined<- readRDS('levine_new/integrated_levine_label_thresh5.rds')
combined<-SetIdent(combined, value = combined[['seurat_clusters']])
combined.avg.all<-as.data.frame(AverageExpression(combined, features = marker_genes, group.by = 'seurat_clusters', slot = 'data', assays = 'RNA'))


# create the subset data 
combined_bs <- subset(combined, subset = old.ident == 'Brainstem Cells')
combined_bs<-SetIdent(combined_bs, value = combined_bs[['seurat_clusters']])
bs.avg<-as.data.frame(AverageExpression(combined_bs, features = marker_genes, group.by = 'seurat_clusters', slot = 'data', assays = 'RNA'))
# since the brainstem dataset is lacking one cluster, we need to remove it from all 
common.clusters<-intersect(colnames(combined.avg.all), colnames(bs.avg))

combined.avg<-combined.avg[common.clusters]
combined.colmax<-apply(combined.avg,1,max)
combined.norm<-sweep(combined.avg,1,combined.colmax,FUN="/")


bs.colmax<-apply(bs.avg,1,max)
bs.norm<-sweep(bs.avg,1,bs.colmax,FUN="/")


combined_sc <- subset(combined, subset = old.ident == 'Spinalcord Cells')
combined_sc<-SetIdent(combined_sc, value = combined_sc[['seurat_clusters']])
sc.avg<-as.data.frame(AverageExpression(combined_sc, features = marker_genes, group.by = 'seurat_clusters', slot = 'data', assays = 'RNA'))
# do the same for spinal cord data
sc.avg<-sc.avg[common.clusters]
sc.colmax<-apply(sc.avg,1,max)
sc.norm<-sweep(sc.avg,1,sc.colmax,FUN="/")


# function to get the majority of cells present in each seurat cluster based on lamina position to add annotation labels 
lamina.label<- function(object){
df.meta<-subset(object@meta.data,select = c('Laminae_fam','seurat_clusters'))
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
label.table$Laminae_fam <- factor(label.table$Laminae_fam, levels = c("1/2o","1/2o/2i","2/3","3","3/4","4/5","5","6","7","8","9","10"))
return(label.table)
}

############## HEATMAP for combined dataset

combined.labels <- lamina.label(combined) # getting lamina labels for each cluster
# remove label that is not needed
extra<-setdiff(colnames(combined.avg.all), colnames(bs.avg))
combined.labels<- combined.labels[!(row.names(combined.labels) %in% extra),]

# adding annotation for laminae 
my_colour = list(lamina = c("1/2o" = "red","1/2o/2i" = "red","2/3" = "red","3" = "red","3/4" = "blue","4/5" = "blue","5" = "blue","6" = "blue","7"="green","8"= "green","9" = "green","10"="green"))

# plotting the heat map
col_fun = colorRamp2(c(0, 1), c("black", "red"))
col_fun(seq(-3, 3)) # palette colors 

# lamina annotation for clusters 
ha = HeatmapAnnotation(lamina = combined.labels,
                       col = my_colour)

# getting the row and column order so that it can be applied to other heatmaps 
o1 = seriate(dist(combined.norm), method = "GW")
o2 = seriate(dist(t(combined.norm)), method = "GW")

col_1<-get_order(o2)

# by default, hclust is applied 
p0<-Heatmap(combined.norm, cluster_rows = F, col = col_fun, top_annotation = ha, name = 'combined', column_order = col_1)
p0


############## HEATMAP for brainstem data 

# plotting the heatmap
p1<-Heatmap(bs.norm, cluster_rows = F, col = col_fun, column_order = col_1, top_annotation = ha, name = 'brainstem')
p1


############## HEATMAP for spinal cord data


# plotting the heatmap
p2<-Heatmap(sc.norm, col = col_fun, column_order = col_1, cluster_rows = F, top_annotation = ha, name = 'spinal cord')
p2



######################### plotting cells in a cluster to find out your type 
# highlight some cells
int.cells<- WhichCells(combined, ident = c(3,50))
p1<- DimPlot(combined, reduction = "umap", cells.highlight = int.cells, label = T) + NoLegend()

# plot final cluster assignments
colors.needed<- CustomColor(combined, cluster.by = 'final_cluster_assignment', seed = T)
p2<- DimPlot(combined, group.by = 'final_cluster_assignment',label = T, cols = colors.needed, label.size = 3.5, repel = TRUE,shuffle = TRUE)

p1+p2


