#########################################################################
## MARKER ANALYSIS FOR THE INHIB DATASET 
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
setwd('Y:/PD/spinal_cord/levine')
# on the outside drive
setwd('U:/eng_research_economo/PD/spinal_cord/levine/')

# load in functions file
# inside drive
source('Y:/PD/spinal_cord/levine/scripts/functions.R')
# outside drive
source('U:/eng_research_economo/PD/spinal_cord/scripts/functions.R')


# source the tree seurat functions
# inside drive
source('Y:/PD/spinal_cord/levine/scripts/tree_functions_seurat.R')

#########################################################################

# Reading in the inhib clustered dataset 

combined_inhib<-readRDS('method1/inhib_integrate1_clustered.rds')

combined_inhib<-SetIdent(combined_inhib, value = combined_inhib[['seurat_clusters']])
######### Hierarchial clustering #########
combined_inhib<-BuildClusterTree(combined_inhib,assay="integrated", dims = 1:40)
PlotClusterTree(combined_inhib, direction = "downwards")

tree <- Tool(combined_inhib, slot = 'BuildClusterTree') # This extracts the phylo tree from our object


######### Group 1 -  node 69 - labels Inhib 1 #########
leaves<-GetLeavesOnly(combined_inhib,tree, 69)
node69.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node69.markers.top10table <- node69.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node69.markers.top10<-rownames(node69.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, features = node69.markers.top10[1:3], cols = c('lightgreen', 'red'), label = T)
p1+p2


######### Group 2 -  node 74 - labels Inhib 6,7 (and 10,11) #########
leaves<-GetLeavesOnly(combined_inhib,tree, 74)
node74.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node74.markers.top10table <- node74.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node74.markers.top10<-rownames(node74.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, features = node74.markers.top10[1:3], cols = c('lightgreen', 'red'), label = T)
p1+p2




######### Group 3 -  node 75 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 75)
node75.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node75.markers.top10table <- node75.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node75.markers.top10<-rownames(node75.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, node75.markers.top10[1:4], cols = c('lightgreen', 'red'), label = T)
p1+p2



######### Group 4 -  node 67 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 67)
node67.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node67.markers.top10table <- node67.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node67.markers.top10<-rownames(node67.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, node67.markers.top10[1:4], cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2


######### Group 5 -  node 68 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 68)
node68.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node68.markers.top10table <- node68.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node68.markers.top10<-rownames(node68.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, node68.markers.top10[2], cols = c('lightgreen', 'red'), label = T)
p1+p2



######### Group 6 -  node 59 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 59)
node59.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node59.markers.top10table <- node59.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node59.markers.top10<-rownames(node59.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, node59.markers.top10[1:4], cols = c('lightgreen', 'red'), label = T)
p1+p2





######### Group 7-  node 63 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 63)
node63.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node63.markers.top10table <- node63.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node63.markers.top10<-rownames(node63.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, node63.markers.top10[1:4], cols = c('lightgreen', 'red'), label = T)
p1+p2




######### Group 8-  node 64 - labels Inhib #########
leaves<-GetLeavesOnly(combined_inhib,tree, 64)
node64.markers<-FindMarkers(combined_inhib,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node64.markers.top10table <- node64.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node64.markers.top10<-rownames(node64.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_inhib,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_inhib, idents = leaves)
p1<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_inhib, c('Gal','Pde11a'), cols = c('lightgreen', 'red'), label = T)
p1+p2


##### cluster 51 #####
p2<-FeaturePlot(combined_inhib, c('Chat'), cols = c('lightgreen', 'red'), label = T)
p2


######### Final annotations of the cluster families ########################################################################

combined_inhib <- RenameIdents(combined_inhib,"53"="Adarb2","51"="Adarb2_Chat","25"="Adarb2","47"="Rorb","18"="Rorb","48"="Rorb","1"="Chrm3","41"="Chrm3","22"="Rorb","0"="Rorb","7"="Rorb","23"="Adarb2","44"="Adarb2","46"="Adarb2","26"="Adarb2","36"="Adarb2","29"="Adarb2","16"="Adarb2","28"="Adarb2","39"="Adarb2","21"="Adarb2","5"="Slit2","43"="Adarb2","33"="Slit2","40"="Slit2","11"="Slit2","32"="Slit2","10"="Rorb","31"="Slit2","45"="Slit2","3"="Chrm3","14"="Chrm3","2"="Chrm3","8"="Chrm3","9"="Chrm3","50"="Adarb2","6"="Adarb2","35"="Adarb2","38"="Slit2","42"="Slit2","27"="Gal","12"="Npy","52"="Adarb2","15"="Npy","34"="Npy","37"="Npy","17"="Rorb_Chrm3","30"="Rorb_Chrm3","13"="Rorb_Chrm3","19"="Rorb_Chrm3","20"="Adarb2_Gal","24"="Rorb_Adarb2_Gal","4"="Rorb_Adarb2_Gal","49"="Adarb2")
combined_inhib <- StashIdent(combined_inhib, save.name = 'Fam.label')

colors.needed<-CustomColor(combined_inhib,seed = T)
p1<-DimPlot(combined_inhib, reduction = "umap", label = F, repel = TRUE,shuffle = TRUE, cols = colors.needed) 
p1



################# Get the cell composition in a seurat cluster################################################################
combined_inhib<-SetIdent(combined_inhib, value = combined_inhib[['seurat_clusters']])
cells.int<-WhichCells(combined_inhib, idents = c(49))
data.int<-combined_inhib[,cells.int]
table(data.int[['Final.clusters']])






# To highlight the cells belonging to a particular ident, and seeing its distribution across seurat_clusters
combined_inhib<-SetIdent(combined_inhib, value = combined_inhib[['Final.clusters']])
cells.int<-WhichCells(combined_inhib, idents = c('Inhib.23'))
combined_inhib<-SetIdent(combined_inhib, value = combined_inhib[['Final.clusters']])
p3<-DimPlot(combined_inhib, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int) + NoLegend()
p3

combined_inhib<-SetIdent(combined_inhib, value = combined_inhib[['seurat_clusters']])




