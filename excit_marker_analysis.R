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

# Reading in the excit clustered dataset 
combined_excit<-readRDS('method1/excit_integrate1_clustered.rds')

combined_excit<-SetIdent(combined_excit, value = combined_excit[['seurat_clusters']])
######### Hierarchial clustering #########
combined_excit<-BuildClusterTree(combined_excit,assay="integrated", dims = 1:40)
PlotClusterTree(combined_excit, direction = "downwards")

tree <- Tool(combined_excit, slot = 'BuildClusterTree') # This extracts the phylo tree from our object


######### Group 1 -  node 66 - labels Excit 8,9,10 #########
leaves<-GetLeavesOnly(combined_excit, tree, 66)
node66.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
               ident.1 = leaves, 
               min.pct = 0.6,
               test.use = 'wilcox',
               min.diff.pct = 0.3, 
               only.pos = TRUE)


node66.markers.top10table <- node66.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node66.markers.top10<-rownames(node66.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Reln', 'Pde11a', 'Prex2', '9530026P05Rik'), cols = c('lightgreen', 'red'), label = T)
p1+p2



######### Group 2 -  node 70 - labels Excit 17,18,19 #########
leaves<-GetLeavesOnly(tree, 70)
node70.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node70.markers.top10table <- node70.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node70.markers.top10<-rownames(node70.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Nmu','Tac2','Sox5'), cols = c('lightgreen', 'red'), label = T)
p1+p2


######### Group 3 -  node 74 - labels Excit 14,15,16 #########
leaves<-GetLeavesOnly(combined_excit,tree, 74)
node74.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
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
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Lmo3','Grik1'), cols = c('lightgreen', 'red'), label = T, min.cutoff = 0)
p1+p2



######### Group 4 -  node 75 - labels Excit 11,12,13 #########
leaves<-GetLeavesOnly(tree, 75)
node75.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
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
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Rreb1','Rgs6','Man1a'), cols = c('lightgreen', 'red'), label = T)
p1+p2
p2



######### Group 5 -  node 64 - labels Excit 1 #########
leaves<-GetLeavesOnly(combined_excit, tree, 64)
node64.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
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
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Cpne4','Bnc2','Dach2','Erbb4'), cols = c('lightgreen', 'red'), label = T)
p1+p2
p2



######### Group 6 -  node 76 - labels Excit 4 #########
leaves<-GetLeavesOnly(tree, 76)
node76.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node76.markers.top10table <- node76.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node76.markers.top10<-rownames(node76.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Nts','Adamtsl1','Zfpm2','Col25a1'), cols = c('lightgreen', 'red'), label = T)
p1+p2
p2


######### Group 7 -  node 77 - labels Excit 2,3 #########
leaves<-GetLeavesOnly(tree, 77)
node77.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node77.markers.top10table <- node77.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node77.markers.top10<-rownames(node77.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Cck','Prkcg'), cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2
p2


######### Group 8 -  Node 72 - labels Excit 20 and many others #########
leaves<-GetLeavesOnly(tree, 72)
node72.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node72.markers.top10table <- node72.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node72.markers.top10<-rownames(node72.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Htr2c','Slit2','Meis2','Esrrg'), cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2
p2




######### Group 9 -  Node 73 - labels Excit 5,6,7 #########
leaves<-GetLeavesOnly(combined_excit,tree, 73)
node73.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                            ident.1 = leaves, 
                            min.pct = 0.2, 
                            logfc.threshold = 0.6,
                            test.use = 'wilcox',
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)


node73.markers.top10table <- node73.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
node73.markers.top10<-rownames(node73.markers.top10table)

# Comparison plot of all types 
options(ggrepel.max.overlaps = 20)
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, cells.highlight = cells.int, group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Maf','Adarb2','Elmo1'), cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2




######### Group 10 -  Node 68 - labels Excit 21,22 #########
leaves<-GetLeavesOnly(tree, 68)
node68.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
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
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE,  cells.highlight = cells.int,group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = node68.markers.top10[1:5], cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2
p2



######### Group 11 -  Node 69 - labels Excit 25,  #########
leaves<-GetLeavesOnly(tree, 69)
node69.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
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
colors.needed<-CustomColor(combined_excit,cluster.by='Final.clusters', seed = T)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE, group.by = 'Final.clusters', cols = colors.needed)
p1

# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = leaves)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE,  cells.highlight = cells.int,group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = c('Tcf4','Nfia'), cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2
p2




################################################################################
# Find markers for individual seurat clusters
clust49.markers<-FindMarkers(combined_excit,          # no need of max cells ident, max is only 1000 in some types
                             ident.1 = 49, 
                             min.pct = 0.2, 
                             logfc.threshold = 0.6,
                             test.use = 'wilcox',
                             min.diff.pct = 0.1, 
                             only.pos = TRUE)

clust49.markers.top10table <- clust49.markers %>% top_n(n = 10, wt = (avg_log2FC * (pct.1/pct.2)))
clust49.markers.top10<-rownames(clust49.markers.top10table)



# Highlight the cells of the cluster you're looking at 
cells.int<-WhichCells(combined_excit, idents = 49)
p1<-DimPlot(combined_excit, reduction = "umap", label = TRUE, repel = TRUE,shuffle = TRUE,  cells.highlight = cells.int,group.by = 'Final.clusters') + NoLegend()
p1

p2<-FeaturePlot(combined_excit, features = clust49.markers.top10[1:4], cols = c('lightgreen', 'red'), label = T, min.cutoff = c(0))
p1+p2
p2



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

