############ making heatmap ############ 


excit.object<-readRDS('levine_new/excit_integrate1_clustered.rds')
DefaultAssay(excit.object)<-'raw'
inhib.object<-readRDS('levine_new/inhib_integrate1_clustered.rds')
DefaultAssay(inhib.object)<-'raw'
cluster_by = 'seurat_clusters'


# checking if all genes exist in the object
E_markers <- c('Slc17a6','Slc17a7')
I_markers_gaba <- c('Gad2', 'Slc6a1','Slc32a1')
I_markers_glyc<-c('Slc6a5')
glia_markers<-c('Trem2','Aqp4')
marker_genes<-c(E_markers, I_markers_gaba, I_markers_glyc, glia_markers)


excit.avg<-AverageExpression(excit.object, features = marker_genes, group.by = 'seurat_clusters')
excit.avg<-excit.avg$raw
colnames(excit.avg)<-paste('E',colnames(excit.avg),sep='.')
# excit.colsums<- colSums(excit.avg) # sum up all the counts in a sample
# excit.norm<-sweep(excit.avg,2,excit.colsums,FUN="/")

excit.colmax<-apply(excit.avg,2,max)
excit.norm<-sweep(excit.avg,2,excit.colmax,FUN="/")

inhib.avg<-AverageExpression(inhib.object, features = marker_genes, group.by = 'seurat_clusters')
inhib.avg<-inhib.avg$raw
colnames(inhib.avg)<-paste('I',colnames(inhib.avg),sep='.')
# inhib.colsums<- colSums(inhib.avg) # sum up all the counts in a sample
# inhib.norm<-sweep(inhib.avg,2,inhib.colsums,FUN="/")

inhib.colmax<-apply(inhib.avg,2,max)
inhib.norm<-sweep(inhib.avg,2,inhib.colmax,FUN="/")

total<- cbind(excit.norm, inhib.norm)
plot.total<-


library(pheatmap)

pheatmap(plot.total, cluster_cols = F,cluster_rows = F, border_color = NA, color=colorRampPalette(c("black", "red"))(100))


pheatmap(combined.maxnorm, cluster_cols = F,cluster_rows = F, border_color = NA, color=colorRampPalette(c("navyblue", "white", "red"))(100))

pheatmap(excit.maxnorm, cluster_cols = F,cluster_rows = F, border_color = NA, color=colorRampPalette(c("navyblue", "white", "red"))(100))
pheatmap(inhib.maxnorm, cluster_cols = F,cluster_rows = F, border_color = NA, color=colorRampPalette(c("navyblue", "white", "red"))(100))


############ Hierarchial clustering ############ 












