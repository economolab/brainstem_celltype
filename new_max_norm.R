
object<-combined_inhib
cluster_by = 'seurat_clusters'


object_data <- object@assays$RNA@counts
object_clusters <- as.vector(t(object[[cluster_by]]))


# checking if all genes exist in the object
marker_genes<-Reduce(intersect,list(marker_genes,rownames(object_data)))

object_matrix <- data.matrix(object_data[marker_genes[1],])



for (i in 2:length(marker_genes)) {
  object_matrix <- cbind(object_matrix, object_data[marker_genes[i],])
}

colnames(object_matrix) <- marker_genes

object_frame <- data.frame(object_matrix,
                           cluster = object_clusters)









object_means <- object_frame %>% 
  group_by(cluster) %>%
  summarize(across(everything(),mean))

min_max_means <- apply(object_means[,2:ncol(object_means)], 2, MaxNorm)
min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
min_max_means <- min_max_means[order(as.numeric(min_max_means$cluster)), ]
return(min_max_means)