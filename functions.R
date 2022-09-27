###### Max normalization FOR HEATMAP IF SELECTED MARKERS TO DETERMINE CLUSTERS - BASED ON NEUROTRANSMITTER STATUS
MinMaxNorm <- function(x, na.rm = TRUE) {
  return((x-min(x))/(max(x)-min(x)))
}

MaxNorm<-function(x, na.rm = TRUE) {
  return(x/max(x))
}

ZScore<-function(x, na.rm = TRUE){
  return(x-mean(x)/sd(x))
}


NormalizeMarkers <- function(object, marker_genes, method = 'max',cluster_by = 'cluster.id', means = T) {
  
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
  
  if (method == 'min_max') {
    
    if (!means) {
      
      object_means <- object_frame
      
      min_max_means <- apply(object_means[,1:ncol(object_means)-1], 2, MinMaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      
      return(min_max_means)
      
    }
    
    else {
      
      object_means <- object_frame %>% 
        group_by(cluster) %>%
        summarize(across(everything(),mean))
      
      min_max_means <- apply(object_means[,2:ncol(object_means)], 2, MinMaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      
      return(min_max_means)
      
    }
    
  }
  
  else if (method == 'max') {
    
    if (!means) {
      
      object_means <- object_frame
      
      min_max_means <- apply(object_means[,1:ncol(object_means)-1], 2, MaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      
      return(min_max_means)
      
    }
    
    else {
      
      object_means <- object_frame %>% 
        group_by(cluster) %>%
        summarize(across(everything(),mean))
      
      min_max_means <- apply(object_means[,2:ncol(object_means)], 2, MaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      min_max_means <- min_max_means[order(as.numeric(min_max_means$cluster)), ]
      return(min_max_means)
      
    }
    
  }
  
  
  
  
  else if (method == 'max_value') {
    
    if (!means) {
      
      object_means <- object_frame
      
      min_max_means <- apply(object_means[,1:ncol(object_means)-1], 2, MaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      
      return(min_max_means)
      
    }
    
    else {
      
      object_means <- object_frame %>% 
        group_by(cluster) %>%
        summarize(across(everything(),mean))
      
      min_max_means <- apply(object_means[,2:ncol(object_means)], 2, MaxNorm)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      min_max_means <- min_max_means[order(as.numeric(min_max_means$cluster)), ]
      return(min_max_means)
      
    }
    
  }
  
  
  
  
  
  else if (method == 'zscore') {
    if (!means) {
      
      object_means <- object_frame
      
      min_max_means <- apply(object_means[,1:ncol(object_means)-1], 2, ZScore)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      
      return(min_max_means)
      
    }
    
    else {
      
      object_means <- object_frame %>% 
        group_by(cluster) %>%
        summarize(across(everything(),mean))
      
      min_max_means <- apply(object_means[,2:ncol(object_means)], 2, ZScore)
      min_max_means <- data.frame(min_max_means, cluster = object_means$cluster)
      min_max_means <- min_max_means[order(as.numeric(min_max_means$cluster)), ]
      return(min_max_means)
      
    }

  }
  
  
  
}

construct.heatmap.maxnorm<-function(cluster_norm_expression, title = "Heatmap of cell markers"){
  og<-cluster_norm_expression  
  cluster_norm_expression<-cluster_norm_expression[1:(length(cluster_norm_expression)-1)]  
  cluster_norm_expression[is.na(cluster_norm_expression)] <- 0
  breaks = seq(min(cluster_norm_expression),max(cluster_norm_expression),length.out=1000)
  
  heatmap.2(as.matrix(t(cluster_norm_expression)),
            main = title, # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            #margins =c(1,11),     # widens margins around plot
            dendrogram='both',     # only draw a row dendrogram
            Colv=F,
            Rowv=F,  
            key=T,
            col=colorRampPalette(c('white','red','dark red')),
            labCol = og[,ncol(og)],  
            breaks=breaks,
            symm=F,
            symbreaks=F,
            cexCol=1,
            srtCol=70,
            margins=c(12,12))
  
}




CustomColor<-function(object, cluster.by = 'active.ident', seed = F) {
  # if you want a seed
  if(seed==T){
    set.seed(1)
  }
  
  if(cluster.by == 'active.ident'){
    x<-levels(object)
  }
  
  else{
    x<-unique(object[[cluster.by]])[[cluster.by]]
  }
  
  colors <- list(cell_type = sample(createPalette(length(x), c('#FF0000', '#00FF00', '#0000FF'), length(x)), length(x), replace = FALSE))
  names(colors$cell_type) <- c(x)
  colors.sort <-colors$cell_type[sort(names(colors$cell_type))]
  return(colors.sort)
}
