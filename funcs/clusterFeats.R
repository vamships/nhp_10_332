# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Performs hierarchical clustering of features
# Cite        : TBD
# ******************************************************************************

clusterFeats = function(feats_scaled,num_clusters,method,dir_res){
  
  feats_dist = dist(t(feats_scaled))
  hier_clust = hclust(feats_dist,method=method)
  clusterCut = cutree(hier_clust,num_clusters)
  
  feats_in_cluster_list = vector("list",length=num_clusters)
  
  counts_in_cluster = numeric(num_clusters)
  
  feats_grouped = NULL
  lcolors_grouped = NULL
  
  for(clusterIdx in 1:num_clusters){
    
    feats_in_cluster_id = which(clusterCut==clusterIdx)
    feats_in_cluster = feats_scaled[,feats_in_cluster_id]
    
    counts_in_cluster[clusterIdx] = length(feats_in_cluster_id)
    feats_grouped = cbind(feats_grouped,feats_in_cluster)
    
    feats_in_cluster_list[[clusterIdx]] = feats_in_cluster
    
  }
  
  pdf(paste(dir_res,'cor_counts.pdf',sep=""))
  bp = barplot(counts_in_cluster,names.arg = 1:num_clusters)
  text(x=bp,y=counts_in_cluster+1,labels=counts_in_cluster,cex=1.5,xpd=TRUE)
  dev.off()
  
  return(list('cluster_list'=feats_in_cluster_list,'feats_grouped'=feats_grouped,'cluster_counts'=counts_in_cluster,'num_feats'=sum(counts_in_cluster)))
  
}