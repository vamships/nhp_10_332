# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Selects features based on correlation with outcome using training set
#               Also minimizes pair-wise correlation between selected features
# Cite        : TBD
# ******************************************************************************

selectFromClusters = function(clustered_list,subj_original,train_idx){
  
  num_clusters = length(clustered_list)
  num_samples = nrow(clustered_list[[1]])
  
  feats_sel = matrix(NA,nrow=num_samples,ncol=num_clusters)
  feats_sel_names = NULL
  
  running_idx = 0
  selected_feat_idx = numeric(num_clusters)
  
  sel_cor = numeric(num_clusters)
  
  for(clusterIdx in 1:num_clusters){
    
    feats_in_cluster = clustered_list[[clusterIdx]]
    num_feats_in_cluster = ncol(feats_in_cluster)
    
    xcor = numeric(num_feats_in_cluster)
    
    for(featIdx in 1:num_feats_in_cluster){
      
      xcor[featIdx] = polyserial(feats_in_cluster[train_idx,featIdx],subj_original[train_idx,'Challenges'])
      
    }
    
    xcor_order = order(abs(xcor),decreasing = T)
    sel_idx = xcor_order[1]
    
    sel_cor[clusterIdx] = abs(xcor[sel_idx])
    
    feats_sel[,clusterIdx] = feats_in_cluster[,sel_idx]
    feats_sel_names = c(feats_sel_names,colnames(feats_in_cluster)[[sel_idx]])
    
    selected_feat_idx[clusterIdx] = sel_idx + running_idx
    running_idx = running_idx + num_feats_in_cluster
    
  }
  
  colnames(feats_sel) = feats_sel_names
  
  # order feats by polyserial
  xcor_order = order(sel_cor)
  feats_sel = feats_sel[,xcor_order]
  selected_feat_idx = selected_feat_idx[xcor_order]

  mark_sel = matrix(0,nrow=1,ncol=ncol(feats_sel))
  rownames(mark_sel) = 'selection'
  colnames(mark_sel) = colnames(feats_sel)

  # minimize pair-wise correlations
  for(featIdx in 1:ncol(feats_sel)){

    if(mark_sel[1,featIdx]==0){

      # Mark as selected
      mark_sel[1,featIdx] = 1

      # Feats
      idx_to_consider = which(mark_sel==0)

      if(length(idx_to_consider)>0){

        cor_consider = cor(feats_sel[train_idx,featIdx],feats_sel[train_idx,idx_to_consider],use='pairwise.complete.obs')
        idx_throw = which(cor_consider>0.7)

        if(length(idx_throw)>0){

          mark_sel[1,idx_to_consider[idx_throw]] = 2

        }

      }

    }

  }

  sub_feat_idx_ovrl = which(mark_sel==1)

  feats_sel = feats_sel[,sub_feat_idx_ovrl]
  selected_feat_idx = selected_feat_idx[sub_feat_idx_ovrl]
  
  return(list('feats_sel'=feats_sel,'selected_feat_idx'=selected_feat_idx))
  
}