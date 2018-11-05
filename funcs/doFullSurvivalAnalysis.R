# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Performs survival analysis using the input features
# Cite        : TBD
# ******************************************************************************

doFullSurvivalAnalysis = function(clustered_list,subjects,num_folds,stop_limit,num_repeat,top_feat_thresh,seed_final,scolors,group_colors,lcolors_original,set_plots,dir_surv){
  
  feats_grouped = data.frame(clustered_list$feats_grouped)
  num_feats = clustered_list$num_feats
  
  # -------------
  # Repeated cross-validation
  # -------------
  
  res_repeat = matrix(0,nrow=num_repeat,ncol=2)
  colnames(res_repeat) = c('cindex_tr','cindex')
  
  feat_counts = numeric(num_feats)
  names(feat_counts) = colnames(feats_grouped)
  
  feat_betas = numeric(num_feats)
  names(feat_betas) = colnames(feats_grouped)
  
  for(repeatIdx in 1:num_repeat){
    
    dir_repeat = paste(dir_surv,"repeat_",repeatIdx,"/",sep="")
    
    if(set_plots){
      
      dir.create(dir_repeat)
      
    }
    
    
    cat(rep("#",30),"\n",sep="")
    cat('Doing repeat test : ',repeatIdx,'\n')
    
    model_rpt = coxSurvivalWithBackSearch(clustered_list,subjects,group_colors,num_folds,stop_limit,set_plots,dir_repeat)
    
    res_repeat[repeatIdx,'cindex_tr'] = model_rpt$cindex_train
    res_repeat[repeatIdx,'cindex'] = model_rpt$cindex_test
    
    feat_counts = feat_counts + colSums(model_rpt$feat_counts,na.rm=TRUE)
    feat_betas = feat_betas + colSums(model_rpt$feat_betas,na.rm=TRUE)
    
  }
  
  write.csv(res_repeat,file=paste(dir_surv,"res_repeat.csv",sep=""),row.names=FALSE,na="")
  
  # -------------
  # Train and test set performance in one plot
  # -------------
  
  train_test = cbind(res_repeat[,c('cindex_tr','cindex')])
  colnames(train_test) = c('train','test')
  write.csv(train_test,file=paste(dir_surv,'train_test.csv',sep=""),row.names = T)
  
  df_test = as.data.frame(cbind(c(rep(1,num_repeat),rep(2,num_repeat)),as.vector(as.matrix(train_test))))
  colnames(df_test) = c('label','C_index')
  df_test$label = as.factor(df_test$label)
  
  pdf(paste(dir_surv,'train_test.pdf',sep=""))
  p = ggplot(df_test,aes(x=label,y=C_index)) + geom_boxplot(width=0.6,notch = F,outlier.shape = NA, na.rm=T, size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Train','Test')) + scale_y_continuous(limits = c(0.5,1), breaks=seq(0.5,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3','deeppink4')) + xlab("") + ylab('Concordance Index\n')
  p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
  p = p + geom_hline(yintercept=median(train_test[,'train'],na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
  p = p + geom_hline(yintercept=median(train_test[,'test'],na.rm=T),colour='deeppink4',size=1.2,linetype='dashed',alpha=0.7)
  print(p)
  dev.off()
  
  # -------------
  # Select the most frequent features
  # -------------
  
  feat_percent = (feat_counts*100)/(num_folds*num_repeat)
  
  betas_avg = feat_betas / feat_counts
  
  feats_order = order(feat_percent,decreasing = T)
  feat_percent = feat_percent[feats_order]
  
  pdf(paste(dir_surv,'cox_feat_freq.pdf',sep=""))
  par(mar=c(12,4,2,0.5))
  barplot(feat_percent[1:10],main=c(" Mean frequencies of features "),ylab='Percent',las=2,cex.names=0.6)
  abline(h = top_feat_thresh,lwd=2,lty=2)
  dev.off()
  
  top_feat_idx = feats_order[which(feat_percent>top_feat_thresh)]
  
  if(length(top_feat_idx)<2){
    
    top_feat_order = order(feat_percent,decreasing = T)
    top_feat_idx = top_feat_order[c(1,2)]
    
  }
  
  feat_counts = feat_counts[feats_order]
  write.csv(feat_counts,paste(dir_surv,'feat_counts.csv',sep=""),row.names = T,col.names = T)
  betas_avg = betas_avg[feats_order]
  write.csv(betas_avg,paste(dir_surv,'feat_betas.csv',sep=""),row.names = T,col.names = T)
  
  pdf(paste(dir_surv,'counts_stats.pdf',sep=""))
  par(mar=c(15,4,2,0.5))
  barplot(feat_percent[1:4],main=c(" Predictor Counts"),beside=TRUE,las=2,cex.names=0.9,ylim=c(0,100))
  dev.off()
  
  pdf(paste(dir_surv,'coeffs_stats.pdf',sep=""))
  par(mar=c(15,4,2,0.5))
  barplot(betas_avg[1:4],main=c(" Predictor Coefficients"),beside=TRUE,las=2,cex.names=0.9,ylim=c(-1,1))
  dev.off()
  
  # -------------
  # Build final model with most frequent features
  # -------------
  
  dir_final = paste(dir_surv,"final/",sep="")
  dir.create(dir_final)
  
  model_final = coxSurvivalFinalModel(feats_grouped[,top_feat_idx],subjects,seed_final,scolors,group_colors,num_folds,plots=TRUE,dir_final)
  cindex_final = model_final$cindex_test
  
  feats_top = feats_grouped[,top_feat_idx]
  write.csv(feats_top,file=paste(dir_final,'feats.csv',sep=""),row.names = T,col.names = T)
  
  return(list('cindex_repeat'=res_repeat[,'cindex'],'cindex_final'=cindex_final,'top_feat_idx'=top_feat_idx))
  
}