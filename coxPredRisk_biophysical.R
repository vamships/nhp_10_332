# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Performs survival analysis using the biophysical measurements
# Cite        : TBD
# ******************************************************************************

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_coxPred_biophysical/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 10 # number of folds for cross-validation
num_clusters = 9 # number of clusters to divide the features into
stop_limit = 0.25 # log-likelihood cutoff for stopping backward search
num_repeat = 100 # number of repetitions of cross-validation
top_feat_thresh = 50 # Frequency cutoff to determine most-frequent features
seed_final = 8357 # seed to determine folds for final cross-validation model
doLookBack = T # To look for correlates of correlates after the final model
lookBackCorThresh = 0.75 # correlation cutoff to determine if features are correlated
set_plots = TRUE


log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nCox-model parameters','\n',sep="",file=log_file,append=T)
cat('Number of folds : ',num_folds,'\n',file=log_file,append=T)
cat('Number of clusters : ',num_clusters,'\n',file=log_file,append=T)
cat('Stop Limit for Backward Elimiation : ',stop_limit,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',num_repeat,'\n',file=log_file,append=T)
cat('Threshold for frequency of features : ',top_feat_thresh,'\n',file=log_file,append=T)
cat('Seed for final model : ',seed_final,'\n',file=log_file,append=T)
cat('Activate lookBack feature : ',doLookBack,'\n',file=log_file,append=T)
cat('Correlation threshold for looking back : ',lookBackCorThresh,'\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

subjects = read.csv('in/subjects.csv', header=TRUE, row.names=1)

# biophysical data
luminex_tp5 = read.csv('in/luminex_tp5.csv',header=T,row.names=1)

# Prepare subject and column colors
scolors = createSubjectColors(subjects,group_colors,challenge_colors)
lcolors_tp5 = createColumnColors(colnames(luminex_tp5),reagent_names,reagent_colors,antigen_names,antigen_colors)

# -------------------------------------------
# Sec 03: Legend
# -------------------------------------------

pdf(paste(dir_res,'legend.pdf',sep=""))
plot.new()
legend("bottomleft",legend=names(group_id),fill=group_colors,cex=0.66)
legend("bottom",legend=c(names(challenge_colors[1:12]),'UI'),fill=challenge_colors,cex=0.66)
legend("bottomright",legend=c('timepoint 5','timepoint 6'),fill=tp_colors,cex=0.66)
legend("center",legend=names(ak_colors[c(2,4,1,3)]),fill=ak_colors[c(2,4,1,3)],cex=0.66)
legend("topright",legend=names(reagent_colors),fill=reagent_colors,cex=0.66)
legend("top",legend=names(antigen_colors_broad),fill=antigen_colors_broad,cex=0.66)
legend("topleft",legend=names(func_colors),fill=func_colors,cex=0.66)
dev.off()

# -------------------------------------------
# Sec 04: Feature Prefiltering
# -------------------------------------------

# Standardize
# if subjects have NA in features, set those features to 0
na_idx = which(is.na(luminex_tp5),arr.ind=TRUE)
feats_scaled = scale(data.matrix(luminex_tp5),center=TRUE,scale=TRUE)
if(length(na_idx) !=0)
  feats_scaled[na_idx] = 0

dir_cluster = paste(dir_res,'cluster/',sep="")
dir.create(dir_cluster)

clustered_list = clusterFeats(feats_scaled,num_clusters,"ward.D2",dir_cluster)

# -------------------------------------------
# Sec 05: Survival Analysis
# -------------------------------------------

dir_surv = paste(dir_res,"surv/",sep="")
dir.create(dir_surv)

surv_performance = doFullSurvivalAnalysis(clustered_list,subjects,num_folds,stop_limit,num_repeat,top_feat_thresh,seed_final,scolors,group_colors,lcolors_original,set_plots,dir_surv)

actual_final = surv_performance$cindex_final
actual_repeat = surv_performance$cindex_repeat

# -------------
# Lookback at features for other correlates
# -------------

if(doLookBack){

  cat('\n\nLookBack feature active\n',file=log_file,append=T)
  dir_lookBack = paste(dir_res,'lookback/',sep="")
  dir.create(dir_lookBack)

  coxSurvivalLookBack(clustered_list$feats_grouped,luminex_tp5,surv_performance$top_feat_idx,lookBackCorThresh,actual_final,seed_final,subjects,scolors,lcolors_tp5,group_colors,num_folds,log_file,dir_lookBack)

}

# -------------------------------------------
# Sec 06: Permutation test
# -------------------------------------------

dir_perm = paste(dir_res,"perm/",sep="")
dir.create(dir_perm)

perm_repeat = matrix(NA,nrow=num_repeat,ncol=1)
colnames(perm_repeat) = 'cindex_test'

perm_log = paste(dir_perm,'perms.txt',sep="")
file.create(perm_log)

for(testIdx in 1:num_repeat){

  cat(rep("#",30),"\n",sep="")
  cat("Permutation test :",testIdx,"\n")
  
  dir_test = paste(dir_perm,"test_",testIdx,"/",sep="")

  # permute rows of outcome matrix
  sample_size = nrow(luminex_tp5)
  perm_order = sample(sample_size)
  subj_x = subjects[perm_order,]

  cat(testIdx,':',perm_order,'\n',file=perm_log,append=T)

  # -------------
  # Survival Analysis
  # -------------
  
  model_perm = coxSurvivalWithBackSearch(clustered_list,subj_x,group_colors,num_folds,stop_limit,FALSE,dir_test)

  perm_repeat[testIdx,'cindex_test'] = model_perm$cindex_test

}

# -------------------------------------------
# Sec 07: Comparing actual and permuted models' performance
# -------------------------------------------

cdf_perm = ecdf(perm_repeat[,'cindex_test'])
actual_avg = mean(actual_repeat)
pval_actual = 1-cdf_perm(actual_avg)
perm_repeat = data.frame(perm_repeat)

# -------------
# Percentile of actual final model in permuted models
# -------------

pdf(paste(dir_perm,'compare_holistic.pdf',sep=""))
p = ggplot(perm_repeat,aes(x=cindex_test)) + geom_density(aes(y=..scaled..)) + scale_x_continuous(limits=c(0.3,1)) + ggtitle('Distribution of Permuted Models\n') + ylab('Density (scaled)\n') + theme(plot.title = element_text(size=25), axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),aspect.ratio=1, legend.position='bottom') + xlab('\nC-index')
p = p + geom_vline(xintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=0.45,y=0.25,size=6,label='Random',color='black')

p = p + geom_vline(xintercept=actual_avg,colour='darkorange3',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=actual_avg+0.06,y=0.8,size=6,label=paste('Actual\n p : ',round(pval_actual,digits=2),sep=""),color='darkorange3')

p = p + geom_vline(xintercept=mean(perm_repeat[,'cindex_test'],na.rm=T),colour='deeppink4',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=mean(perm_repeat[,'cindex_test'],na.rm=T),y=0.5,size=6,label='Permuted\nMean',color='deeppink4')

print(p)
dev.off()

# -------------
# Comparing repeated cross-validation between actual and permuted
# -------------

df_test = as.data.frame(cbind(c(rep(1,num_repeat),rep(2,num_repeat)),as.vector(as.matrix(cbind(actual_repeat,perm_repeat[,'cindex_test'])))))
colnames(df_test) = c('label','C_index')
df_test$label = as.factor(df_test$label)

write.csv(df_test,file=paste(dir_perm,'robust.csv',sep=""),row.names = T)

diff_test = wilcox.test(actual_repeat,perm_repeat[,"cindex_test"],alternative="two.sided")

pdf(paste(dir_perm,'robust_test.pdf',sep=""))
p = ggplot(df_test,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0.25,1), breaks=seq(0.25,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=12,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('darkorange3','deeppink4')) + xlab("") + ylab('Concordance Index\n')
p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(actual_repeat,na.rm=T),colour='darkorange3',size=1.2,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(perm_repeat[,'cindex_test'],na.rm=T),colour='deeppink4',size=1.2,linetype='dashed',alpha=0.7)
p = p + geom_point(aes(x=1,y=actual_avg),shape=18,size=5,color='darkblue')
p = p + geom_hline(yintercept=actual_avg,colour='darkblue',size=0.7,linetype='dashed',alpha=0.7)
p = p + annotate("segment",x=1,xend=2,y=0.95,yend=0.95,size=2)
p = p + annotate("text",x=1.5,y=0.93,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T)))
print(p)
dev.off()

cat('\nModel Evaluation results','\n',file=log_file,append=T)
cat('Concordance Index\n',file=log_file,append=T)
cat('Repeated Cox : ',actual_final,'(Final model)',actual_avg,'(Test)','\n',file=log_file,append=T)
cat('Permuted Cox : ',mean(perm_repeat[,1],na.rm=T),'\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)