# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Performs survival analysis using the functional measurements
# Cite        : TBD
# ******************************************************************************

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_coxPred_functional/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 10 # number of folds for cross-validation
num_repeat = 100 # number of repetitions of cross-validation
seed_final = 8357 # seed to determine folds for final cross-validation model

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nCox-model parameters','\n',sep="",file=log_file,append=T)
cat('Number of folds : ',num_folds,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',num_repeat,'\n',file=log_file,append=T)
cat('Seed for final model',seed_final,'\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

subjects = read.csv('in/subjects.csv', header=TRUE, row.names=1)

# Functional Data
functions = read.csv('in/functions.csv', header=TRUE, row.names=1)

# Prepare subject and column colors
scolors = createSubjectColors(subjects,group_colors,challenge_colors)
funcolors = createFuncColors(colnames(functions),func_names,func_colors,tp_names,tp_colors)

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
# Sec 04: Feature Setup
# -------------------------------------------

feats = functions
lcolors_original = funcolors

# -------------------------------------------
# Sec 05: Survival Analysis
# -------------------------------------------

dir_surv = paste(dir_res,"surv/",sep="")
dir.create(dir_surv)

# -------------
# Visualize selected features
# -------------

ldata = scale(feats)
ldata[ldata>3] = 3; ldata[ldata < -3] = -3
lr = ceiling(max(abs(min(ldata,na.rm=TRUE)), max(ldata,na.rm=TRUE)))
lbreaks = seq(-lr,lr,0.1)

pdf(paste(dir_surv,'final-selection-raw.pdf',sep=""))
heatmap.4(ldata, col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.65, margin=c(10,8), breaks=lbreaks, symkey=FALSE, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors,NumRowSideColors=4,ColSideColors=lcolors_original,NumColSideColors=2,Rowv=FALSE, Colv=FALSE, na.color='black',lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3)), lhei=c(2,0.4,6.0),lwid=c(0.3,0.1,0.3))
dev.off()

pdf(paste(dir_surv,'final-selection-by-challenge.pdf',sep=""))
chall_sort = sort(subjects[,'Challenges'], index.return=TRUE)
heatmap.4(ldata[chall_sort$ix,], col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.6, margin=c(8,5), breaks=lbreaks, symkey=FALSE, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors[chall_sort$ix,],NumRowSideColors=4,ColSideColors=lcolors_original,NumColSideColors=2,Rowv=FALSE, Colv=FALSE, na.color='black',lmat=rbind(c(6,0,5),c(0,0,2),c(4,1,3)), lhei=c(2,0.4,6.0),lwid=c(0.3,0.1,0.3))
dev.off()

# -------------
# Repeated cross-validation
# -------------

cindex_repeat = matrix(NA,nrow=num_repeat,ncol=2)
colnames(cindex_repeat) = c('cindex_tr','cindex')

for(testIdx in 1:num_repeat){
  
  cat(rep("#",30),"\n",sep="")
  cat('Doing repeat test : ',testIdx,'\n')
  
  dir_repeat = paste(dir_surv,"repeat_",testIdx,"/",sep="")
  dir.create(dir_repeat)
  
  model_repeat = coxSurvivalFinalModel(feats,subjects,testIdx,scolors,group_colors,num_folds,plots=TRUE,dir_repeat)
  
  cindex_repeat[testIdx,1] = model_repeat$cindex_train
  cindex_repeat[testIdx,2] = model_repeat$cindex_test
  
}

# -------------
# Build final model
# -------------

dir_final = paste(dir_surv,"final/",sep="")
dir.create(dir_final)
model_final = coxSurvivalFinalModel(feats,subjects,seed_final,scolors,group_colors,num_folds,plots=TRUE,dir_final)

cindex_final = model_final$cindex_test

write.csv(feats,file=paste(dir_final,'feats.csv',sep=""),row.names = T,col.names = T)

# -------------
# Plot performance
# -------------

train_test = cindex_repeat
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

# -------------------------------------------
# Sec 06: Permutation test
# -------------------------------------------

dir_perm = paste(dir_res,"perm/",sep="")
dir.create(dir_perm)

cindex_perm = matrix(NA,nrow=num_repeat,ncol=1)
for(testIdx in 1:num_repeat){
  
  sample_size = nrow(feats)
  
  perm_order = sample(sample_size)
  subj_x = subjects[perm_order,]
  
  cat(rep("#",30),"\n",sep="")
  cat("Permutation test :",testIdx,"\n")
  
  dir_repeat = paste(dir_perm,"repeat_",testIdx,"/",sep="")
  dir.create(dir_repeat)
  
  model_perm = coxSurvivalFinalModel(feats,subj_x,testIdx,scolors,group_colors,num_folds,plots=TRUE,dir_repeat)
  
  cindex_perm[testIdx,] = model_perm$cindex_test
  
}

# -------------------------------------------
# Sec 07: Comparing actual and permuted models' performance
# -------------------------------------------

cdf_perm = ecdf(cindex_perm[,1])
actual_avg = mean(cindex_repeat[,2])
pval_actual = 1-cdf_perm(actual_avg)
cindex_perm = data.frame(cindex_perm)

# -------------
# Percentile of actual final model in permuted models
# -------------

pdf(paste(dir_res,'compare_holistic.pdf',sep=""))
p = ggplot(cindex_perm,aes(x=cindex_perm)) + geom_density(aes(y=..scaled..)) + scale_x_continuous(limits=c(0.3,1)) + ggtitle('Distribution of Permuted Models\n') + ylab('Density (scaled)\n') + theme(plot.title = element_text(size=25), axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),aspect.ratio=1, legend.position='bottom') + xlab('\nC-index')
p = p + geom_vline(xintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=0.45,y=0.25,size=6,label='Random',color='black')

p = p + geom_vline(xintercept=actual_avg,colour='darkorange3',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=actual_avg+0.06,y=0.8,size=6,label=paste('Actual\n p : ',round(pval_actual,digits=2),sep=""),color='darkorange3')

p = p + geom_vline(xintercept=mean(cindex_perm[,1],na.rm=T),colour='deeppink4',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("text",x=mean(cindex_perm[,1],na.rm=T),y=0.5,size=6,label='Permuted\nMean',color='deeppink4')

print(p)
dev.off()

# -------------
# Comparing repeated cross-validation between actual and permuted
# -------------

df_test = as.data.frame(cbind(c(rep(1,num_repeat),rep(2,num_repeat)),as.vector(as.matrix(cbind(cindex_repeat[,2],cindex_perm)))))
colnames(df_test) = c('label','C_index')
df_test$label = as.factor(df_test$label)

write.csv(df_test,file=paste(dir_res,'robust.csv',sep=""),row.names = T)

cindex_repeat = as.matrix(cindex_repeat)
diff_test = wilcox.test(cindex_repeat[,2],cindex_perm[,1],alternative="two.sided")
eff_test = cliff.delta(cindex_repeat[,2],cindex_perm[,1])
eff_interp = as.character(eff_test$magnitude)

pdf(paste(dir_res,'robust_test.pdf',sep=""))
p = ggplot(df_test,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted')) + scale_y_continuous(limits = c(0.25,0.85), breaks=seq(0.25,0.85,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=25,colour='black') ,axis.text.x = element_text(size=25,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='bottom') + scale_fill_manual(values=c('deepskyblue','mistyrose3')) + xlab("") + ylab('Concordance Index')
p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + annotate("segment",x=1,xend=2,y=0.75,yend=0.75,size=3,colour=effect_colors[eff_interp])
p = p + annotate("rect",xmin=c(1.43), xmax=c(1.57), ymin=c(0.73) , ymax=c(0.77), color=effect_colors[eff_interp], fill=effect_colors[eff_interp])
p = p + annotate("text",x=1.5,y=0.79,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T)))
print(p)
dev.off()

cat('\nModel Evaluation results','\n',file=log_file,append=T)
cat('Concordance Index\n',file=log_file,append=T)
cat('Repeated Cox : ',cindex_final,'(Final model)',mean(cindex_repeat[,2]),'(Test)','\n',file=log_file,append=T)
cat('Permuted Cox : ',mean(cindex_perm[,1],na.rm=T),'\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)