# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Performs group classification using the biophysical measurements
# Cite        : TBD
# ******************************************************************************

rm(list = ls())

source('funcsToImport.R')

plot_font = 'Helvetica'

dir_res = paste('results_predClassM_biophysical/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

# glmnet parameters
alphas = 1
cvFolds = 10
repeatRun = 100
intc=TRUE
weights_bal = TRUE
grpType = "ungrouped"

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nGlmnet Multinomial Classification parameters','\n',sep="",file=log_file,append=T)
cat('Alpha range : ',paste(alphas,'',sep=','),'\n',file=log_file,append=T)
cat('Grouped regularization type : ',grpType,'\n',file=log_file,append=T)
cat('Balancing classes',weights_bal,'\n',file=log_file,append=T)
cat('Intercept for logistic regression',intc,'\n',file=log_file,append=T)
cat('Number of folds : ',cvFolds,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',repeatRun,'\n',file=log_file,append=T)

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

pdf(paste(dir_res,'legend.pdf',sep=""),family=plot_font)
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
# Sec 04: Glmnet classification
# -------------------------------------------

# Create labels for classification
classes = subjects[,'groupID']
feats = luminex_tp5

dir_class = paste(dir_res,"class/",sep="")
dir.create(dir_class)

numFeat = ncol(feats)

# -------------
# Visualize selected features
# -------------

ldata = scale(feats)
ldata[ldata>3] = 3; ldata[ldata < -3] = -3
lr = ceiling(max(abs(min(ldata,na.rm=TRUE)), max(ldata,na.rm=TRUE)))
lbreaks = seq(-lr,lr,0.1)

pdf(paste(dir_res,'final-selection-raw.pdf',sep=""))
heatmap.4(ldata, col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.6, margin=c(8,5), breaks=lbreaks, symkey=FALSE, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors,NumRowSideColors=4,Rowv=FALSE, Colv=FALSE, na.color='black')
dev.off()

pdf(paste(dir_res,'final-selection-by-challenge.pdf',sep=""))
chall_sort = sort(subjects[,'Challenges'], index.return=TRUE)
heatmap.4(ldata[chall_sort$ix,], col=bluered, scale='none', trace='none', cexRow=0.5, cexCol=0.6, margin=c(8,5), breaks=lbreaks, symkey=FALSE, dendrogram='none',hclust=hclust.ward,RowSideColors=scolors[chall_sort$ix,],NumRowSideColors=4,Rowv=FALSE, Colv=FALSE, na.color='black')
dev.off()

# -------------
# Scale features and set NAs to 0
# -------------

feats = scale(feats)
na_idx = which(is.na(feats),arr.ind=TRUE)
if(length(na_idx)!=0){
  
  feats[na_idx] = 0
  
}

cat('\n\nClassification Results ','\n',file=log_file,append=T)
cat('\nPredictions for groups classes\n',file=log_file,append=T)

label = classes

weights = rep(1,length(label))

class_model = glmnetMultiClass(feats,label,weights,numFeat,grpType,intc,alphas,cvFolds,repeatRun)

# -------------------------------------------
# Sec 05: Visualize prediction performance
# -------------------------------------------

pdf(paste(dir_class,'best_model.pdf',sep=""))
plot(class_model$final_fit,main=paste('alpha: ',class_model$best_alpha,'\n',sep=""))
dev.off()

# -------------
# Plot log-odds
# -------------

pred_prob = class_model$best_model$preval
pred_class = apply(pred_prob,1,which.max)

group_set = unique(label)
label_tform = numeric(length(pred_class))
for(chooseIdx in group_set){
  
  label_tform[label==chooseIdx] = chooseIdx-2
  
}

log_odds = numeric(length(pred_class))
for(sampleIdx in 1:length(pred_class)){
  
  pred_class[sampleIdx] = group_set[pred_class[sampleIdx]]
  other_prob = max(pred_prob[sampleIdx,-label_tform[sampleIdx]])
  log_odds[sampleIdx] = log2(pred_prob[sampleIdx,label_tform[sampleIdx]]/other_prob)
  
}

# -------------
# Plot confusion matrix
# -------------

confMat = confusionMatrix(pred_class,label)
plotConfusion(confMat,label,dir_class)

df_odds = as.data.frame(cbind(as.vector(log_odds),label,pred_class))
colnames(df_odds) = c('odds','label','label_pred')
df_odds$label = as.factor(df_odds$label)
df_odds$label_pred = as.factor(df_odds$label_pred)

pdf(paste(dir_class,'box_odds_gg.pdf',sep=""),family=plot_font)
p1 = ggplot(df_odds,aes(x=label,y=odds)) + geom_boxplot(width=0.5,notch = F,coef=1.58,outlier.shape = NA,size=1,colour="black") + geom_point(position = position_jitter(w=0.1),size=4,aes(colour=label_pred)) + ylab('Log Odds') + scale_x_discrete(name="",labels=names(group_id)) + theme(plot.title = element_text(size=20), axis.line = element_line(colour = "black",size=1.5), axis.title.y=element_text(size=20), axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none') + scale_colour_manual(values=group_colors) + scale_fill_manual(values=group_colors)
p1 = p1 + geom_hline(yintercept=0,colour='black',size=0.78,linetype='dashed',alpha=0.7)
print(p1)
dev.off()

# -------------
# Plot coefficients
# -------------

coeff_min = class_model$best_model$coeff
coeff_min_nz_idx = which(rowSums(abs(coeff_min))!=0)

pdf(paste(dir_class,'coeffs_min.pdf',sep=""))
par(mar=c(12,4,2,0.5))
yrange = range(coeff_min)
barplot(t(coeff_min[coeff_min_nz_idx,]),main=c(" Predictor Coefficients"),beside=TRUE,las=2,cex.names=0.9,col=group_colors)
dev.off()

coeffs_sum = rowSums(abs(coeff_min[coeff_min_nz_idx,]))
coeffs_order = order(coeffs_sum, decreasing=T)
coeffs_sel = coeff_min_nz_idx[coeffs_order[c(1,2)]]

feats_sel = luminex_tp5[,coeffs_sel]

pdf(paste(dir_class,"biplot.pdf",sep=""))
df = data.frame(cbind(feats_sel,subjects[,'groupID']))
colnames(df) = c('f1','f2','label')
df$label = as.factor(df$label)

p = ggplot(df,aes(x=f2,y=f1,colour=label)) + geom_point(size=4,aes(colour=label)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_color_manual(values=group_colors) + scale_fill_manual(values=group_colors) + xlab(colnames(feats_sel)[2]) + ylab(colnames(feats_sel)[1])
print(p)
dev.off()

# frequency of features in CV
feat_counts = rowSums(class_model$feat_counts)
#feat_order = order(feat_counts,decreasing = T)
#feat_counts = feat_counts[feat_order]
feat_betas = class_model$feat_betas
feat_betas = feat_betas / feat_counts
feat_stats = cbind(feat_counts,feat_betas)
rownames(feat_stats) = colnames(luminex_tp5)
colnames(feat_stats) = c('Percent','IM Mosaic','IM 239','AE 239')

write.csv(feat_stats,paste(dir_class,'feat_stats.csv',sep=""),row.names = T,col.names = T)

pdf(paste(dir_class,'counts_stats.pdf',sep=""))
par(mar=c(15,4,2,0.5))
yrange = range(coeff_min)
barplot(feat_stats[coeff_min_nz_idx,1],main=c(" Predictor Counts"),beside=TRUE,las=2,cex.names=0.9)
dev.off()

pdf(paste(dir_class,'coeffs_stats.pdf',sep=""))
par(mar=c(15,4,2,0.5))
yrange = range(coeff_min)
barplot(t(feat_stats[coeff_min_nz_idx,2:4]),main=c(" Predictor Coefficients"),beside=TRUE,las=2,cex.names=0.9,col=group_colors)
dev.off()

#
cat('Rept : Mean of Classification Error',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n')
cat('Perm : Mean of Classification Error',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n')

robust_test = as.data.frame(cbind(1-class_model$cv_repeat[,'min'],1-class_model$cv_permut[,'min']))
colnames(robust_test) = c('Luminex','Permuted')

df_test = as.data.frame(cbind(c(rep(1,repeatRun),rep(2,repeatRun)),as.vector(as.matrix(robust_test))))
colnames(df_test) = c('label','Acc')
df_test$label = as.factor(df_test$label)

cdf_perm = ecdf(robust_test[,2])
actual_avg = mean(robust_test[,1])
pval_actual = 1-cdf_perm(actual_avg)

# -------------------------------------------
# Sec 06: Size-matched permutation test
# -------------------------------------------

# Size matched robustness
cv_permut = matrix(numeric(1),repeatRun,2)
colnames(cv_permut) = c('min','se1')
randSelSize = length(coeff_min_nz_idx)
sample_size = ncol(feats)

for(testIdx in seq(1,repeatRun)){
  
  cv_perm = cv.glmnet(feats[,sample(sample_size,randSelSize)],label,nfolds=cvFolds,family="multinomial",standardize=FALSE,weights=weights,alpha=1, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
  
  cv_permut[testIdx,'min'] = cv_perm$cvm[match(cv_perm$lambda.min,cv_perm$lambda)]
  cv_permut[testIdx,'se1'] = cv_perm$cvm[match(cv_perm$lambda.1se,cv_perm$lambda)]
  
}

# -------------------------------------------
# Sec 07: Comparing actual and permuted models' performance
# -------------------------------------------

robust_randSel = as.data.frame(cbind(rep(3,repeatRun),cv_permut[,1]))
colnames(robust_randSel) = c('label','Acc')
robust_randSel$label = as.factor(robust_randSel$label)

robust = rbind(df_test,robust_randSel)
robust$label = as.factor(robust$label)

cdf_randSel = ecdf(robust_randSel[,2])
pval_actual_2 = 1-cdf_randSel(actual_avg)

eff_test_1 = cliff.delta(robust_test[,1],robust_test[,2])
eff_test_2 = cliff.delta(robust_test[,1],robust_randSel[,2])
eff_test_3 = cliff.delta(robust_test[,2],robust_randSel[,2])

eff_interp_1 = as.character(eff_test_1$magnitude)
eff_interp_2 = as.character(eff_test_2$magnitude)
eff_interp_3 = as.character(eff_test_3$magnitude)

pdf(paste(dir_class,'robust_test_3v.pdf',sep=""),family=plot_font)
p = ggplot(robust,aes(x=label,y=Acc)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted','Random\nSelection')) + scale_y_continuous(limits = c(0,1), breaks=seq(0,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_fill_manual(values=c('deepskyblue','mistyrose3','mistyrose3')) + xlab("") + ylab('Balanced Accuracy')

p = p + geom_hline(yintercept=0.33,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_test[,1],na.rm=T),colour='deepskyblue',size=1.2,linetype='dashed')
p = p + geom_hline(yintercept=mean(robust_test[,2],na.rm=T),colour='mistyrose3',size=1.2,linetype='dashed')
p = p + geom_hline(yintercept=mean(robust_randSel[,2],na.rm=T),colour='mistyrose3',size=1.2,linetype='dashed')

p = p + annotate("segment",x=1,xend=2,y=0.12,yend=0.12,size=2,colour=effect_colors[eff_interp_1])
p = p + annotate("rect",xmin=c(1.5), xmax=c(1.6), ymin=c(0.1) , ymax=c(0.14), color=effect_colors[eff_interp_1], fill=effect_colors[eff_interp_1])
p = p + annotate("text",x=1.5,y=0.17,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T)))

p = p + annotate("segment",x=1,xend=3,y=0.03,yend=0.03,size=2,colour=effect_colors[eff_interp_2])
p = p + annotate("rect",xmin=c(1.9), xmax=c(2), ymin=c(0.01) , ymax=c(0.05), color=effect_colors[eff_interp_2], fill=effect_colors[eff_interp_2])
p = p + annotate("text",x=2,y=0.08,size=6,label=paste('P : ',format(pval_actual_2,digits=3,scientific=T)))


print(p)
dev.off()

# print to log file
cat('Repeated Balanced Accuracy',class_model$repeat_mse_min,'(',class_model$repeat_mse_min_sd,')\n',file=log_file,append=T)
cat('Permuted Balanced Accuracy',class_model$permut_mse_min,'(',class_model$permut_mse_min_sd,')\n',file=log_file,append=T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
