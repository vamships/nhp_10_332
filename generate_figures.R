# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Generates figures for the manuscript
# Cite        : TBD
# ******************************************************************************

# Copyright (C) <2018>  <Srivamshi Pittala>

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

rm(list = ls())

source('funcsToImport.R')

dir_res = paste('results_figures/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 10
num_repeat = 100 # number of repetitions of cross-validation
set_plots = TRUE

dir_biophysical = 'results_coxPred_biophysical/'
dir_network = 'results_networkPlot/'
dir_randomSel = 'results_coxPred_randomSel/'
dir_functional = 'results_coxPred_functional/'
dir_classify = 'results_predClassM_biophysical/'

plot_font = 'Helvetica'

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

# -------------------------------------------
# Sec 02: Data
# -------------------------------------------

subjects = read.csv('in/subjects.csv', header=TRUE, row.names=1)

# biophysical data
luminex_tp5 = read.csv('in/luminex_tp5.csv',header=T,row.names=1)

# Functional Data
functions = read.csv('in/functions.csv', header=TRUE, row.names=1)

# Prepare subject and column colors
scolors = createSubjectColors(subjects,group_colors,challenge_colors)
lcolors_tp5 = createColumnColors(colnames(luminex_tp5),reagent_names,reagent_colors,antigen_names,antigen_colors)
funcolors = createFuncColors(colnames(functions),func_names,func_colors,tp_names,tp_colors)

grp_3 = which(subjects[,'groupID']==3)
grp_4 = which(subjects[,'groupID']==4)
grp_5 = which(subjects[,'groupID']==5)

# -------------------------------------------
# Sec 03: Legend
# -------------------------------------------

dir_legends = paste(dir_res,'legends/',sep="")
dir.create(dir_legends)

pdf(paste(dir_legends,'groups.pdf',sep=""),family = plot_font)
plot.new()
legend("center",legend=names(group_id),fill=group_colors,cex=4)
dev.off()

pdf(paste(dir_legends,'ovp.pdf',sep=""),family = plot_font)
plot.new()
legend("center",legend=c('Predicted','Observed'),lty=c(1,2),cex=4,lwd=6)
dev.off()

pdf(paste(dir_legends,'main.pdf',sep=""),family = plot_font)
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
# Sec 04: Figure 2
# -------------------------------------------

dir_f2 = paste(dir_res,'fig_2/',sep="")
dir.create(dir_f2)

file.copy(paste(dir_biophysical,'lookback/1.R2A.4.high.C1.TR/26.R2A.4.high.G49/km_compare.pdf',sep=""),paste(dir_f2,'fig_f2a.pdf',sep=""))
file.copy(paste(dir_biophysical,'lookback/1.R2A.4.high.C1.TR/26.R2A.4.high.G49/CvR_compare_test_2.pdf',sep=""),paste(dir_f2,'fig_f2b.pdf',sep=""))
file.copy(paste(dir_biophysical,'lookback/1.R2A.4.high.C1.TR/26.R2A.4.high.G49/final_feat_selection_sorted.pdf',sep=""),paste(dir_f2,'fig_f2c.pdf',sep=""))

# -------------------------------------------
# Sec 05: Figure 3
# -------------------------------------------

dir_f3 = paste(dir_res,'fig_3/',sep="")
dir.create(dir_f3)

feats_sel = read.csv(paste(dir_functional,'surv/final/feats.csv',sep=""),row.names = 1)

for(featIdx in 1:ncol(feats_sel)){
  
  df = data.frame(cbind(feats_sel[,featIdx],subjects[,'groupID']))
  colnames(df) = c('feat','groupID')
  df$groupID = as.factor(df$groupID)
  
  featName = colnames(feats_sel)[featIdx]
  
  yrange = range(df$feat,na.rm = T)
  
  yrange[2] = yrange[2]*2
  ydiff = (yrange[2]-yrange[1])
  
  ystep = ydiff/8
  
  test_g3g4 = wilcox.test(df[grp_3,1],df[grp_4,1],alternative = 'two.sided')
  pval_g3g4 = test_g3g4$p.value
  
  test_g3g5 = wilcox.test(df[grp_3,1],df[grp_5,1],alternative = 'two.sided')
  pval_g3g5 = test_g3g5$p.value
  
  test_g4g5 = wilcox.test(df[grp_4,1],df[grp_5,1],alternative = 'two.sided')
  pval_g4g5 = test_g4g5$p.value
  
  pval_set = c(pval_g3g4,pval_g3g5,pval_g4g5)
  pval_set = format(pval_set,digits=3,scientific=T)
  
  pdf(paste(dir_f3,'fig_f3gh_',featIdx,'.pdf',sep=""))
  p1 = ggplot(df,aes(x=groupID,y=feat)) + geom_boxplot(width=0.4,notch = F,coef=0,outlier.shape = NA,size=1,colour="black") + geom_point(position = position_jitter(w=0.1),size=3.5,aes(colour=groupID)) + ylab('MFI') + ggtitle(featName) + scale_x_discrete(name="",labels=c('IM Mosaic','IM 239','AE 239')) + scale_y_continuous(limits=yrange) + theme(plot.title = element_text(size=15), axis.line = element_line(colour = "black",size=1.5), axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=20,colour='black'), axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none')  + scale_colour_manual(values=group_colors) + scale_fill_manual(values=group_colors)
  
  combIdx = 1
  p1 = p1 + annotate("segment",x=1,xend=2,y=yrange[1]+0.45*ydiff+ystep*combIdx,yend=yrange[1]+0.45*ydiff+ystep*combIdx,size=1,colour='black')
  #p1 = p1 + annotate("point",x=c(1,2),y=yrange[1]+0.45*ydiff+ystep*combIdx,size=4,colour='black')
  p1 = p1 + annotate("text",x=1.5,y=yrange[1]+0.48*ydiff+ystep*combIdx,size=5,label=paste('P : ',pval_set[combIdx]))
  
  combIdx = 2
  p1 = p1 + annotate("segment",x=1,xend=3,y=yrange[1]+0.45*ydiff+ystep*combIdx,yend=yrange[1]+0.45*ydiff+ystep*combIdx,size=1,colour='black')
  #p1 = p1 + annotate("point",x=c(1,3),y=yrange[1]+0.45*ydiff+ystep*combIdx,size=4,colour='black')
  p1 = p1 + annotate("text",x=2,y=yrange[1]+0.48*ydiff+ystep*combIdx,size=5,label=paste('P : ',pval_set[combIdx]))
  
  combIdx = 3
  p1 = p1 + annotate("segment",x=2,xend=3,y=yrange[1]+0.45*ydiff+ystep*combIdx,yend=yrange[1]+0.45*ydiff+ystep*combIdx,size=1,colour='black')
  #p1 = p1 + annotate("point",x=c(2,3),y=yrange[1]+0.45*ydiff+ystep*combIdx,size=4,colour='black')
  p1 = p1 + annotate("text",x=2.5,y=yrange[1]+0.48*ydiff+ystep*combIdx,size=5,label=paste('P : ',pval_set[combIdx]))
  
  print(p1)
  dev.off() 
  
}

file.copy(paste(dir_functional,'surv/final/km_compare.pdf',sep=""),paste(dir_f3,'fig_f3i.pdf',sep=""))
file.copy(paste(dir_functional,'surv/final/CvR_compare_test_2.pdf',sep=""),paste(dir_f3,'fig_f3j.pdf',sep=""))
file.copy(paste(dir_functional,'surv/final/final_feat_selection_sorted.pdf',sep=""),paste(dir_f3,'fig_f3k.pdf',sep=""))

# -------------------------------------------
# Sec 06: Figure S1
# -------------------------------------------

dir_s1 = paste(dir_res,'fig_s1/',sep="")
dir.create(dir_s1)

file.copy(paste(dir_classify,'class/confusion.pdf',sep=""),paste(dir_s1,'fig_s1c.pdf',sep=""))
file.copy(paste(dir_classify,'class/box_odds_gg.pdf',sep=""),paste(dir_s1,'fig_s1d.pdf',sep=""))
file.copy(paste(dir_classify,'class/robust_test_3v.pdf',sep=""),paste(dir_s1,'fig_s1e.pdf',sep=""))
file.copy(paste(dir_classify,'class/biplot.pdf',sep=""),paste(dir_s1,'fig_s1f.pdf',sep=""))

# -------------------------------------------
# Sec 07: Figure S2
# -------------------------------------------

dir_s2 = paste(dir_res,'fig_s2/',sep="")
dir.create(dir_s2)

file.copy(paste(dir_biophysical,'lookback/1.R2A.4.high.C1.TR/26.R2A.4.high.G49/CvR_test.pdf',sep=""),paste(dir_s2,'fig_s2a.pdf',sep=""))

robust_perm = read.csv(paste(dir_biophysical,'perm/robust.csv',sep=""),row.names = 1)
robust_randomSel = read.csv(paste(dir_randomSel,'perm/robust.csv',sep=""),row.names = 1)

robust_randomSel[,1] = 3
robust_all = as.data.frame(rbind(robust_perm,robust_randomSel))

robust_all$label = as.factor(robust_all$label)

diff_test_1 = wilcox.test(robust_all[1:100,2],robust_all[101:200,2],alternative="two.sided")
diff_test_2 = wilcox.test(robust_all[1:100,2],robust_all[201:300,2],alternative="two.sided")
diff_test_3 = wilcox.test(robust_all[101:200,2],robust_all[201:300,2],alternative="two.sided")

cdf_perm = ecdf(robust_all[101:200,2])
cdf_randSel = ecdf(robust_all[201:300,2])
actual_avg = mean(robust_all[1:100,2])
pval_actual = 1-cdf_perm(actual_avg)
pval_actual_2 = 1-cdf_randSel(actual_avg)

eff_test_1 = cliff.delta(robust_all[1:100,2],robust_all[101:200,2])
eff_test_2 = cliff.delta(robust_all[1:100,2],robust_all[201:300,2])
eff_test_3 = cliff.delta(robust_all[101:200,2],robust_all[201:300,2])

eff_interp_1 = as.character(eff_test_1$magnitude)
eff_interp_2 = as.character(eff_test_2$magnitude)
eff_interp_3 = as.character(eff_test_3$magnitude)

pdf(paste(dir_s2,'fig_s2b.pdf',sep=""),family = plot_font)
par(mar=c(12,7,3,1.5))
p = ggplot(robust_all,aes(x=label,y=C_index)) + geom_violin(size=1,colour="black",aes(fill=label)) + scale_x_discrete(labels=c('Actual','Permuted','Random\nSelection')) + scale_y_continuous(limits = c(0.25,1), breaks=seq(0.25,1,0.1)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=15,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=20,colour='black'), axis.text.y = element_text(size=17,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_fill_manual(values=c('deepskyblue','mistyrose3','mistyrose3')) + xlab("") + ylab('Concordance Index')

p = p + geom_hline(yintercept=0.5,colour='black',size=0.78,linetype='dashed',alpha=0.7)
p = p + geom_hline(yintercept=mean(robust_all[1:100,2],na.rm=T),colour='deepskyblue',size=1.5,linetype='dashed')
p = p + geom_hline(yintercept=mean(robust_all[201:300,2],na.rm=T),colour='mistyrose3',size=1.2,linetype='dashed')
p = p + geom_hline(yintercept=mean(robust_all[101:200,2],na.rm=T),colour='mistyrose3',size=1.2,linetype='dashed')

p = p + annotate("segment",x=1,xend=2,y=0.76,yend=0.76,size=2,colour=effect_colors[eff_interp_1])
p = p + annotate("rect",xmin=c(1.4), xmax=c(1.6), ymin=c(0.74) , ymax=c(0.78), color=effect_colors[eff_interp_1], fill=effect_colors[eff_interp_1])
p = p + annotate("text",x=1.5,y=0.8,size=6,label=paste('P : ',format(pval_actual,digits=3,scientific=T),sep=""))

p = p + annotate("segment",x=1,xend=3,y=0.85,yend=0.85,size=2)
p = p + annotate("rect",xmin=c(1.9), xmax=c(2.1), ymin=c(0.83) , ymax=c(0.87), color=effect_colors[eff_interp_2], fill=effect_colors[eff_interp_2])
p = p + annotate("text",x=2,y=0.89,size=6,label=paste('P : ',format(pval_actual_2,digits=3,scientific=T),sep=""))

# p = p + annotate("segment",x=2,xend=3,y=0.78,yend=0.78,size=2,colour=effect_colors[eff_interp_3])
# p = p + annotate("rect",xmin=c(2.4), xmax=c(2.6), ymin=c(0.76) , ymax=c(0.80), color=effect_colors[eff_interp_3], fill=effect_colors[eff_interp_3])
# p = p + annotate("text",x=2.5,y=0.82,size=6,label=paste('P : ',format(diff_test_3$p.value,digits=3,scientific=T),sep=""))

print(p)
dev.off()

file.copy(paste(dir_network,'Features.png',sep=""),paste(dir_s2,'fig_s2c.png',sep=""))

# -------------------------------------------
# Sec 08: Figure S3
# -------------------------------------------

dir_s3 = paste(dir_res,'fig_s3/',sep="")
dir.create(dir_s3)

feats_sel = read.csv(paste(dir_functional,'surv/final/feats.csv',sep=""),row.names = 1)

pdf(paste(dir_s3,"fig_s3c.pdf",sep=""))
df = data.frame(cbind(feats_sel,subjects[,'groupID']))
colnames(df) = c('f1','f2','label')
df$label = as.factor(df$label)

p = ggplot(df,aes(x=f2,y=f1,colour=label)) + geom_point(size=4,aes(colour=label)) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=20,colour='black'), axis.title.y = element_text(size=20,colour='black') ,axis.text.x = element_text(size=15,colour='black'), axis.text.y = element_text(size=15,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_color_manual(values=group_colors) + scale_fill_manual(values=group_colors) + xlab(colnames(feats_sel)[2]) + ylab(colnames(feats_sel)[1])
print(p)
dev.off()

file.copy(paste(dir_functional,'robust_test.pdf',sep=""),paste(dir_s3,'fig_s3d.pdf',sep=""))

for(featIdx in 1:ncol(feats_sel)){
  
  df = data.frame(cbind(feats_sel[,featIdx],subjects[,'groupID'],subjects[,'Challenges']))
  colnames(df) = c('feat','groupID','Challenges')
  df$groupID = as.factor(df$groupID)
  df$Challenges = as.factor(df$Challenges)
  
  featName = colnames(feats_sel)[featIdx]
  
  spea_all = cor.test(df[,1],subjects[,'Challenges'],method="spearman")
  spea_3 = cor.test(df[grp_3,1],subjects[grp_3,'Challenges'],method="spearman")
  spea_4 = cor.test(df[grp_4,1],subjects[grp_4,'Challenges'],method="spearman")
  spea_5 = cor.test(df[grp_5,1],subjects[grp_5,'Challenges'],method="spearman")
  
  spea_func = paste(round(spea_all$estimate,2),', p:',format(spea_all$p.value,scientific=T,digits=3),sep=" ")
  spea_func_3 = paste(round(spea_3$estimate,2),'p:',format(spea_3$p.value,scientific=T,digits=3),sep=" ")
  spea_func_4 = paste(round(spea_4$estimate,2),'p:',format(spea_4$p.value,scientific=T,digits=3),sep=" ")
  spea_func_5 = paste(round(spea_5$estimate,2),'p:',format(spea_5$p.value,scientific=T,digits=3),sep=" ")
  
  pdf(paste(dir_s3,"fig_s3efg_",featIdx,".pdf",sep=""))
  p1 = ggplot(df,aes(x=Challenges,y=feat)) + geom_point(aes(size=groupID,color=groupID,alpha=groupID)) + ylab(featName) + xlab('Challenges') + scale_x_discrete(labels=c(as.character(1:12),'UI')) + theme(axis.line = element_line(colour = "black",size=1), axis.title.x = element_text(size=25,colour='black'), axis.title.y = element_text(size=15,colour='black') ,axis.text.x = element_text(size=17,colour='black'), axis.text.y = element_text(size=17,colour='black'),panel.background = element_blank(), panel.grid.major.x = element_line(colour='gray65',size=0.3,linetype = 'dashed'), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none') + scale_colour_manual(name='Groupwise Coeff\nSpearman',values=group_colors,labels=c(spea_func_3,spea_func_4,spea_func_5)) + scale_size_manual(name='Overall\nSpearman',values=c(3.5,3.5,3.5),labels=c(spea_func,NA,NA)) + scale_alpha_manual(values=c(0.65,0.65,0.65))
  print(p1)
  dev.off()
  
}


cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
