# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Performs survival analysis using the biophysical measurements
#               and permutation test with size-matched random selection
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

dir_res = paste('results_coxPred_randomSel/',sep="")
dir.create(dir_res)

# -------------------------------------------
# Sec 01: Hyper-parameters
# -------------------------------------------

num_folds = 10 # number of folds for cross-validation
num_repeat = 100 # number of repetitions of cross-validation
seed_final = 8357 # seed to determine folds for final cross-validation model
randSelSize = 4 # number of features in the final model of biophysical

log_file = paste(dir_res,'log_file',sep="")
file.create(log_file)
cat(rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)

cat('\nCox-model parameters','\n',sep="",file=log_file,append=T)
cat('Number of folds : ',num_folds,'\n',file=log_file,append=T)
cat('Number of repeated evaluations : ',num_repeat,'\n',file=log_file,append=T)
cat('Seed for final model : ',seed_final,'\n',file=log_file,append=T)
cat('Number of features to randomly select : ',randSelSize,'\n',file=log_file,append=T)

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
# Sec 04: Robustness test
# -------------------------------------------

dir_perm = paste(dir_res,"perm/",sep="")
dir.create(dir_perm)

perm_final = matrix(NA,nrow=num_repeat,ncol=1)
colnames(perm_final) = 'cindex_test'

perm_log = paste(dir_perm,'perms.txt',sep="")
file.create(perm_log)

for(testIdx in 1:num_repeat){
  
  cat(rep("#",30),"\n",sep="")
  cat("Permutation test :",testIdx,"\n")
  
  feats_sample = luminex_tp5
  sample_size = ncol(feats_sample)
  
  perm_order = sample(sample_size,randSelSize)
  feats_perm = data.frame(feats_sample[,perm_order])
  
  cat(testIdx,':',perm_order,'\n',file=perm_log,append=T)
  
  dir_test = paste(dir_perm,"test_",testIdx,"/",sep="")
  dir.create(dir_test)
  
  model_perm = coxSurvivalFinalModel(feats_perm,subjects,seed_final,scolors,group_colors,num_folds,plots=TRUE,dir_test)
  
  perm_final[testIdx,'cindex_test'] = model_perm$cindex_test
  
}

# -------------------------------------------
# Sec 05: Comparing actual and permuted models' performance
# -------------------------------------------

# -------------
# Comparing repeated cross-validation between actual and permuted
# -------------

df_test = as.data.frame(cbind(rep(1,num_repeat),perm_final[,'cindex_test']))
colnames(df_test) = c('label','C_index')

write.csv(df_test,file=paste(dir_perm,'robust.csv',sep=""),row.names = T)

cat('\n\n',rep('#',8),format(Sys.time(),"%d/%m/%Y %X"),rep('#',8),'\n\n',file=log_file,append=T)
