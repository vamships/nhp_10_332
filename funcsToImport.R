# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Imports the necessary packages, functions, and defines global variables
# Cite        : TBD
# ******************************************************************************

library(glmnet)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(survival)
library(corrplot)
library(caret)
library(survcomp)
library(polycor)
library(elasticnet)
library(effsize)
library(e1071)

# -----------------------------------------------------
source('funcs/clusterFeats.R')
source('funcs/convertKMtoEvents.R')
source('funcs/coxSurvivalWithBackSearch.R')
source('funcs/coxSurvivalFinalModel.R')
source('funcs/coxSurvivalLookBack.R')
source('funcs/createColumnColors.R')
source('funcs/createSubjectColors.R')
source('funcs/doCoxBackSearch.R')
source('funcs/doFullSurvivalAnalysis.R')
source('funcs/extractProbabilityFromKM.R')
source('funcs/glmnetMultiClass.R')
source('funcs/heatmap4.R')
source('funcs/plotConfusion.R')
source('funcs/selectFromClusters.R')
source('funcs/takeOffOneFeat.R')
# -----------------------------------------------------
# Define colors

# group type

group_id = c(3,4,5)
names(group_id) = c('IM Mosaic','IM 239','AE 239')

group_colors = c("#FF8000","#008000","#400080")
#names(group_colors) = c('Control','Gag','IM Mosaic','IM 239','AE 239')
names(group_colors) = group_id

reagent_names = c('rR2A','R2A','R2A.2','R2A.3','R2A.4.low','R2A.4.high','rR3A.1','R3A.1','R3A.3','C1q','MBL','aRhIgG.PE.low','aRhIgG.PE.high','hFcgRIIIa','hFcgRIIIb','hFcgRIIa','rIgG','IgA')

reagent_colors = c(rR2A='#690B65',R2A='#690B65',R2A.2='#690B65',R2A.3='#690B65',R2A.4.low='#690B65',R2A.4.high='#690B65',rR3A.1="#10E88C",R3A.1="#10E88C",R3A.3="#10E88C",C1q="#DC28D6",MBL="plum2",aRhIgG.PE.low="orange2",aRhIgG.PE.high="orange2",hFcgRIIIa="darkseagreen",hFcgRIIIb="olivedrab",hFcgRIIa='mediumblue',rIgG="orange2",IgA='#2F587B')

antigen_names = c('C1.Ak','C1.TR','G73','G119','gp120','gp130','gp140','Gag','V1a','G146','G49','V1V2','2A5','84','His')

antigen_colors = c(C1.Ak='antiquewhite',C1.TR='antiquewhite',G73='antiquewhite',G119='antiquewhite',gp120='burlywood4',gp130='lightseagreen',gp140='khaki3',Gag='#F4A895',V1a='mediumorchid',G146='mediumorchid',G49='mediumorchid',V1V2='mediumblue','2A5'='mediumblue','84'='mediumblue',His='khaki3')

antigen_colors_broad = c('C peptides'='antiquewhite',gp120='lightskyblue',gp130='lightseagreen',gp140='khaki2',Gag='#F4A895','V peptides'='mediumorchid','Probes'='mediumblue')

ak_colors = c(Both="gray50",None="gray95",Resistant="black",Sensitive="gray75")

# colors based on how many challenges they survived
challenge_colors = colorRampPalette(c('beige','bisque4'))(13);
challenge_colors = colorRampPalette(c('honeydew2','grey34'))(13);
challenge_colors[13] = '#5DFC0A'
names(challenge_colors) = c(1:13)

# colors based on protection
interest_colors = c('purple','orange')

effect_colors = c('grey83','grey63','grey43','grey23')
names(effect_colors) = c('negligible','small','medium','large')

func_names = c('ADNP','ADCP','ADCC','ADCD','CD107','IFNg','MIP1b')
func_colors = c(ADNP='#F57D3A',ADCP='#2D7BBB',ADCC='#B49495',ADCD='#16412B',CD107='#95C665',IFNg='#690B65',MIP1b='#D2201A')

tp_names = c('tp5','tp6')
tp_colors = c(tp5='#FE746A',tp6='#FC8AD2')
