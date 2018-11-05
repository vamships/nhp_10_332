# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Performs multinomial logistic classification using elastic net
# Cite        : TBD
# ******************************************************************************

glmnetMultiClass = function(feats,predVec,weights,numFeat,grpType,intc,alpha,cvFolds,repeatRun){
  
  augMatrix = cbind(feats,predVec)
  
  rmList = which(is.na(predVec))
  
  if(length(rmList)!=0){
    
    if(cvFolds==nrow(augMatrix))
      cvFolds = cvFolds - length(rmList)
    
    augMatrix = augMatrix[-rmList,]
  }
  
  cat("removed : ",length(rmList)," subjects.\n")
  cat("num folds : ",cvFolds,"\n")
  
  feats = augMatrix[,1:numFeat]
  func = augMatrix[,(numFeat+1)]
  
  numClasses = length(unique(func))
  
  # Cross-validation
  
  cv_repeat = matrix(numeric(1),repeatRun,2)
  cv_permut = matrix(numeric(1),repeatRun,2)
  
  colnames(cv_repeat) = c('min','se1')
  colnames(cv_permut) = c('min','se1')
  
  feat_betas = matrix(0,nrow=numFeat,ncol=numClasses)
  feat_counts = matrix(0,nrow=numFeat,ncol=repeatRun)
  
  for(testIdx in seq(1,repeatRun)){
    
    cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,family="multinomial",standardize=FALSE,weights=weights,alpha=alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
    
    cv_perm = cv.glmnet(feats,func[sample(nrow(feats))],nfolds=cvFolds,family="multinomial",standardize=FALSE,weights=weights,alpha=alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
    
    cv_repeat[testIdx,'min'] = cv_glm$cvm[match(cv_glm$lambda.min,cv_glm$lambda)]
    cv_repeat[testIdx,'se1'] = cv_glm$cvm[match(cv_glm$lambda.1se,cv_glm$lambda)]
    
    cv_permut[testIdx,'min'] = cv_perm$cvm[match(cv_perm$lambda.min,cv_perm$lambda)]
    cv_permut[testIdx,'se1'] = cv_perm$cvm[match(cv_perm$lambda.1se,cv_perm$lambda)]
    
    cv_glm_coef_min = coef(cv_glm,s="lambda.min")
    betas_min = matrix(0,nrow=numFeat,ncol=numClasses)
    
    for(classIdx in 1:numClasses){
    betas_min[,classIdx] = cv_glm_coef_min[[classIdx]][-1]
    }
    
    feat_betas = feat_betas + betas_min
    
    feat_counts[,testIdx] = ((1:numFeat) %in% which(rowSums(abs(betas_min))!=0))*1
    
  }
  
  # final model
  set.seed(3454)
  folds_list = createFolds(func,cvFolds,list=F)
  cv_glm = cv.glmnet(feats,func,nfolds=cvFolds,foldid=folds_list,family="multinomial",standardize=FALSE,weights=weights,alpha=alpha, type.multinomial = grpType,type.measure="class",keep=TRUE,intercept=intc)
  rm(.Random.seed,envir = globalenv())
  
  #save the prevalidated array for viz
  fit_preval_min = cv_glm$fit.preval[,,match(cv_glm$lambda.min,cv_glm$lambda)]
  fit_preval_1se = cv_glm$fit.preval[,,match(cv_glm$lambda.1se,cv_glm$lambda)]
  
  #get feature coefficients for the model
  cv_glm_coef_min = coef(cv_glm,s="lambda.min")
  cv_glm_coef_1se = coef(cv_glm,s="lambda.1se")
  
  mse_min = cv_glm$cvm[match(cv_glm$lambda.min,cv_glm$lambda)]
  mse_1se = cv_glm$cvm[match(cv_glm$lambda.1se,cv_glm$lambda)]
  
  betas_min = matrix(0,nrow=numFeat,ncol=numClasses)
  rownames(betas_min) = colnames(feats)
  betas_1se = matrix(0,nrow=numFeat,ncol=numClasses)
  rownames(betas_1se) = colnames(feats)
  
  for(classIdx in 1:numClasses){
    
    betas_min[,classIdx] = cv_glm_coef_min[[classIdx]][-1]
    betas_1se[,classIdx] = cv_glm_coef_1se[[classIdx]][-1]
    
  }
  
  min_model = list("lambda"=cv_glm$lambda.min,"mse"=mse_min,"preval"=fit_preval_min,"coeff"=betas_min)
  se1_model = list("lambda"=cv_glm$lambda.1se,"mse"=mse_1se,"preval"=fit_preval_1se,"coeff"=betas_1se)
  
  best_model = min_model
  
  return(list("removed"=rmList,"final_func"=func,"best_alpha"=alpha,"final_fit"=cv_glm,"best_model"=best_model,"min_model"=min_model,"se1_model"=se1_model,"repeat_mse_min"=mean(cv_repeat[,'min']),"feat_betas"=feat_betas,"feat_counts"=feat_counts,"repeat_mse_min_sd"=sd(cv_repeat[,'min']),"repeat_mse_se1"=mean(cv_repeat[,'se1']),"repeat_mse_se1_sd"=sd(cv_repeat[,'se1']),"cv_repeat"=cv_repeat,"permut_mse_min"=mean(cv_permut[,'min']),"permut_mse_min_sd"=sd(cv_permut[,'min']),"permut_mse_se1"=mean(cv_permut[,'se1']),"permut_mse_se1_sd"=sd(cv_permut[,'se1']),"cv_permut"=cv_permut))
  
}