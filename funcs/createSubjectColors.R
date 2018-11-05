# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Create row colors based on sample's group ID and survival time
# Cite        : TBD
# ******************************************************************************

createSubjectColors = function(subjects,group_colors,challenge_colors){
  
  # create row colors
  scolors = matrix(nrow=nrow(subjects), ncol=2, dimnames=list(rownames(subjects), c('groupID','challenges')))
  
  for (i in 1:nrow(scolors)) {
    scolors[i,'groupID'] = group_colors[toString(subjects[i,'groupID'])]
    scolors[i,'challenges'] = challenge_colors[subjects[i,'Challenges']]
    if(subjects[i,'censor']==0){
      scolors[i,'challenges'] = '#5DFC0A'
    }
  }
  
  return(scolors)
  
}