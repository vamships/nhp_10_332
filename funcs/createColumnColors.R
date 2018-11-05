# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Decription  : Create column color panel based on the reagents and antigens
# Cite        : TBD
# ******************************************************************************

createColumnColors = function(featNames,reagent_names,reagent_colors,antigen_names,antigen_colors){
  
  lcolors = matrix(nrow=length(featNames),ncol=2,dimnames=list(featNames,c('reagent','antigen')))
  
  for(reagent_name in reagent_names){
    
    lcolors[grep(reagent_name,featNames),'reagent'] = reagent_colors[reagent_name]
    
  }
  
  for(antigen_name in antigen_names){
    
    lcolors[grep(antigen_name,featNames),'antigen'] = antigen_colors[antigen_name]
    
  }
  
  return(lcolors)
  
}

createFuncColors = function(featNames,func_names,func_colors,tp_names,tp_colors){
  
  lcolors = matrix(nrow=length(featNames),ncol=2,dimnames=list(featNames,c('func','timepoint')))
  
  for(func_name in func_names){
    
    lcolors[grep(func_name,featNames),'func'] = func_colors[func_name]
    
  }
  
  for(tp_name in tp_names){
    
    lcolors[grep(tp_name,featNames),'timepoint'] = tp_colors[tp_name]
    
  }
  
  return(lcolors)
  
}