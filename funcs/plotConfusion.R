# ******************************************************************************
# Author      : Srivamshi Pittala
# Advisor     : Prof. Chris Bailey-Kellogg
# Project     : NIH 10-332
# Description : Plot confusion matrix for multinomial classification
# Cite        : TBD
# ******************************************************************************

plotConfusion = function(conf_mat,labels,dir_class){
  
  actual = as.data.frame(table(labels))
  names(actual) = c("Reference","ReferenceFreq")
  
  df = data.frame(conf_mat$table)
  df = merge(df,actual,by=c('Reference'))
  df$Rate = df$Freq/df$ReferenceFreq
  
  pdf(paste(dir_class,'confusion.pdf',sep=''),family=plot_font)
  p = ggplot(df,aes(x=Reference, y=Prediction,fill=Rate),color='black',size=1) + geom_tile(color="black",size=0.5) + labs(x="\nReference",y="Prediction\n") + scale_x_discrete(labels=c('IM Mosaic','IM 239','AE 239')) + scale_y_discrete(labels=c('IM Mosaic','IM 239','AE 239')) + theme(axis.title.x = element_text(size=16,colour='black'), axis.title.y = element_text(size=18,colour='black'), axis.text.y=element_text(angle=45,size=16,colour = 'black'),axis.text.x=element_text(size=18,colour = 'black'),axis.ticks=element_blank(),panel.background = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(), legend.position='none')
  p = p + geom_text(aes(x=Reference, y=Prediction, label=sprintf("%.2f", Rate)),data=df, size=10, colour="black") + scale_fill_gradient(low="white",high="#458B00",limits=c(0,1))
  print(p)
  dev.off()
  
}