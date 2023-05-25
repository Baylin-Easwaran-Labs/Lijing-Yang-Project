
wd=c('C:/Users/86137/OneDrive/Desktop/nature_genetics/figures/codes for figure 4/input_data')
setwd(wd)
getwd()
list.files()


library(dplyr)
library(stringr)
temp=list.files(path=wd,pattern="*.txt")
#列出文件夹中的txt文件。也可以列出其他文件，例如csv文件。
temp
length(temp)

data=list()
for(i in 1:length(temp)){
  data[[i]]=read.table(temp[i],header=T,sep='',na.strings = c("NA"))
}
temp

for(i in 1:length(temp)){
  dataset=data[[i]]
  cut_off_pvalue=0.05
  cut_off_logFC=1
  dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, ifelse(dataset$log2FoldChange >= cut_off_logFC,'Up','Down'),'Stable' )
  #abs() 是取绝对值的意思。
  head(dataset)
  
  library(ggplot2)
  p = ggplot(dataset,aes(x=log2FoldChange,y=-log10(padj),colour=change))+
    geom_point(alpha=0.4,size=3.5)+
    scale_color_manual(values=c("blue","grey", "red"))+
    geom_vline(xintercept=c(-1,1),lty=4,col='black',lwd=0.8)+
    geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col='black',lwd=0.8)+
    labs(title=str_c(substr(temp[i],14,nchar(temp[i])-4), "_volcano"),x='log2FoldChange',y='-log10(padj)')+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position='right',
          legend.title = element_blank())
  print(i)
  pdf(str_c(as.character(i),'_',substr(temp[i],14,nchar(temp[i])-4), "_volcano.pdf"),8,6)
  print(p)
  dev.off()
  
}

