wd='C:/Users/86137/OneDrive/Desktop/nature_genetics/supplemental figures/sup_figure_5/panel_f'
setwd(wd)

file=list.files(wd, pattern='.txt')
length(file)
file

i=1
data=list()

for(i in 1:length(file)){
  data[[i]]=read.table(file[i],header=T,sep='',na.strings = c("NA"))
  
}
length(data)



names(data)=file
names(data)


table=data.frame(file)
table
i=2

for(i in 1:length(file)){
a=data[[i]]
names(data[i])
names(a)
a=a[,c("ENTREZID", "SYMBOL", "log2FoldChange", "padj" )]
head(a)

nrow(a)
up=length(which(a$log2FoldChange >= 1 & a$padj <0.05))
down=length(which(a$log2FoldChange <= (-1) & a$padj <0.05))
up
down

table[i,2]=up
table[i,3]=down

}
table

names(table)=c('group','up','down')

write.csv(table,file="table.csv",quote=F,row.names = F)



###############################
library(ggplot2)
setwd(wd)
data <- read.csv(file="table.csv", header=T)
data$down

y=c(data$up,data$down)
x=c(data$group,data$group)
fill=c(rep('up',4),rep('down',4))
data=data.frame(x,y,fill)
data$x
nchar('annotatedRes_')
data$x=substr(data$x,14,nchar(data$x)-7)
data

for(i in 1:nrow(data)){
  if(data$fill[i]=='down'){data$y[i]= (-1)*data$y[i]}
}

unique(data$x)
data$x <- factor(data$x, levels = c("FMVsBM_Prox","Prox.FMVsBM_Cdx2KO","FMVsBM_Dist",'Dist.FMVsBM_Cdx2KO'))


figure=ggplot(data, aes(x = x, y = y, fill = fill)) +
  geom_bar(stat = "identity",width = 0.9, position = position_dodge(0.9))+
  geom_text(aes(x=x,y=y+100*y/abs(y),label=abs(y)),position=position_dodge(width=0.9),color='black')+
  ggtitle('number_geneExpression')+
  scale_fill_discrete(name='change')+
  #scale_fill_manual(name='change',values=c('blue','red'))+
  #可以修改柱子的颜色。不能和scale_fill_discrete 同时使用。
  #scale_y_continuous(name='???',breaks=seq(-120,40,20),limits=c(-120,40)) #y??????????????????,??????????????????
  theme_bw() 
#theme(panel.grid.major=element_line(colour=NA)) #????????????
figure

figure_1=ggplot(data, aes(x = x, y = y, fill = fill)) +
  geom_bar(stat = "identity",width = 0.9, position = position_dodge(0.9))+
  geom_text(aes(x=x,y=y+100*y/abs(y),label=abs(y)),position=position_dodge(width=0.9),color='black')+
  ggtitle('number_geneExpression')+
  #scale_fill_discrete(name='change')+
  scale_fill_manual(name='change',values=c('blue','red'))+
  #可以修改柱子的颜色。不能和scale_fill_discrete 同时使用。
  #scale_y_continuous(name='???',breaks=seq(-120,40,20),limits=c(-120,40)) #y??????????????????,??????????????????
  theme_bw() 



pdf('number_geneExpression.pdf',10,6)
print(figure)
print(figure_1)
dev.off() 

