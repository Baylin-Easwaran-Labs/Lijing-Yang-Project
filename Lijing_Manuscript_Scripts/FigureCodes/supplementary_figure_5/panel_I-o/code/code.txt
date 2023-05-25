wd='C:/Users/86137/OneDrive/Desktop/nature_genetics/supplemental figures/sup_figure_5/panel_I-o'
setwd(wd)
getwd()

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(org.Mm.eg.db)
#管他有没有用，先library了再说。
list.files()


x <- readLines("Parikh_CellCluster_GeneSetList.gmt")
res <- strsplit(x, "\t")
#将x按"\t"分割符分割，生成一列表。
names(res) <- vapply(res, function(y) {y[1]}, character(1))
#function(y) {y[1]}自定义的功能函数，输入y，输出y[1]
#character(1) 字符""，表示一个字符。
#y[1]的意义有点难理解，就是列表中每个亚列表的第一个元素。
res <- lapply(res, "[", -c(1:2))
#"["这个具体意义我不是很懂。可能是针对元素删减的操作符。删掉第1,2行。
length(res)



gene_set=res

gene_set[1]
gene_set[2]
#############################
#BiocManager::install('homologene')
library(homologene)


for(i in 1:length(res)){
  genelist<-gene_set[[i]]
  homologene(genelist, inTax = 9606, outTax = 10090)
  #人转小鼠，函数里面的数值是物种的代码。9606表示人，10090表示小鼠。
  a=homologene(genelist, inTax = 9606, outTax = 10090)
  a
  nrow(a)
  length(unique(a$`10090`))
  a$label_1=rep(names(res[i]))
  a$label_2=rep(names(res[i]))
  gene_set[[i]]=dplyr::select(a,c(4,label_1,label_2))
  #可能因为10090_ID以数值开头，有些函数不识别。所以换为第4行。
  gene_set[[i]][,1]=as.character(gene_set[[i]][,1])
  #患者字符类型才能行后续的对比。
  gene_set[[i]][,1]
}

gene_set[[i]][,1]
i
#正确运行的话，此时i应该等于50.
#上述代码主要预先处理geneset。将人的转换成老鼠的。



gene_set_1=data.frame()
for(i in 1:length(res)){
  gene_set_1=rbind(gene_set_1,gene_set[[i]])
}
#cbind()不能运行。rbind()能运行。gene_set_1为空数据时，不能运行cbind()。但能运行rbind()。奇怪。
gene_set_1
go2gene  <- gene_set_1[, c(2, 1)]
go2name  <- gene_set_1[, c(2, 3)]
#主要的得到go2gene,和go2name。后边的分析需要这两个变量。下一步是输入数据的预处理。
#geneset文件已经整理成代码要求的格式。下一步就是要准备gene表达数据了。



###################
temp=list.files(path=wd,pattern="*.txt")
#列出文件夹中的txt文件。也可以列出其他文件，例如csv文件。
temp
length(temp)
temp

data_1=list()
for(i in 1:length(temp)){
  data_1[[i]]=read.table(temp[i],header=T,sep='',na.strings = c("NA"))
}
#现在输入数据储存在data_1这个数据框中了。
data_1


#i=1
for(i in 1:8){
  DEG_mtx=data_1[[i]]
  df=dplyr::select(DEG_mtx,c(SYMBOL,log2FoldChange))
  df_id=bitr(df$SYMBOL,
             fromType = 'SYMBOL',
             toType = "ENTREZID",
             OrgDb = 'org.Mm.eg.db'
  )
  #小鼠的是org.Mm.eg.db，人的是org.Hs.eg.db
  df_all=dplyr::inner_join(df,df_id,by ="SYMBOL")
  
  
  df_all_sort=df_all[order(df_all$log2FoldChange,decreasing = T),]
  head(df_all_sort)
  
  gene_fc=df_all_sort$log2FoldChange
  head(gene_fc)
  names(gene_fc)=df_all_sort$ENTREZID
  head(gene_fc)
  names(gene_fc)
  
  x <- GSEA(
    gene_fc,
    TERM2GENE = go2gene,
    TERM2NAME = go2name,
    pvalueCutoff=1)
  #pvalueCutoff=1 可以输出所有结果，但并不是所有结果有统计学意义。
  
  #准备一个给保存文件取名的变量。
  library(stringr)
  file_name=temp[i]
  test=file_name
  test=substr(test,1,nchar(test)-4)
  test=substr(test,14,nchar(test))
  test
  file_name=test
  
  #paths=names(res)
  #paths=c('HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT')
  #paths=names(x@geneSets)
  paths=c("Colonocytes","Enteroendocrines","Goblets","Undifferentiated #1","Crypt Top Colonocytes")
  figure=gseaplot2(x,paths,color='firebrick',rel_heights=c(1,.1,.15),title=file_name)
  figure_1=gseaplot2(x,paths,color='firebrick',pvalue_table=TRUE,rel_heights=c(1,.1,.15),title=file_name)
  print(file_name)
  
  #保存分析的数据，包括NES和pvalue
  sortKEGG = x[order(x$enrichmentScore,decreasing =T),]
  head(sortKEGG)
  dim(sortKEGG)
  write.csv(sortKEGG,file=str_c(file_name,'_gsea','.csv'),quote=F,row.names = T)
  #把统计数据保存起来
  
  #将分析的图以PDF文件形式保存起来。
  pdf(str_c(file_name, ".pdf"),10,10)
  print(figure)
  print(figure_1)
  dev.off() 
  print('finish')
  print(i)
  
}



