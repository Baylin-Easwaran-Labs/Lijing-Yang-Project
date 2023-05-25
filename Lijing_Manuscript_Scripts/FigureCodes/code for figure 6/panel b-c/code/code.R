#设置工作路径
wd='C:/Users/86137/OneDrive/Desktop/R_test'
setwd(wd)
getwd()
list.files()
#setwd("C:/Users/86137/OneDrive/Desktop/R data")

#数据下载
rm(list = ls())  # make the environment empty. 
options(stringsAsFactors = F) #避免系统自动将数据框内的相关变量转换成因子变量。
library(GEOquery)
gse = "GSE39582"
eSet <- getGEO(gse, 
               destdir = '.', 
               getGPL = F)

list.files()

#(1)提取表达矩阵exp
exp <- exprs(eSet[[1]])
exp[1:4,1:4]
#exp就是一个表达矩阵



#exp = log2(exp+1)  
#我并不明白exp为什么要加1. 担心有0的数据，对数转换后会出现无穷小的数。
#不同来源的数据，基因表达的形式可能不一样，要注意鉴别啊。最好是学会自己分析原始数据。
#(2)提取临床信息
pd <- pData(eSet[[1]])  
colnames(pd)
write.table(pd,file='pd.csv',sep=',',row.names=FALSE)

pd=pd[,c('geo_accession','source_name_ch1','characteristics_ch1.26','dependancy sample:ch1','tumor.location:ch1')]
write.table(pd,file='pd_1.csv',sep=',',row.names=FALSE)




#感觉eSet就是一个复合文件。
# pd包含了样本临床先关的信息。

#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p   #逻辑函数，先查看是不是一致。
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl <- eSet[[1]]@annotation
gpl
save(gse,pd,exp,gpl,file = "step1output.Rdata")

load(file = "step1output.Rdata")

#先将CDX2的表达水平与临床相关信息整合。
#再筛选所需要的样本。
#去除正常组织；去除非独立样本；去除无肿瘤位置的样本。
#剩下的样本可以分为两组。近端 和 远端 结肠癌组。

rm(list = ls())   
#清除r里面的变量环境

load(file = "step1output.Rdata")


#先给基因名添加symbol.
#http://www.bio-info-trainee.com/1399.html
if(!require(hgu133plus2.db))BiocManager::install("hgu133plus2.db") #require()是逻辑函数
library(hgu133plus2.db)
ls("package:hgu133plus2.db")  #会有一系列的文件出来，应该是不同的命名对应关系。
ids <- toTable(hgu133plus2SYMBOL)
which(ids$symbol == 'CDX2')
ids[c(12736,34953),]

#找cdx2在基因表达矩阵中的数据
rownames(exp)
which(rownames(exp) == '206387_at')
which(rownames(exp) == '231606_at')

exp[c(15834,40861),]
x=exp[c(15834,40861),]
write.csv(x,file="exp_cdx2_1.csv",quote=F,row.names = T)
#对于表达矩阵有对应于同一基因的不同探针，取平均表达值高的那个探针。
#所以这里我们要选取的CDX2对应的探针是'206387_at'


cdx2=as.data.frame(exp[15834,])
cdx2
cdx2=data.frame(cdx2,rownames(cdx2))
names(cdx2)[2]='sample'
#现在把CDX2的表达与样本信息整合
nrow(pd)
names(pd)[1]='sample'

library(dplyr)
df_all=dplyr::inner_join(pd,cdx2,by ="sample")  #根据指定的共同列，合并表格。


data.frame(colnames(df_all))

colnames(df_all)[6]='cdx2'

#现在数据已经整合，下一步就是筛选输入数据了。
#需要剔除正常组织，非独立样本。

unique(df_all$source_name_ch1)
x=which(df_all$source_name_ch1 =='Frozen tissue of non tumoral colorectal mucosa')
df_all=df_all[-x,]
which(df_all$source_name_ch1 =='Frozen tissue of non tumoral colorectal mucosa')

##分析左右结肠cdx2的表达，理论上并不需要是否braf突变信息。
#但因为后续有与braf是否突变的卡方分析，所以我在这还是剔除了无BRAF突变信息的样本。
unique(df_all$characteristics_ch1.26)
x=which(df_all$characteristics_ch1.26 == "braf.mutation: N/A")
nrow(df_all)
length(x)
df_all=df_all[-x,]
nrow(df_all)


unique(df_all$`dependancy sample:ch1`)
length(unique(df_all$sample))


#筛除没有左 右结肠信息的样本
unique(df_all$`tumor.location:ch1`)
#都有

nrow(df_all)
#最后合计有512个样本纳入的分析。

#现在就可以统计数据出图了。x轴为左右结肠信息，y轴为cdx2基因表达信息


library(ggplot2)
data3=df_all
colnames(data3)[5]='LOCATION'
data3$LOCATION <- factor(data3$LOCATION, levels = c("proximal","distal"))

p <- ggplot(data3, aes(x=LOCATION, y=cdx2)) 
p
#加统计量进去。
p + geom_violin(aes(fill = LOCATION), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")


# 改指定颜色，想怎么调就怎么调。
p2 <- p + geom_violin(aes(fill = LOCATION), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))

p4

pdf('abc.pdf',10,6)
print(p4)
dev.off()

pdf('abc_1.pdf',10,10)
print(p4)
dev.off()


#做个T-test。
t.test(data3$cdx2~data3$LOCATION)



##现在开始分析cdx2表达与左右结肠 或 braf突变的相关性。
#当然是选用卡方检验啦。
df_all
df_all$cdx2_1='normal'



for(i in 1:nrow(df_all)){
  if(df_all$cdx2[i] <= 6.5){df_all$cdx2_1[i] = 'negative'}
  else {df_all$cdx2_1[i] = 'positive'}
}


############################
#Just focus on the correlation between Cdx2 and location
names(df_all)[7]='CDX2_1'
names(df_all)[5]='LOCATION'
names(df_all)[3]='BRAF'


library(ggstatsplot)
library(ggplot2)
library(dplyr)



data_2=df_all
head(data_2)
data_6=dplyr::select(data_2, c(CDX2_1,LOCATION))
data_6
table(data_6)
nrow(data_6)


data_6$CDX2_1=factor(data_6$CDX2_1, levels=c('positive','negative'))
data_6$CDX2_1
data_6$LOCATION=factor(data_6$LOCATION, levels=c('distal','proximal'))
data_6$LOCATION
ggbarstats(data_6,LOCATION,CDX2_1,palette = 'Set2')
abc=ggbarstats(data_6,LOCATION,CDX2_1,palette = 'Set2')
pdf('abc_4.pdf',10,10)
print(abc)
dev.off()


test=table(data_6)
test
chisq.test(test)
chisq.test(test)$expected





##########################
#Just focus on the correlation between Cdx2 and braf mutation
data_7=df_all
head(data_7)
data_7=dplyr::select(data_7, c(CDX2_1,BRAF))
data_7
table(data_7)
data_7$CDX2_1=factor(data_7$CDX2_1, levels=c('positive','negative'))
data_7$CDX2_1
data_7$BRAF=factor(data_7$BRAF, levels=c('braf.mutation: WT','braf.mutation: M'))
data_7$BRAF
ggbarstats(data_7,BRAF,CDX2_1,palette = 'Set2')
abc=ggbarstats(data_7,BRAF,CDX2_1,palette = 'Set2')

pdf('abc_5.pdf',10,10)
print(abc)
dev.off()

table(data_7)
test=table(data_7)
test
chisq.test(test)
chisq.test(test)$expected

#################

