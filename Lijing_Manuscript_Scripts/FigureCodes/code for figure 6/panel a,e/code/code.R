#设置工作路径
wd='C:/Users/86137/OneDrive/Desktop/R_test'
setwd(wd)
getwd()
list.files()

#download the files of TCGA data that I used in the following analysis
url <- 'https://data.cyverse.org/dav-anon/iplant/home/lijingyang23/TCGA_raw_data/rawdata.zip'

path='C:/Users/86137/OneDrive/Desktop/R_test'

download.file(url,'rawdata.zip')
list.files()
unzip('rawdata.zip')
list.files()



#在工作路径下建立一个文件夹
dir.create('data_in_one')
#显示rawdata/ 路径下的文件。
dir("rawdata/")

#利用代码，将*FPKM.txt文件拷贝到data_in_one中。
for (dirname in dir("rawdata/")){
  file <- list.files(paste0(getwd(),"/rawdata/",dirname),pattern = "*FPKM.txt")
  file.copy(paste0(getwd(),"/rawdata/",dirname,"/",file),"data_in_one") 
}
#这个循环有意思。理解后其实蛮简单的。

#网站上下载的与测序文件相关的信息文件，其上有病人的索引信息。
#但这里面并没有病人的临床信息。
metadata <- jsonlite::fromJSON('metadata.cart.2021-12-19.json')
metadata$associated_entities
head(metadata$associated_entities)


naid_df <- data.frame()
for (i in 1:nrow(metadata)){
  naid_df[i,1] <- metadata$file_name[i]
  naid_df[i,2] <- metadata$associated_entities[i][[1]]$entity_submitter_id
  naid_df[i,3] <- metadata$associated_entities[i][[1]]$case_id
}
colnames(naid_df) <- c("filename","TCGA_id",'patient_uuid')
length(naid_df$patient_uuid)
length(unique(naid_df$patient_uuid))
#通过比较我们能发现，patient_uuid 是有重复名字的。
#可以理解为 1个病人可能会取多个标本。

length(naid_df$TCGA_id)
length(unique(naid_df$TCGA_id))
##通过比较我们能发现，TCGA_id 是没有重复名字的。 
#这是标本信息，所以标本标识是唯一的。但不同的标本可以来源于同一个病人。

head(naid_df)
write.table(naid_df,file='naid_df.csv',sep=',',row.names=FALSE)
#metadata是一个列表形式。注意metadata[[3]]里面的数据类型，其里面还套了一个表格。
#metadata$associated_entities 这里面是数据套数据，连环套连环。理解后就简单了。

test <- data.table::fread(paste0("data_in_one/",naid_df$filename[1]))
head(test)
#naid_df$filename[1]是读取的文件名。"data_in_one/" 是文件所在的路径。说白了，
#就是用函数读取文件。

expr_df <- data.frame(matrix(NA,nrow(test),nrow(naid_df)))
#我觉得这是一种很好的构建数据框的办法。
head(expr_df)
#做一个数据框，以基因的数据为行，样本的数据为列。


#下面这个循环时间漫长，慎重运行。
######################################
for (i in 1:nrow(naid_df)) {
  print(i)
  expr_df[,i]= data.table::fread(paste0("data_in_one/",naid_df$filename[i]))[,2]
}
head(expr_df)
#把每个样本的基因表达数据就整到一起了。这个运行时间有点长啊。
######################################



colnames(expr_df) <- naid_df$TCGA_id
head(expr_df)
#TCGA_id是每个样本的信息。是唯一的。

gene_id <- data.table::fread(paste0("data_in_one/",naid_df$filename[1]))$V1
expr_df <- cbind(gene_id,expr_df)
head(expr_df)


#tail(expr_df$gene_id,10)
#expr_df <- expr_df[1:(nrow(expr_df)-5),]
#这行命令是针对count数据类型而言的。其数据后面5行没什么用，需要删除掉。
#但针对FPKM数据，后面5行是有用的。


expr_df
save(expr_df,file = "expr_df.Rdata")
load(file = "expr_df.Rdata")
#write.table(expr_df,file='expr_df.csv',sep=',',row.names=FALSE)
#write.csv(expr_df,file="expr_df_2.csv",quote=F,row.names = F)
# 每个样本的rna表达量我们清楚了。下一步就是要设置组做分析了啊。
#所以首先要准备分组信息了。


#可以从这读入expr_df的数据
#load(file = "expr_df.Rdata")
expr_df_copy=expr_df


#去除重复的样本信息。但问题是，这些样本信息需要去除么？先去除试下吧。
a=colnames(expr_df)[-1]
length(a)
length(unique(a))
#可以判断下是否有重名，如果有的话可以跑下面代码去重。
#经分析没有重复
#因为列名是每个样本的信息。所以没有重复很正常。
#如果是患者信息的话，就会有重复信息了。
head(a)
b=0
c=data.frame()
for(i in 1:(length(a)-1)) {
  for(j in (i+1):length(a)){
    if(a[i] == a[j]){b = b+1
    c[b,1:2]=c(i,j)
    }
  }
}
b
c
#gcta_id没用重复数据。



######################################################
#expr_df中的01表示肿瘤，11表示癌旁正常组织。
expr_df=expr_df_copy
length(names(expr_df))
#注意：第一列是基因的名字。
#所以实际上是514个样本。


test=names(expr_df)
test=substr(test,start=14,stop=15)

test
unique(test)
# 11 表示癌旁正常组织。
# 01Primary Soild Tumor
# 02Recurrent Soild Tumor
# 06Metastatic

expr_df_01=expr_df[,-which(test =='11')]
length(names(expr_df_01))
expr_tumour=expr_df_01
#表示只取——肿瘤组织--的基因表达信息。


#正常组织的基因表达信息
expr_df=expr_df_copy
length(names(expr_df))
test=names(expr_df)
test=substr(test,start=14,stop=15)
expr_df_02=expr_df[,c(1,which(test =='11'))]
length(names(expr_df_02))
expr_normal=expr_df_02
#表示取癌旁（也就是正常组织）标本的表达信息。






#为了方便，后面的代码保持一致。但要注意，输入的基因表达信息是癌旁正常组织
#还是肿瘤组织呢。
#读入患者编码及肿瘤位置信息
data <- read.csv(file="patient_tumorLocation.csv", header=T)
nrow(data)
#这里是病人的信息。有459个病人信息。
#因为有的病人取取了多个标本，所以样本的数目是514个。
data
head(data)
unique(data$anatomic_neoplasm_subdivision)
#去掉data中肿瘤位置信息包含NA的数据

data$bcr_patient_uuid=tolower(data$bcr_patient_uuid)
data$anatomic_neoplasm_subdivision=tolower(data$anatomic_neoplasm_subdivision)
head(data)
unique(data$anatomic_neoplasm_subdivision)
data=data[which(data$anatomic_neoplasm_subdivision != '[discrepancy]'),]
data=data[which(data$anatomic_neoplasm_subdivision != '[not available]'),]
unique(data$anatomic_neoplasm_subdivision)
nrow(data)
data
nrow(data)
#对data里面的数据先清洗下，因为我们只搜集左右结肠的位置信息。

naid_df$patient_uuid=tolower(naid_df$patient_uuid)

# 就是用两个数据中的共同信息来整合这两个数据。
a=naid_df
b=data
for(i in 1:nrow(a)){
  for (j in 1:nrow(b)){
    if(a$patient_uuid[i] == b$bcr_patient_uuid[j]){a[i,4]=b$anatomic_neoplasm_subdivision[j]}
  }
}
#就是为了给TCGA_id后面加上肿瘤/正常组织的取材位置信息
nrow(a)
unique(a$V4)
which(is.na(a$V4))
a=a[-which(is.na(a$V4)),]
nrow(a)
unique(a$V4)

n1=which(a$V4=="ascending colon")
n2=which(a$V4=="sigmoid colon")
n3=which(a$V4=="cecum")
n4=which(a$V4=="hepatic flexure")
n5=which(a$V4=="splenic flexure")
n6=which(a$V4=="descending colon")
a=a[c(n1,n2,n3,n4,n5,n6),]
nrow(a)

n1=which(a$V4=="ascending colon")
n2=which(a$V4=="sigmoid colon")
n3=which(a$V4=="cecum")
n4=which(a$V4=="hepatic flexure")
n5=which(a$V4=="splenic flexure")
n6=which(a$V4=="descending colon")
a[c(n1,n3,n4),5]='proximal'
a[c(n2,n5,n6),5]='distal'
nrow(a)
naid_df_1=a
#naid_df_1包含患者id，tcga_id，病人的肿瘤发生位置
#如果我们后续补上基因表达信息后，就可以分析出图了。





########################
#现在分析--肿瘤--中的基因表达作图。
expr_df=expr_tumour
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

#################################
gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000165556'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831
head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在cdx在--左右结肠癌--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_tumour=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('cdx2_tumor.pdf',10,10)
print(figure)
dev.off()
#############################


##################################
#正常组织中基因表达做图。
expr_df=expr_normal
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000165556'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831

head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在gene_name在左右结肠--正常组织--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_normal=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('cdx2_normal.pdf',10,10)
print(figure)
dev.off()




##################################
#现在把tumor 和 normal 的数据整合在一起画图
gene_name_tumour$v4='tumour'
gene_name_normal$v4='normal'

gene_name_total=rbind(gene_name_tumour,gene_name_normal)
nrow(gene_name_total)

library(stringr)
gene_name_total$v5=str_c(gene_name_total$location,'+',gene_name_total$v4)
gene_name=gene_name_total[,c(2,5)]
names(gene_name)[2]='location'
nrow(gene_name)
##保存gene_name在各个样本中的表达数据
write.csv(gene_name,file='CDX2_expression.csv',quote=F,row.names = T)


library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
unique(gene_name$location)
gene_name$location <- factor(gene_name$location, levels = c("proximal+normal",'proximal+tumour',"distal+normal",'distal+tumour'))
levels(gene_name$location)
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9",'99FFFF','7bbb5e'))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('cdx2_normalAndTumor.pdf',10,10)
print(figure)
dev.off()



#######################################################
#图画完了，是不是要做方差分析了呢？
gene_name
nrow(gene_name)
data=gene_name
names(data)=c('weight','group')
data$weight
data$group


#data <- PlantGrowth
library("ggpubr")
ggline(data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       #order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")


group_data <- split(data[,1], data[,2])
group_data
unlist(lapply(group_data, function(x){
  shapiro.test(x)$p.value
}))
#shapiro.test(x)是检验数据是否服从正态分布。H0为服从，H1为不服从。
#如果p值小于0.05，数据就是不服从。
#x可以是数值型向量，允许存在NA，但是非丢失数据需要在3-5000内
#proximal+tumour组的p值0.0001307421。这组数据不服从。
#问题来了，如果服从，能用方差分析么？
#我们暂时假设满足方差分析的常规条件，按常规方法来进行方差分析。

library(car)
qqPlot(group_data[[1]])


leveneTest(weight~group, data = data)
#leveneTest()方差齐性检验。
#H0方差齐，H1方差不齐。p=0.03086. p<0.05,方差不齐。

outlierTest(lm(weight~group, data=data))
#检测异常值。Nevada被判定为离群点。
#好吧，现在发现数据不全满足正态分布，方差不齐，有离群数据。
#怎么破？ 
#你要相信数据的真实性，你要根据数据选择合适的统计学方法。


#假设上述数据能用方差分析
aov=aov(weight~group, data = data)
#aov就是方差分析。
#aov是一个list数据，里面有各种各样的信息。
aov
summary(aov)
#可以查看方差分析的结果
#但这个方法还不能大组内的两两检验分析


#install.packages("agricolae")
library(agricolae)
bon <- LSD.test(aov,"group", p.adj="bonferroni")
#多组分析中的两两检验
bon$groups
plot(bon)
#你要意识到，数据的类型及想回答的问题，决定了什么样的分析方法才是
#合适的分析方法。



#####################
data
proximal_normal=data$weight[which(data$group == 'proximal+normal')]
proximal_normal
distal_normal=data$weight[which(data$group == 'distal+normal')]
proximal_tumour=data$weight[which(data$group == 'proximal+tumour')]
distal_tumour=data$weight[which(data$group == 'distal+tumour')]

t.test(proximal_normal,distal_normal)
t.test(proximal_normal,proximal_tumour) #*
t.test(distal_normal,distal_tumour)
t.test(proximal_tumour,distal_tumour) #*
#T检验分析，仅做测试用。

###上述分析方法是基于数据服从正态分布，方差齐。
#但数据实际上是不服从正态分布，且方差不齐。
#解决办法，数据不服从正态分布，仍可以行方差分析。
#welch's anova分析
oneway.test(weight~group, data = data, var.equal = F)
#但此办法不能做两两检验


#方差不齐，组内两两比较 可以用Games-Howell 分析。需要自建一个函数：
# Tukey の方法による多重比較
# Games-Howell の方法も選択できるように拡張 2009/08/03
oneway.test(weight~group, data = data, var.equal = F)

tukey <- function(	data,					# 観察値ベクトル
                   group,					# 群変数ベクトル
                   method=c("Tukey", "Games-Howell"))	# 手法の選択
{
  OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)			# 各群のケース数
  a <- length(n)						# 群の数
  phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
  Mean <- tapply(data, group, mean)			# 各群の平均値
  Variance <- tapply(data, group, var)			# 各群の不偏分散
  result1 <- cbind(n, Mean, Variance)			# 各群の統計量
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey の方法
    v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
    t <- combn(a, 2, function(ij)			# 対比較
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
    Tukey <- cbind(t, p)					# 対比較の結果
    rownames(Tukey) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
  }
  else {							# Games-Howell の方法
    t.df <- combn(a, 2, function(ij) {		# 対比較
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	# 有意確率を計算する
    Games.Howell <- cbind(t, df, p)			# 対比較の結果
    rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Games.Howell=Games.Howell))
  }
}

#进行两两比较
tukey(data$weight, data$group, method="G")
#注意数据的输入格式。

#在这个时间，一套完美的数据分析过程已顺利完成。

test=tukey(data$weight, data$group, method="G")
test
test1=test[["Games.Howell"]]
test1 = as.data.frame(test1)
test1
write.csv(test1,file="cdx2.csv",quote=F,row.names = T)
data.frame(levels(gene_name$location))
write.csv(test1,file="cdx2_group_index.csv",quote=F,row.names = T)



######################################



###############---HNF4A----#################
########################
#现在分析--肿瘤--中的基因表达作图。
expr_df=expr_tumour
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

#################################
gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000101076'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831
head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在cdx在--左右结肠癌--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_tumour=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('HNF4A_tumor.pdf',10,10)
print(figure)
dev.off()
#############################


##################################
#正常组织中基因表达做图。
expr_df=expr_normal
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000101076'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831

head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在gene_name在左右结肠--正常组织--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_normal=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('HNF4A_normal.pdf',10,10)
print(figure)
dev.off()




##################################
#现在把tumor 和 normal 的数据整合在一起画图
gene_name_tumour$v4='tumour'
gene_name_normal$v4='normal'

gene_name_total=rbind(gene_name_tumour,gene_name_normal)
nrow(gene_name_total)

library(stringr)
gene_name_total$v5=str_c(gene_name_total$location,'+',gene_name_total$v4)
gene_name=gene_name_total[,c(2,5)]
names(gene_name)[2]='location'
nrow(gene_name)


library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
unique(gene_name$location)
gene_name$location <- factor(gene_name$location, levels = c("proximal+normal",'proximal+tumour',"distal+normal",'distal+tumour'))
levels(gene_name$location)
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9",'99FFFF','7bbb5e'))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure

write.csv(gene_name,file='HNF4A_expression.csv',quote=F,row.names = T)
pdf('HNF4A_normalAndTumor.pdf',10,10)
print(figure)
dev.off()



#######################################################
#图画完了，是不是要做方差分析了呢？
gene_name
nrow(gene_name)
data=gene_name
names(data)=c('weight','group')
data$weight
data$group


#data <- PlantGrowth
library("ggpubr")
ggline(data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       #order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")


group_data <- split(data[,1], data[,2])
group_data
unlist(lapply(group_data, function(x){
  shapiro.test(x)$p.value
}))
#shapiro.test(x)是检验数据是否服从正态分布。H0为服从，H1为不服从。
#如果p值小于0.05，数据就是不服从。
#x可以是数值型向量，允许存在NA，但是非丢失数据需要在3-5000内
#proximal+tumour, distal+tumour 的p值 < 0.05。这组数据不服从。
#一个有意思的发现。正常组织HNF4A的表达都满足正态分布
#但肿瘤组织的表达都不满足正态分布。
#问题来了，如果服从，能用方差分析么？
#我们暂时假设满足方差分析的常规条件，按常规方法来进行方差分析。

library(car)
qqPlot(group_data[[1]])


leveneTest(weight~group, data = data)
#leveneTest()方差齐性检验。
#H0方差齐，H1方差不齐。p=0.0001792. p<0.05,方差不齐。

outlierTest(lm(weight~group, data=data))
#检测异常值。Nevada被判定为离群点。
#好吧，现在发现数据不全满足正态分布，方差不齐，有离群数据。
#怎么破？ 
#你要相信数据的真实性，你要根据数据选择合适的统计学方法。


#假设上述数据能用方差分析。
#虽然，实际上这里并不能用。
aov=aov(weight~group, data = data)
#aov就是方差分析。
#aov是一个list数据，里面有各种各样的信息。
aov
summary(aov)
#可以查看方差分析的结果
#但这个方法还不能大组内的两两检验分析


#install.packages("agricolae")
library(agricolae)
bon <- LSD.test(aov,"group", p.adj="bonferroni")
#多组分析中的两两检验
bon$groups
plot(bon)
#你要意识到，数据的类型及想回答的问题，决定了什么样的分析方法才是
#合适的分析方法。



#####################
data
proximal_normal=data$weight[which(data$group == 'proximal+normal')]
proximal_normal
distal_normal=data$weight[which(data$group == 'distal+normal')]
proximal_tumour=data$weight[which(data$group == 'proximal+tumour')]
distal_tumour=data$weight[which(data$group == 'distal+tumour')]

t.test(proximal_normal,distal_normal)
t.test(proximal_normal,proximal_tumour) #*
t.test(distal_normal,distal_tumour)
t.test(proximal_tumour,distal_tumour) #*
#T检验分析，仅做测试用。

###上述分析方法是基于数据服从正态分布，方差齐。
#但数据实际上是不服从正态分布，且方差不齐。
#解决办法，数据不服从正态分布，仍可以行方差分析。
#welch's anova分析
oneway.test(weight~group, data = data, var.equal = F)
#但此办法不能做两两检验


#方差不齐，组内两两比较 可以用Games-Howell 分析。需要自建一个函数：
# Tukey の方法による多重比較
# Games-Howell の方法も選択できるように拡張 2009/08/03
oneway.test(weight~group, data = data, var.equal = F)

tukey <- function(	data,					# 観察値ベクトル
                   group,					# 群変数ベクトル
                   method=c("Tukey", "Games-Howell"))	# 手法の選択
{
  OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)			# 各群のケース数
  a <- length(n)						# 群の数
  phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
  Mean <- tapply(data, group, mean)			# 各群の平均値
  Variance <- tapply(data, group, var)			# 各群の不偏分散
  result1 <- cbind(n, Mean, Variance)			# 各群の統計量
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey の方法
    v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
    t <- combn(a, 2, function(ij)			# 対比較
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
    Tukey <- cbind(t, p)					# 対比較の結果
    rownames(Tukey) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
  }
  else {							# Games-Howell の方法
    t.df <- combn(a, 2, function(ij) {		# 対比較
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	# 有意確率を計算する
    Games.Howell <- cbind(t, df, p)			# 対比較の結果
    rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Games.Howell=Games.Howell))
  }
}

#进行两两比较
tukey(data$weight, data$group, method="G")
#注意数据的输入格式。

#在这个时间，一套完美的数据分析过程已顺利完成。

test=tukey(data$weight, data$group, method="G")
test
test1=test[["Games.Howell"]]
test1 = as.data.frame(test1)
test1
write.csv(test1,file="HNF4A.csv",quote=F,row.names = T)
data.frame(levels(gene_name$location))
write.csv(test1,file="HNF4A_group_index.csv",quote=F,row.names = T)




###############---SATB2----#################
########################
#现在分析--肿瘤--中的基因表达作图。
expr_df=expr_tumour
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

#################################
gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000119042'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831
head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在cdx在--左右结肠癌--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_tumour=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('SATB2_tumor.pdf',10,10)
print(figure)
dev.off()
#############################


##################################
#正常组织中基因表达做图。
expr_df=expr_normal
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000119042'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831

head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在gene_name在左右结肠--正常组织--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_normal=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('SATB2_normal.pdf',10,10)
print(figure)
dev.off()




##################################
#现在把tumor 和 normal 的数据整合在一起画图
gene_name_tumour$v4='tumour'
gene_name_normal$v4='normal'

gene_name_total=rbind(gene_name_tumour,gene_name_normal)
nrow(gene_name_total)

library(stringr)
gene_name_total$v5=str_c(gene_name_total$location,'+',gene_name_total$v4)
gene_name=gene_name_total[,c(2,5)]
names(gene_name)[2]='location'
nrow(gene_name)


library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
unique(gene_name$location)
gene_name$location <- factor(gene_name$location, levels = c("proximal+normal",'proximal+tumour',"distal+normal",'distal+tumour'))
levels(gene_name$location)
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9",'99FFFF','7bbb5e'))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure

write.csv(gene_name,file='SATB2_expression.csv',quote=F,row.names = T)
pdf('SATB2_normalAndTumor.pdf',10,10)
print(figure)
dev.off()



#######################################################
#图画完了，是不是要做方差分析了呢？
gene_name
nrow(gene_name)
data=gene_name
names(data)=c('weight','group')
data$weight
data$group


#data <- PlantGrowth
library("ggpubr")
ggline(data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       #order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")


group_data <- split(data[,1], data[,2])
group_data
unlist(lapply(group_data, function(x){
  shapiro.test(x)$p.value
}))
#shapiro.test(x)是检验数据是否服从正态分布。H0为服从，H1为不服从。
#如果p值小于0.05，数据就是不服从。
#x可以是数值型向量，允许存在NA，但是非丢失数据需要在3-5000内
#proximal+tumour, distal+tumour 的p值 < 0.05。这组数据不服从。
#一个有意思的发现。正常组织SATB2的表达都满足正态分布
#但肿瘤组织的表达都不满足正态分布。
#问题来了，如果服从，能用方差分析么？
#我们暂时假设满足方差分析的常规条件，按常规方法来进行方差分析。

library(car)
qqPlot(group_data[[1]])


leveneTest(weight~group, data = data)
#leveneTest()方差齐性检验。
#H0方差齐，H1方差不齐。p=0.0001792. p<0.05,方差不齐。

outlierTest(lm(weight~group, data=data))
#检测异常值。Nevada被判定为离群点。
#好吧，现在发现数据不全满足正态分布，方差不齐，有离群数据。
#怎么破？ 
#你要相信数据的真实性，你要根据数据选择合适的统计学方法。


#假设上述数据能用方差分析。
#虽然，实际上这里并不能用。
aov=aov(weight~group, data = data)
#aov就是方差分析。
#aov是一个list数据，里面有各种各样的信息。
aov
summary(aov)
#可以查看方差分析的结果
#但这个方法还不能大组内的两两检验分析


#install.packages("agricolae")
library(agricolae)
bon <- LSD.test(aov,"group", p.adj="bonferroni")
#多组分析中的两两检验
bon$groups
plot(bon)
#你要意识到，数据的类型及想回答的问题，决定了什么样的分析方法才是
#合适的分析方法。



#####################
data
proximal_normal=data$weight[which(data$group == 'proximal+normal')]
proximal_normal
distal_normal=data$weight[which(data$group == 'distal+normal')]
proximal_tumour=data$weight[which(data$group == 'proximal+tumour')]
distal_tumour=data$weight[which(data$group == 'distal+tumour')]

t.test(proximal_normal,distal_normal)
t.test(proximal_normal,proximal_tumour) #*
t.test(distal_normal,distal_tumour)
t.test(proximal_tumour,distal_tumour) #*
#T检验分析，仅做测试用。

###上述分析方法是基于数据服从正态分布，方差齐。
#但数据实际上是不服从正态分布，且方差不齐。
#解决办法，数据不服从正态分布，仍可以行方差分析。
#welch's anova分析
oneway.test(weight~group, data = data, var.equal = F)
#但此办法不能做两两检验


#方差不齐，组内两两比较 可以用Games-Howell 分析。需要自建一个函数：
# Tukey の方法による多重比較
# Games-Howell の方法も選択できるように拡張 2009/08/03
oneway.test(weight~group, data = data, var.equal = F)

tukey <- function(	data,					# 観察値ベクトル
                   group,					# 群変数ベクトル
                   method=c("Tukey", "Games-Howell"))	# 手法の選択
{
  OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)			# 各群のケース数
  a <- length(n)						# 群の数
  phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
  Mean <- tapply(data, group, mean)			# 各群の平均値
  Variance <- tapply(data, group, var)			# 各群の不偏分散
  result1 <- cbind(n, Mean, Variance)			# 各群の統計量
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey の方法
    v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
    t <- combn(a, 2, function(ij)			# 対比較
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
    Tukey <- cbind(t, p)					# 対比較の結果
    rownames(Tukey) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
  }
  else {							# Games-Howell の方法
    t.df <- combn(a, 2, function(ij) {		# 対比較
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	# 有意確率を計算する
    Games.Howell <- cbind(t, df, p)			# 対比較の結果
    rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Games.Howell=Games.Howell))
  }
}

#进行两两比较
tukey(data$weight, data$group, method="G")
#注意数据的输入格式。

#在这个时间，一套完美的数据分析过程已顺利完成。

test=tukey(data$weight, data$group, method="G")
test
test1=test[["Games.Howell"]]
test1 = as.data.frame(test1)
test1
write.csv(test1,file="SATB2.csv",quote=F,row.names = T)
data.frame(levels(gene_name$location))
write.csv(test1,file="SATB2_group_index.csv",quote=F,row.names = T)


#########################################################################
###############---VIL1----#################
########################
#现在分析--肿瘤--中的基因表达作图。
expr_df=expr_tumour
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

#################################
gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000127831'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831
head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在cdx在--左右结肠癌--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_tumour=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('VIL1_tumor.pdf',10,10)
print(figure)
dev.off()
#############################


##################################
#正常组织中基因表达做图。
expr_df=expr_normal
expr_df$gene_id=substr(expr_df$gene_id,1,15)
#处理一下expr_df$gene_id的字符串。
head(expr_df$gene_id)

gene_name=expr_df[which(expr_df$gene_id[] == 'ENSG00000127831'),]
#'ENSG00000165556'对应的是CDX2
# HNF4A ENSG00000101076 
# SATB2 ENSG00000119042
# VIL1 ENSG00000127831

head(gene_name)
gene_name=t(gene_name)
head(gene_name)
rownames(gene_name)
gene_name=cbind(rownames(gene_name),gene_name)
colnames(gene_name)=c('patient_id','gene_name')
head(gene_name)
gene_name=gene_name[-1,]
head(gene_name)
gene_name=data.frame(gene_name)
head(gene_name)
#把gene_name的表达量调整成数据框的格式，以方便后续的操作。


for(i in 1: nrow(gene_name)){
  for(j in 1:nrow(naid_df_1)){
    if(gene_name$patient_id[i] == naid_df_1$TCGA_id[j]){gene_name$location[i] = naid_df_1$V5[j]}
  }
}
head(gene_name)
#现在gene_name在左右结肠--正常组织--的表达的数据表已经准备好了。现在就可以做图了啊。
gene_name_normal=gene_name

library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
gene_name$location <- factor(gene_name$location, levels = c("proximal","distal"))
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9"))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure


pdf('VIL1_normal.pdf',10,10)
print(figure)
dev.off()




##################################
#现在把tumor 和 normal 的数据整合在一起画图
gene_name_tumour$v4='tumour'
gene_name_normal$v4='normal'

gene_name_total=rbind(gene_name_tumour,gene_name_normal)
nrow(gene_name_total)

library(stringr)
gene_name_total$v5=str_c(gene_name_total$location,'+',gene_name_total$v4)
gene_name=gene_name_total[,c(2,5)]
names(gene_name)[2]='location'
nrow(gene_name)


library(ggplot2)
#gene_name$location <- as.factor(gene_name$location)
#levels(gene_name$location)=c('1','2')
#这个设置proximal 和 distal顺序变后，结果可能会相反。一定要慎重考虑。
unique(gene_name$location)
gene_name$location <- factor(gene_name$location, levels = c("proximal+normal",'proximal+tumour',"distal+normal",'distal+tumour'))
levels(gene_name$location)
#上述方法不靠谱，换这个。运行成功。

gene_name$gene_name=as.numeric(gene_name$gene_name)
p <- ggplot(gene_name, aes(x = location, y = gene_name))
p
p + geom_violin ()
p + geom_violin(aes(fill = location), trim = FALSE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "red")

p2 <- p + geom_violin(aes(fill = location), trim = FALSE) + theme_minimal()
p2
p3=p2 + scale_fill_manual(values=c("#66CCCC", "#56B4E9",'99FFFF','7bbb5e'))
p3
p4=p3 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+ 
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.5,colour="black"))
p4
figure=p4
figure

write.csv(gene_name,file='VIL1_expression.csv',quote=F,row.names = T)
pdf('VIL1_normalAndTumor.pdf',10,10)
print(figure)
dev.off()



#######################################################
#图画完了，是不是要做方差分析了呢？
gene_name
nrow(gene_name)
data=gene_name
names(data)=c('weight','group')
data$weight
data$group


#data <- PlantGrowth
library("ggpubr")
ggline(data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       #order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")


group_data <- split(data[,1], data[,2])
group_data
unlist(lapply(group_data, function(x){
  shapiro.test(x)$p.value
}))
#shapiro.test(x)是检验数据是否服从正态分布。H0为服从，H1为不服从。
#如果p值小于0.05，数据就是不服从。
#x可以是数值型向量，允许存在NA，但是非丢失数据需要在3-5000内
#proximal+tumour, distal+tumour 的p值 < 0.05。这组数据不服从。
#一个有意思的发现。正常组织VIL1的表达都满足正态分布
#但肿瘤组织的表达都不满足正态分布。
#问题来了，如果服从，能用方差分析么？
#我们暂时假设满足方差分析的常规条件，按常规方法来进行方差分析。

library(car)
qqPlot(group_data[[1]])


leveneTest(weight~group, data = data)
#leveneTest()方差齐性检验。
#H0方差齐，H1方差不齐。p=0.1017. p>0.05,方差是齐的。哈哈。

outlierTest(lm(weight~group, data=data))
#检测异常值。Nevada被判定为离群点。
#好吧，现在发现数据不全满足正态分布，方差齐，有离群数据。
#怎么破？ 
#你要相信数据的真实性，你要根据数据选择合适的统计学方法。


#假设上述数据能用方差分析。
#虽然，实际上这里并不能用。
aov=aov(weight~group, data = data)
#aov就是方差分析。
#aov是一个list数据，里面有各种各样的信息。
aov
summary(aov)
#可以查看方差分析的结果
#但这个方法还不能大组内的两两检验分析


#install.packages("agricolae")
library(agricolae)
bon <- LSD.test(aov,"group", p.adj="bonferroni")
#多组分析中的两两检验
bon$groups
plot(bon)
#你要意识到，数据的类型及想回答的问题，决定了什么样的分析方法才是
#合适的分析方法。



#####################
data
proximal_normal=data$weight[which(data$group == 'proximal+normal')]
proximal_normal
distal_normal=data$weight[which(data$group == 'distal+normal')]
proximal_tumour=data$weight[which(data$group == 'proximal+tumour')]
distal_tumour=data$weight[which(data$group == 'distal+tumour')]

t.test(proximal_normal,distal_normal)
t.test(proximal_normal,proximal_tumour) #*
t.test(distal_normal,distal_tumour)
t.test(proximal_tumour,distal_tumour) #*
#T检验分析，仅做测试用。

###上述分析方法是基于数据服从正态分布，方差齐。
#但数据实际上是不服从正态分布，但--方差是齐的--。
#解决办法，数据不服从正态分布，仍可以行方差分析。
#welch's anova分析
oneway.test(weight~group, data = data, var.equal = F)
#但此办法不能做两两检验



#方差不齐，组内两两比较 可以用Games-Howell 分析。需要自建一个函数：
# Tukey の方法による多重比較
# Games-Howell の方法も選択できるように拡張 2009/08/03
oneway.test(weight~group, data = data, var.equal = F)

tukey <- function(	data,					# 観察値ベクトル
                   group,					# 群変数ベクトル
                   method=c("Tukey", "Games-Howell"))	# 手法の選択
{
  OK <- complete.cases(data, group)			# 欠損値を持つケースを除く
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)			# 各群のケース数
  a <- length(n)						# 群の数
  phi.e <- sum(n)-a					# 誤差分散（群内不偏分散）の自由度
  Mean <- tapply(data, group, mean)			# 各群の平均値
  Variance <- tapply(data, group, var)			# 各群の不偏分散
  result1 <- cbind(n, Mean, Variance)			# 各群の統計量
  rownames(result1) <- paste("Group", 1:a, sep="")
  method <- match.arg(method)
  if (method == "Tukey") {				# Tukey の方法
    v.e <- sum((n-1)*Variance)/phi.e		# 誤差分散（群内不偏分散）
    t <- combn(a, 2, function(ij)			# 対比較
      abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
    p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)	# 有意確率を計算する
    Tukey <- cbind(t, p)					# 対比較の結果
    rownames(Tukey) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
  }
  else {							# Games-Howell の方法
    t.df <- combn(a, 2, function(ij) {		# 対比較
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)	# 有意確率を計算する
    Games.Howell <- cbind(t, df, p)			# 対比較の結果
    rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
    return(list(result1=result1, Games.Howell=Games.Howell))
  }
}

#进行两两比较
tukey(data$weight, data$group, method="G")
#注意数据的输入格式。

#在这个时间，一套完美的数据分析过程已顺利完成。

test=tukey(data$weight, data$group, method="G")
test
test1=test[["Games.Howell"]]
test1 = as.data.frame(test1)
test1
write.csv(test1,file="VIL1.csv",quote=F,row.names = T)
data.frame(levels(gene_name$location))
write.csv(test1,file="VIL1_group_index.csv",quote=F,row.names = T)

##################################################
#这套代码也太长了。1600多行。





##########################correlationship analysis################
#setwd('E:/R_study/用TCGA数据库做分析_20211219/figure_4')
getwd()
temp=list.files(pattern='*.csv')
temp
temp=c("CDX2_expression.csv","HNF4A_expression.csv" , "SATB2_expression.csv",  "VIL1_expression.csv" )

a=list()
for(i in 1:length(temp)){
  a[[i]]=read.csv(file=temp[i], header=T,check.names = FALSE)
}
names(a)=temp
a

for(i in 1:length(temp)){
  names(a[[i]])[1]='sample'
  names(a[[i]])[2]=temp[i]
}

names(a[[1]])


library(dplyr)
data=a[[1]]
for(i in 2:length(temp)){
  data=dplyr::inner_join(data,a[[i]],by ="sample")  
}
data
class(data)
names(data)

#column=c('sample',"ENSG00000165556_CDX2.csv", "ENSG00000101076_HNF4A.csv","ENSG00000135100_HNF1A.csv",
#"ENSG00000119042_SATB2.csv","ENSG00000105996_HOXA2.csv",'group')

column=c('sample',"CDX2_expression.csv","HNF4A_expression.csv","SATB2_expression.csv","VIL1_expression.csv","location.y.y")

column
names(data)


data_1=dplyr::select(data,column)

head(data_1)

####################################################
#筛选数据，做相关性分析了。
#所有代码，只需要改下面的分组参数就可以了。其他的不变。
data_2=data_1[which(data_1$"location.y.y" == 'proximal+tumour'),]
data_2
names(data_2)
data_2=dplyr::select(data_2,names(data_2)[2:5])
data_2

names(data_2)=c('CDX2','HNF4A','SATB2','VIL1')

#names(data_2)=c('Cdx2','Hnf4a','Hnf1a','Satb2')



#cor() contained by R language itself
#cor(data_2,y=NULL, use='everything',method=c('pearson','kendall','spearman'))
#default format
cor(data_2,y=NULL, use='everything',method=c('pearson'))
cor(data_2,method=c('pearson'))
cor(data_2)
#The results analyzed by these three commond are the same. 


library(Hmisc)
cor=rcorr(as.matrix(data_2))
#This function applies 'pearson method'
cor$r
cor$n
cor$P


#对相关性作图
#BiocManager::install("tidyverse")
#BiocManager::install("circlize")
library(tidyverse)
library(corrplot)
library(circlize)
corrplot(cor$r,method='circle',type='upper')

pdf('correlation-ship.pdf',10,10)
print(corrplot(cor$r,method='circle',type='upper'))
dev.off()


#BiocManager::install("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(data_2, histogram=TRUE, pch=19)


