# Code for Differential ChIP-seq peak analyses by DiffBind
## The following code is for comparing Proximal vs. Dista organoids in base medium. The same code is used for the various comparisons by selecting the appropriate files.

a='C:/Users/lijingyang/Desktop/R_test'
setwd(a)
getwd()
list.files()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DiffBind")
#browseVignettes("DiffBind")
library(DiffBind)

samples <- read.csv(file="Comparison_Groups.csv", header=T)
samples
names(samples)

Comparison_Groups=dba(sampleSheet="Comparison_Groups.csv", dir=a)
Comparison_Groups

plot(Comparison_Groups)


peakdata <- dba.show(Comparison_Groups)$Intervals

peakdata

Comparison_Groups <- dba.blacklist(Comparison_Groups, blacklist=DBA_BLACKLIST_MM10,
                           greylist=FALSE)

peakdata.BL <- dba.show(Comparison_Groups)$Intervals
peakdata.BL

peakdata - peakdata.BL

Comparison_Groups <- dba.count(Comparison_Groups)
backup_1=Comparison_Groups

# Reads, FRiP.
info <- dba.show(Comparison_Groups)
info

libsizes= cbind(LibReads=info$Reads, FRiP=info$FRip, PeakReads=round(info$Reads*info$FRiP))
rownames(libsizes)=info$ID
libsizes

Comparison_Groups=dba.normalize(Comparison_Groups)

Comparison_Groups=dba.normalize(Comparison_Groups)
norm=dba.normalize(Comparison_Groups, bRetrieve=TRUE)
norm
normlibs=cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,NormLibSize=round(norm$lib.sizes/norm$norm.factors))

rownames(normlibs)=info$ID

normlibs

Comparison_Groups=dba.contrast(Comparison_Groups,reorderMeta=list(Condition='distal'))


Comparison_Groups
Comparison_Groups=dba.analyze(Comparison_Groups)


dba.show(Comparison_Groups, bContrasts=TRUE)
plot(Comparison_Groups)

Comparison_Groups.DB=dba.report(Comparison_Groups)
Comparison_Groups.DB

sum(Comparison_Groups.DB$Fold>0)
sum(Comparison_Groups.DB$Fold<0)



dba.plotVenn(Comparison_Groups,contrast=1,bDB=TRUE,bGain=TRUE,bLoss=TRUE,bAll=FALSE)

figure_1 = dba.plotPCA(Comparison_Groups,DBA_TREATMENT,label=DBA_CONDITION)
figure_2 = dba.plotPCA(Comparison_Groups,DBA_TREATMENT,label=DBA_REPLICATE)
figure_1
figure_2


dba.plotPCA(Comparison_Groups,contrast=1,label=DBA_TISSUE)
dba.plotPCA(Comparison_Groups,contrast=1, b3D=TRUE)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("rgl")


dba.plotPCA(Comparison_Groups,contrast=1, b3D=TRUE)

dba.plotVolcano(Comparison_Groups)

corvals <- dba.plotHeatmap(Comparison_Groups)


hmap <- colorRampPalette(c("red", "black", "green"))(n = 6)

readscores <- dba.plotHeatmap(Comparison_Groups, contrast=1, correlations=FALSE,scale="row", colScheme = hmap)

Comparison_Groups.DB

write.table(Comparison_Groups.DB,file='abc.csv',sep=',',row.names=FALSE)



#BiocManager::install("profileplyr")
library(profileplyr)

profiles <- dba.plotProfile(Comparison_Groups, merge=c(DBA_REPLICATE))
dba.plotProfile(profiles)

figure=dba.plotProfile(profiles)

pdf('abc.pdf',10,10)
dba.plotProfile(profiles)
dev.off()
