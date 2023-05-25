###################################
#Load libaries and data
##Order of loading libraries and attaching objects important
####################################
library(biomaRt)
library(ComplexHeatmap)
library(data.table)
library(AnnotationDbi)
library(org.Mm.eg.db)


attach("/RData/Lijing_RNA_seq_Analyses_062922.RData", name="Lijing_RNA_seq_Analyses_062922")
ls("Lijing_RNA_seq_Analyses_062922")

attach("/RData/Peak_Calls_Analysis.RData", name="Peak_Calls_Analysis")
ls("Peak_Calls_Analysis")


dim(annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig)
head(annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig)

annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig[annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig$SYMBOL %in% "Satb2", ]
annotatedRes_Prox.Cdx2KOVsWT_FM_DF_sig[annotatedRes_Prox.Cdx2KOVsWT_FM_DF_sig$SYMBOL %in% "Satb2", ]

annotatedRes_Prox.Cdx2KOVsWT_BM_DF[annotatedRes_Prox.Cdx2KOVsWT_BM_DF$SYMBOL %in% "Satb2", ]

annotatedRes_Prox.Cdx2KOVsWT_BM_DF[annotatedRes_Prox.Cdx2KOVsWT_BM_DF$SYMBOL %in% "Hnf4a", ]


dim(annotatedRes_Prox.Cdx2KOVsWT_BM_DF[annotatedRes_Prox.Cdx2KOVsWT_BM_DF$padj <= 0.05, ])

dim(annotatedRes_Prox.Cdx2KOVsWT_BM_DF[annotatedRes_Prox.Cdx2KOVsWT_BM_DF$newPadj <= 0.05, ])

dim(annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig)

identical(annotatedRes_Prox.Cdx2KOVsWT_BM_DF[annotatedRes_Prox.Cdx2KOVsWT_BM_DF$padj <= 0.05, ], annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig)

############################
# Complexheatmap
############################

## Get normalized counts (from normCounts) values for genes that are identified as significantly differential expressed in proximal Cdx2-KO

normCounts.DF.X <- as.data.frame(normCounts) #normCounts.DF.X instead of normCounts.DF is used here because DESEq2 does not filter low count values for diff. expr. analyses
normCounts.DF.X <- normCounts.DF.X + 1 #because many counts values are 0
dim(normCounts.DF.X)
normCounts.DF.X <- log2(normCounts.DF.X) #log transform counts

eToSym.X <- select(org.Mm.eg.db,
                 keys = rownames(normCounts.DF.X),
                 keytype = "ENTREZID",
                 columns="SYMBOL")
normCounts.DF.X <- merge(eToSym.X,normCounts.DF.X,
                       by.x=1,
                       by.y=0,
                       all.x=FALSE,
                       all.y=TRUE)

#colnames(normCounts.DF)
colnames(normCounts.DF.X) <- gsub(pattern="_quant.sf", replacement="",colnames(normCounts.DF.X))

head(normCounts.DF.X)
dim(normCounts.DF.X)

normCounts_Prox.Cdx2KOVsWT_FM.BM_sig <- normCounts.DF.X[
  normCounts.DF.X$ENTREZID %in% union(annotatedRes_Prox.Cdx2KOVsWT_FM_DF_sig$ENTREZID, annotatedRes_Prox.Cdx2KOVsWT_BM_DF_sig$ENTREZID), ]

cbind(samples$Samples ,samples$Concise_Description, colnames(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig[,-c(1,2)]))
rownames(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig) <- normCounts_Prox.Cdx2KOVsWT_FM.BM_sig[,1]
colnames(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig)[3:18] <- samples$Concise_Description
head(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig)
dim(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig)

## Convert cdx2PositivelyCorrelatedGenes to mouse ID

cdx2PosCorrProx_mouse.geneID <- convertHumanToMouseEntrezGeneIDList(cdx2PositivelyCorrelatedGenes_Proximal.geneID)
cdx2PosCorrDist_mouse.geneID <- convertHumanToMouseEntrezGeneIDList(cdx2PositivelyCorrelatedGenes_Distal.geneID)

intersect(rownames(normCounts_Prox.Cdx2KOVsWT_FM.BM_sig), cdx2PosCorrProx_mouse.geneID[,"mouse_entrezgene_id"])


## Get genes that are significantly downregulated in Prox.Cdx2KOVsWT in BM/FM and that are postiviely correlated with Cdx2 expression in proximal colon cancers.
m <- normCounts_Prox.Cdx2KOVsWT_FM.BM_sig
goi <- c("Satb2", "Vil1", "Hnf4a")
goi.id <- m[m$SYMBOL %in% goi, "ENTREZID"]
goi.id <- intersect(goi.id, as.character(cdx2PosCorrProx_mouse.geneID[,"mouse_entrezgene_id"])) #get those genes that are also positively correlated with CDX2 in proximal colon cancers
m[goi.id, c("ENTREZID", "SYMBOL")]

#m <- m[,3:18]
m.scaled <- m[,grep("cdx2", colnames(m))] #retain only counts values of Cdx2-KO columns.
m.scaled <- t(scale(t(m.scaled)))
m.scaled <-  na.omit(m.scaled)
m <- m[rownames(m.scaled), ] #recreate m so that the rows are exactly same in m and m.scaled because this has to be same to create anno.mark and anno.labels.
identical(rownames(m.scaled), rownames(m))
anno.mark <- which(rownames(m.scaled) %in% goi.id)
anno.labels <- m[anno.mark, "SYMBOL"]
m[anno.mark, c("SYMBOL", "ENTREZID")]

which(is.na(m.scaled) == T)
ha = rowAnnotation(foo = anno_mark(at = anno.mark, labels = anno.labels))
pdf("Fig6.pdf", height=8, width=8)
Heatmap(m.scaled, name = "mat", cluster_rows = T, clustering_distance_rows = "euclidean", right_annotation = ha,
        row_names_side = "left", row_names_gp = gpar(fontsize = 4))
dev.off()


