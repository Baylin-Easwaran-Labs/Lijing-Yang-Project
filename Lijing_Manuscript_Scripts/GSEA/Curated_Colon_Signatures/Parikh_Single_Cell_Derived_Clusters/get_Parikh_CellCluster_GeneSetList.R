rm(list=ls())
library(readxl)
library(GSEABase)

Parikh_Single_Cell_Derived_CellCluster_Signature <- read_xlsx("D:/HariEaswaran/R_Working_Directory_onc-cbio/Required_Files/GSEA/Curated_Colon_Signatures/Parikh_Single_Cell_Derived_Clusters/Parikh_41586_2019_992_MOESM3_ESM.xlsx", sheet=2, skip=1, col_names = TRUE)
Parikh_Single_Cell_Derived_CellCluster_Signature <- Parikh_Single_Cell_Derived_CellCluster_Signature[,1:7]
colnames(Parikh_Single_Cell_Derived_CellCluster_Signature) <- c("AUC", "power", "avg_logFC", "pct.1", "pct.2", "cluster", "gene")

dim(Parikh_Single_Cell_Derived_CellCluster_Signature)
unique(Parikh_Single_Cell_Derived_CellCluster_Signature$cluster)

Parikh_Single_Cell_Derived_CellCluster_Signature <- Parikh_Single_Cell_Derived_CellCluster_Signature[-which(Parikh_Single_Cell_Derived_CellCluster_Signature$cluster %in% NA), ]

Parikh_CellCluster_GeneSetList <- lapply(unique(Parikh_Single_Cell_Derived_CellCluster_Signature$cluster), FUN=function(l){
  xx <- as.data.frame(Parikh_Single_Cell_Derived_CellCluster_Signature)
  xx <- xx[xx$avg_logFC > 0, ] # get positive markers of cell clusters
  return(as.vector(xx[xx$cluster %in% l, "gene"]))
})
names(Parikh_CellCluster_GeneSetList) <- unique(Parikh_Single_Cell_Derived_CellCluster_Signature$cluster)
sapply(Parikh_CellCluster_GeneSetList, length)

########################
# Create gmt file
########################
require(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Parikh_CellCluster_GeneSetList <- lapply(Parikh_CellCluster_GeneSetList, FUN=function(l){
 xx <- getBM(attributes = c('entrezgene_id'),
        filters = 'hgnc_symbol',
        values = l, 
        mart = human)
 return(as.character(xx$entrezgene_id))
})

sapply(Parikh_CellCluster_GeneSetList, class)
sapply(Parikh_CellCluster_GeneSetList, length)


# Convert the list to a GeneSetCollection using GSEABase functionality
myGeneSetList <- list()
for (i in 1:length(Parikh_CellCluster_GeneSetList)) {
  myGeneSetList[[i]] <- GeneSet(Parikh_CellCluster_GeneSetList[[i]], setName = names(Parikh_CellCluster_GeneSetList)[i])
}

Parikh_CellCluster_GeneSetList <- GeneSetCollection(myGeneSetList)

# Convert GeneSetCollection to gmt format
toGmt(Parikh_CellCluster_GeneSetList, con = "D:/HariEaswaran/R_Working_Directory_onc-cbio/Required_Files/GSEA/Curated_Colon_Signatures/Parikh_Single_Cell_Derived_Clusters/Parikh_CellCluster_GeneSetList.gmt")
