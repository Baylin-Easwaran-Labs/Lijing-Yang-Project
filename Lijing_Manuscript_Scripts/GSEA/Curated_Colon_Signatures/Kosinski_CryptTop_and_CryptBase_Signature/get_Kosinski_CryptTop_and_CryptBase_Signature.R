rm(list=ls())
library(readxl)
require(biomaRt)
library(GSEABase)
################################################
# Get the CryptBottom signature ("367 Positive Significant Genes (Highly expressed in colon crypt)
################################################
Kosinski_CryptBottom_Signature <- read_xls("D:/HariEaswaran/R_Working_Directory_onc-cbio/Required_Files/GSEA/Curated_Colon_Signatures/Kosinski_CryptTop_and_CryptBase_Signature/Kosinski_CryptTop_and_CryptBase_Signature.xls", sheet=1, range = cell_rows(23:390), col_names = TRUE)

#Get gene Symbol from Kosinski data
Kosinski_CryptBottom.geneSymbol <- unname(
  sapply(Kosinski_CryptBottom_Signature$"Gene Name", FUN=function(i){
    xx <- strsplit(i, split="||", fixed=T)[[1]][1]
    xx <- strsplit(xx, split=" ")[[1]][1]
  })
)
Kosinski_CryptBottom.geneSymbol <- Kosinski_CryptBottom.geneSymbol[-which(Kosinski_CryptBottom.geneSymbol %in% NA)] 

#Get Entrez gene ID from Kosinski data
Kosinski_CryptBottom.geneID <- unname(
  sapply(Kosinski_CryptBottom_Signature$"Gene Name", FUN=function(i){
    xx <- strsplit(i, split="||", fixed=T)[[1]][6]
    xx <- strsplit(xx, split=" ")[[1]][2]
  })
)

#Get the entrezgene_id for the gene Symbols. This will be later merged with the geneIDs obtained from the Kosinski data (Kosinski_CryptBottom.geneID). This is done because some gene symbols do not have associated gene IDs, and others have gene IDs but not asocaitd gene symbols, in the Kosinski data. 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Kosinski_CryptBottom.geneSymbolTogeneID <- getBM(attributes = c('entrezgene_id'),
            filters = 'hgnc_symbol',
            values = Kosinski_CryptBottom.geneSymbol, 
            mart = human)
Kosinski_CryptBottom.geneSymbolTogeneID <- as.character(Kosinski_CryptBottom.geneSymbolTogeneID$entrezgene_id)


Kosinski_CryptBottom.conmbined_geneID <- unique(c(Kosinski_CryptBottom.geneID, Kosinski_CryptBottom.geneSymbolTogeneID))

length(Kosinski_CryptBottom.conmbined_geneID)

################################################
# Get the CryptTop signature (602 Negative Significant Genes (Highly expressed in colon top))
################################################
Kosinski_CryptTop_Signature <- read_xls("D:/HariEaswaran/R_Working_Directory_onc-cbio/Required_Files/GSEA/Curated_Colon_Signatures/Kosinski_CryptTop_and_CryptBase_Signature/Kosinski_CryptTop_and_CryptBase_Signature.xls", sheet=1, range = cell_rows(393:995), col_names = TRUE)

#Get gene Symbol from Kosinski data
Kosinski_CryptTop.geneSymbol <- unname(
  sapply(Kosinski_CryptTop_Signature$"Gene Name", FUN=function(i){
    xx <- strsplit(i, split="||", fixed=T)[[1]][1]
    xx <- strsplit(xx, split=" ")[[1]][1]
  })
)
Kosinski_CryptTop.geneSymbol <- Kosinski_CryptTop.geneSymbol[-which(Kosinski_CryptTop.geneSymbol %in% NA)] 

#Get Entrez gene ID from Kosinski data
Kosinski_CryptTop.geneID <- unname(
  sapply(Kosinski_CryptTop_Signature$"Gene Name", FUN=function(i){
    xx <- strsplit(i, split="||", fixed=T)[[1]][6]
    xx <- strsplit(xx, split=" ")[[1]][2]
  })
)

#Get the entrezgene_id for the gene Symbols. This will be later merged with the geneIDs obtained from the Kosinski data (Kosinski_CryptTop.geneID). This is done because some gene symbols do not have associated gene IDs, and others have gene IDs but not asocaitd gene symbols, in the Kosinski data. 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
Kosinski_CryptTop.geneSymbolTogeneID <- getBM(attributes = c('entrezgene_id'),
                                              filters = 'hgnc_symbol',
                                              values = Kosinski_CryptTop.geneSymbol, 
                                              mart = human)
Kosinski_CryptTop.geneSymbolTogeneID <- as.character(Kosinski_CryptTop.geneSymbolTogeneID$entrezgene_id)


Kosinski_CryptTop.conmbined_geneID <- unique(c(Kosinski_CryptTop.geneID, Kosinski_CryptTop.geneSymbolTogeneID))

length(Kosinski_CryptTop.conmbined_geneID)


########################
# Create gmt file
########################
Kosinski_CryptTopAndBottom_GeneSetList <- list(Kosinski_CryptBottom.signature=Kosinski_CryptBottom.conmbined_geneID, Kosinski_CryptTop.signature=Kosinski_CryptTop.conmbined_geneID)

sapply(Kosinski_CryptTopAndBottom_GeneSetList, length)

# Convert the list to a GeneSetCollection using GSEABase functionality
myGeneSetList <- list()
for (i in 1:length(Kosinski_CryptTopAndBottom_GeneSetList)) {
  myGeneSetList[[i]] <- GeneSet(Kosinski_CryptTopAndBottom_GeneSetList[[i]], setName = names(Kosinski_CryptTopAndBottom_GeneSetList)[i])
}

Kosinski_CryptTopAndBottom_GeneSetList <- GeneSetCollection(myGeneSetList)

# Convert GeneSetCollection to gmt format
toGmt(Kosinski_CryptTopAndBottom_GeneSetList, con = "D:/HariEaswaran/R_Working_Directory_onc-cbio/Required_Files/GSEA/Curated_Colon_Signatures/Kosinski_CryptTop_and_CryptBase_Signature/Kosinski_CryptTopAndBottom_GeneSetList.gmt")
