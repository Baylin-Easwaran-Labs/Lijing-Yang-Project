# Lijing-Yang-Project
Please edit the file paths as required.

The overall pipeline and the files used for the analyses are shown in AnalysesOverview.tiff

The accession number in GEO database for the RNAs-seq data reported in this paper is GSE218480, Chip-seq data is GSE218479.

RNA-seq analyses
1) The file Salmon_alignment.sh has codes for aligning the RNA-seq data.
2) RNA_seq_Analyses.Rmd processes aligned files and summarizes  Salmon aligned transcripts to genes, normalize RNA-seq data, and performs differential expression analyses. The file RNA_seq_Analyses.R initiates RNA_seq_Analyses.Rmd run. The output of this file is Lijing_RNA_seq_Analyses_062922.Rdata where all the analyzed files are stored. The RData file can be used for various downstream analyses.
3) RNA_seq_Analyses.html has the output R Markdown output from executing RNA_seq_Analyses.Rmd.  
4) Gene set enrichment analysis of the RNA-seq data can be reproduced by running the codes in the folder FigureCodes.



ChIP-seq analyses
1) ChIP-seq reads are aligned and peaks are called using ChIPSeq_Alignment_MACSPeak.sh
2) Differential peak analyses is performed using DiffBind_Analyses.R
3) Analyses of the peak calls is done using run_Peak_Calls_Analysis.R which executes Peak_Calls_Analysis.Rmd
4) The file Peak_Calls_Analysis_062922.html contains the R Markdown output of Peak_Calls_Analysis.Rmd.
5) Peak_Calls_Analysis.RData contains all R objects from running Peak_Calls_Analysis.Rmd. This file can be used for further analyses of the ChIP-seq data.

The Folder FigureCodes contains the codes for reproducing the figures.
