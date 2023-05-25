#Note:
# KEGG.db is deprecated and works with Biocondictor 3.12 (and 3.11) which is available in with R 4.0. So do "module load conda_R/4.0" to load R 4.0.

# RNA-seq analyses
rm(list=ls())

#Name these file paths and fine names as required
output_dir <- "/Analysis/Diff_Expr_Analyses_Results/"
RData_filename <- "Lijing_RNA_seq_Analyses_062922.RData"


library(rmarkdown)
library(BiocStyle)
render("RNA_seq_Analyses/RNA_seq_Analyses.Rmd", output_dir = output_dir, output_file = "RNA_seq_Analyses.html")

