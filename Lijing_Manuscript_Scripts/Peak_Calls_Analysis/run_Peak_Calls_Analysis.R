rm(list = ls())


library(rmarkdown)
#BiocManager::install("BiocStyle")
library(BiocStyle)


render("..Peak_Calls_Analysis/Peak_Calls_Analysis.Rmd", output_dir = outputDir, output_file = "Peak_Calls_Analysis.html")
