#installing packages (tools)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("limma")
BiocManager::install("GenomeInfoDb")
BiocManager::install("SummarizedExperiment", version = "3.5")
install.packages("SummarizedExperiment")

#importing the downloaded packages(tools) 
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)


#Downloading GE (gene expression data)

#hg 38 RNA se
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
GDCdownload(query.exp.hg38,files.per.chunk = 1)

exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name

#  DEGs cal 
