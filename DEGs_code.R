#installing packages (tools)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("limma")


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
<<<<<<< HEAD
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name

#  DEGs cal 
=======
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
>>>>>>> 0053f2451b502d6024244d067a08a08e5e305535
