library("devtools")
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")
#BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
######transcript data downlaod or hg38
#hg 38 RNA se
query.exp.hg38 <- GDCquery(project =CancerName, 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
#GDCdownload(query.exp.hg38,files.per.chunk = 1)

rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
#