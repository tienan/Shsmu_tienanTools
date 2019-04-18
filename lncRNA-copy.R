setwd("/media/tienan/00006784000048231/R/")

library("devtools")
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")
#BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
library("biomaRt")
######transcript data downlaod or hg38
#hg 38 RNA se
query.exp.hg38 <- GDCquery(project ="TCGA-LUSC", #TCGA-LUAD;TCGA-LUSC
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
LUADMatrix <- assay(LUADRnaseqSE ,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
head(LUADRnaseqSE)
#GDCdownload(query.exp.hg38,files.per.chunk = 1)

rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
rownames(exp.hg38.values )<- values(LUADRnaseqSE)$external_gene_name
nrow(exp.hg38.values)
rownames(LUADRnaseqSE)[rownames(LUADRnaseqSE)=="PDIA3P1"]
#
LUADRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(LUADRnaseqSE)

dataNorm <- TCGAanalyze_Normalization(tabDF = LUADRnaseqSE , geneInfo =  geneInfo)
rownames(dataNorm)[rownames(dataNorm)=="PDIA3P1"]
rownames(geneInfo)[rownames(geneInfo)=="PDIA3P1"]


dataFilt <- TCGAanalyze_Filtering(tabDF = exp.hg38.values ,
                                  method = "quantile", 
                                  qnt.cut =  0.01)

samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.1 ,
                            logFC.cut = 0.05,
                            method = "glmLRT")
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])
nrow(dataDEGsFiltLevel)

write.csv(file = "LUAD_DEGs.csv",x = dataDEGsFiltLevel)
dataDEGsFiltLevel$mRNA
dataDEGsFiltLevel$mRNA[dataDEGsFiltLevel$mRNA=="PDIA3P1"]


affyids=dataDEGsFiltLevel$mRNA
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)



refseqids = dataDEGsFiltLevel$mRNA
refseqids[refseqids=="PDIA3P1"]
refseqids = c("PDIA3P1","LOC146880")
refseqids = values(LUADRnaseqSE)$external_gene_name

ipro = getBM(attributes=c("refseq_mrna","hgnc_symbol","transcript_biotype",
                          "refseq_ncrna"), 
             filters="hgnc_symbol",
             values=refseqids, 
             mart=ensembl)
ipro
ipro$hgnc_symbol[ipro$hgnc_symbol=="PDIA3P1"]
#write.csv(file = "ALLLUAD_DEGs2Transcript.csv",x =ipro)

write.csv(file = "LUAD_DEGs2Transcript.csv",x =ipro)


############################################LUSC
query.exp.hg38 <- GDCquery(project ="TCGA-LUSC", #TCGA-LUAD;TCGA-LUSC
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
LUSCRnaseqSE <- GDCprepare(query.exp.hg38)

exp.hg38.values <- assay(LUSCRnaseqSE)
rownames(exp.hg38.values )<- values(LUSCRnaseqSE)$external_gene_name
nrow(exp.hg38.values)
rownames(exp.hg38.values)[rownames(exp.hg38.values)=="PDIA3P1"]

dataFilt <- TCGAanalyze_Filtering(tabDF = exp.hg38.values ,
                                  method = "quantile", 
                                  qnt.cut =  0.01)

samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.1 ,
                            logFC.cut = 0.05,
                            method = "glmLRT")
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])
nrow(dataDEGsFiltLevel)

dataDEGsFiltLevel$mRNA
dataDEGsFiltLevel$mRNA[dataDEGsFiltLevel$mRNA=="PDIA3P1"]

write.csv(file = "LUSC_DEGs.csv",x = dataDEGsFiltLevel)



affyids=dataDEGsFiltLevel$mRNA
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)



refseqids = dataDEGsFiltLevel$mRNA
refseqids[refseqids=="PDIA3P1"]
refseqids = c("PDIA3P1","LOC146880")
refseqids = values(LUADRnaseqSE)$external_gene_name

ipro = getBM(attributes=c("refseq_mrna","hgnc_symbol","transcript_biotype",
                          "refseq_ncrna"), 
             filters="hgnc_symbol",
             values=refseqids, 
             mart=ensembl)
ipro
ipro$hgnc_symbol[ipro$hgnc_symbol=="PDIA3P1"]
#write.csv(file = "ALLLUAD_DEGs2Transcript.csv",x =ipro)

write.csv(file = "LUSC_DEGs2Transcript.csv",x =ipro)


