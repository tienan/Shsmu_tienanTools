#################################DL dif gene###########################################
getwd()
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")

#diff_H1975 = genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]
#diff_A549 = genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]

#Union
diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.05,]$gene)
#intersect
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)

# extract the foldchange of diff exp gene
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$q_value,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$q_value
)

# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,5])>0,]

# gene name order
gene_name = as.data.frame(sort(tolower(diff_gene)))
colnames(diff_gene_filer_1 )=c("gene_name","gene_id","H1975_foldChange","p-value","A549_foldChange","p-value")
diff_gene_filer_1_union = diff_gene_filer_1

write.csv(x = diff_gene_filer_1_union,file = "diff_gene_filer_1_union_DL.csv")



#################################tumor vs benign dif gene########################################
getwd()
# query <- GDCquery(project = "TCGA-LUAD", 
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   experimental.strategy = "RNA-Seq",
#                   platform = "Illumina HiSeq",
#                   file.type = "results",
#                   legacy = TRUE)
#LUADMatrix <- assay(LUADRnaseqSE ,"raw_count")
# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2  hg19
query <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
#############################################up hg38
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LUADRnaseqSE <- GDCprepare(query)
rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
LUADMatrix <- assay(LUADRnaseqSE)
library(SummarizedExperiment)
 # or LUADMatrix <- assay(BRCARnaseqSE,"raw_count")
dataNorm <- TCGAanalyze_Normalization(tabDF = LUADMatrix, geneInfo =  geneInfo)
dataNorm_DL = LUADMatrix[rownames(LUADMatrix)%in%diff_gene_filer_1[,2],]
dataNorm_DL_sort = dataNorm_DL[order(rownames(dataNorm_DL)),]


tabSurvKM<-
  TCGAanalyze_SurvivalKM(clinical_patient_Cancer,dataNorm_DL_sort,
                         Genelist = rownames(dataNorm_DL_sort),
                         p.cut = 0.05,
                         Survresult = T,
                         ThreshTop=0.5,
                         ThreshDown=0.5)

write.csv(tabSurvKM,"sur_DL_union.csv")

?TCGAanalyze_SurvivalKM
dataNorm[rownames(dataNorm) == "CTSV"]

write.csv(rownames(dataNorm),"tmp.txt")


diff_gene_filer_1

dataNorm_DL_sort[,1]

fivenum = apply(dataNorm_DL, 2, fivenum)

gene_name_exp_carcer_sign = dataNorm_DL_sort

# nrow(DL_sign)
# nrow(dataNorm_DL)
# rownames(dataNorm_DL_sort)

for (i in 1:nrow(gene_name_exp_carcer_sign)){
  
  if (DL_sign[i]==1)# DL increase
  {
    gene_name_exp_carcer_sign[i,] = 
      dataNorm_DL_sort[i,] - 
      fivenum( dataNorm_DL_sort[i,])[4] 
    gene_name_exp_carcer_sign[i,] = 
      ifelse(gene_name_exp_carcer_sign[i,]<0,0,1)
  }else# DL decrease
  {
    gene_name_exp_carcer_sign[i,] = 
      dataNorm_DL_sort[i,] - 
      fivenum( dataNorm_DL_sort[i,])[2] 
    gene_name_exp_carcer_sign[i,] = 
      ifelse(gene_name_exp_carcer_sign[i,]<0,1,0)
  }
}


apply(gene_name_exp_carcer_sign, 2, sum)

dataNorm_DL_sort_DLSign = as.data.frame(
  rbind(dataNorm_DL_sort,apply(gene_name_exp_carcer_sign, 2, sum)))


dataNorm_DL_sort[23,] = apply(gene_name_exp_carcer_sign, 2, sum)

tabSurvKM<-
  TCGAanalyze_SurvivalKM(clinical_patient_Cancer,dataNorm_DL_sort,
                         Genelist = rownames(dataNorm_DL_sort),
                         p.cut = 0.05,
                         Survresult = T,
                         ThreshTop=0.5,
                         ThreshDown=0.5)


rownames(dataNorm_DL_sort_DLSign) = c(rownames(dataNorm_DL_sort),"DL")
###############
# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
LUADRnaseq_CorOutliers <- TCGAanalyze_Preprocessing(LUADRnaseqSE)

library(edgeR)
library(limma)
patient_id = colnames(LUADMatrix )
for (i in 1:length(patient_id)){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}


gene_name_exp_dif = LUADMatrix[,sign == "01A"|sign == "11A"]

t1 = edgeR::DGEList(gene_name_exp_dif,group = as.factor(sign[sign == "01A"|sign == "11A"]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
t3
?edgeR::exactTest
?TCGAanalyze_Normalization
dataNorm <- TCGAanalyze_Normalization(tabDF = LUADMatrix, geneInfo =  geneInfo)
# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
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
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

###############################################gene sur TCGA##################################
library(TCGAbiolinks)
diff_gene_filer_1
# Survival Analysis SA
tabSurvKMcomplete <- NULL
clinical_patient_Cancer <- GDCquery_clinic("TCGA-LUAD","clinical")

LUAD_rnaseqv2 <- LUADMatrix
dataNorm <- TCGAanalyze_Normalization(tabDF = LUAD_rnaseqv2, geneInfo =  geneInfo)
dataLUADcomplete <- dataNorm

tokenStop<- 1

tabSurvKMcomplete <- NULL

#for( i in 1: round(nrow(dataLUADcomplete)/100)){
#  tokenStart <- tokenStop
#  tokenStop <-100*i
  tabSurvKM<-
    TCGAanalyze_SurvivalKM(clinical_patient_Cancer,dataLUADcomplete,
                                    Genelist = rownames(dataLUADcomplete),
                                    Survresult = F,
                                    ThreshTop=0.5,
                                    ThreshDown=0.5)
#  tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
#}
abSurvKMcomplete = tabSurvKM
tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.05,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]
  
write.csv(tabSurvKMcomplete,file = "LUADsur.csv")



tmp_id = rownames(tabSurvKM) 
tmp=NULL
tmp_name = NULL
for (i in 1:length(tmp_id)){
  tmp = strsplit2(tmp_id[i],split = "\\|")
  tmp_name[i] = tmp[1] 
}

rownames(tabSurvKM) = tmp_name 


# the intersection with dif gene 
tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
  rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
  ]
write.csv(tabSurvKMcompleteDEGs,file = "tabSurvKMcompleteDEGsLUAD.csv")



################################### DL Cancer ######################################