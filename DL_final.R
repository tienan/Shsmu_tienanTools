setwd("../Shsmu_tienanTools/")
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")

#diff_H1975 = genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]
#diff_A549 = genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]

#Union
diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.05,]$gene)

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



# intersect of gene name 
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                      genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)# extract the foldchange of diff exp gene
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
diff_gene_filer_1_intersection = diff_gene_filer_1
write.csv(x = diff_gene_filer_1_intersection ,file = "diff_gene_filer_1_intersection_DL.csv")

#################################################diff analysis###########################################
getwd()
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)
# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LUADRnaseqSE <- GDCprepare(query)

LUADMatrix <- assay(LUADRnaseqSE ,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
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

# diff gene analysis

####mutation data##############################################################################################
laml_LUAD <- read.maf(maf="../TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf.gz")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top = 100)
mut_genes = subsetMaf(maf = laml)

#target_dl_deg_mut$mRNA
diff_gene_filer_1_intersection = as.data.frame(diff_gene_filer_1_intersection)

dl_interaction_cancer = merge(dataDEGsFiltLevel,diff_gene_filer_1_intersection,by.x="mRNA",by.y="gene_id")

dl_interaction_cancer[dl_interaction_cancer$logFC*as.numeric(as.character(dl_interaction_cancer$H1975_foldChange))<0,]


diff_gene_filer_1_union = as.data.frame(diff_gene_filer_1_union)

dl_union_cancer = merge(dataDEGsFiltLevel,diff_gene_filer_1_union,by.x="mRNA",by.y="gene_id")

dl_union_cancer[dl_union_cancer$logFC*as.numeric(as.character(dl_union_cancer$H1975_foldChange))<0,]

write.csv(x = dl_interaction_cancer[dl_interaction_cancer$logFC*as.numeric(as.character(dl_interaction_cancer$H1975_foldChange))<0,]
 ,file = "dl_interaction_cancer.csv")

write.csv(x = dl_union_cancer[dl_union_cancer$logFC*as.numeric(as.character(dl_union_cancer$H1975_foldChange))<0,]
 ,file = "dl_union_cancer.csv")

#mutation###############################################################
setwd("../DL/")
dl_union_cancer$mRNA

mut_genes = subsetMaf(maf = laml_LUAD)

mut_genes[mut_genes$Hugo_Symbol%in%dl_interaction_cancer$mRNA,c(1:20)]

write.csv(x = mut_genes[mut_genes$Hugo_Symbol%in%dl_interaction_cancer$mRNA,c(1:20)]
          ,file = "dl_interaction_cancer_mutation.csv")

mut_genes[mut_genes$Hugo_Symbol%in%dl_union_cancer$mRNA,c(1:20)]

write.csv(x = mut_genes[mut_genes$Hugo_Symbol%in%dl_union_cancer$mRNA,c(1:20)]
          ,file = "dl_union_cancer_mutation.csv")

##### DL condition ###################################################################

setwd("../DL/")

# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")
#cutoff < 0.1 to get diff exp gene data
diff_H1975 = genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]
diff_A549 = genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]
# intersect of gene name 
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                      genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)
# extract the foldchange of diff exp gene
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.)
# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,4])>0,]
# data of TCGA 

library("limma")
####TCGA RNA data 
#install.packages("devtools")
library("devtools")
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")

BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
######transcript data downlaod or hg38
#hg 38 RNA se
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
#GDCdownload(query.exp.hg38,files.per.chunk = 1)

exp.hg38.values <- assay(LUADRnaseqSE )
rownames(LUADRnaseqSE ) <- values(LUADRnaseqSE)$external_gene_name


#write.csv(exp.hg38.values,file = "stad_exp_hg38_FPKM.csv")

# extract the targeted gene
rownames(exp.hg38.values) = tolower(rownames(exp.hg38.values))


# gene collection
exp.hg38.values_targeted_gene = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]

# patient_id tidy
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "_"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}

colnames(exp.hg38.values) = patient_id

#gene_name_exp = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
#rownames(gene_name_exp)
gene_name_exp_carcer = exp.hg38.values_targeted_gene[,sign == "01A"]


