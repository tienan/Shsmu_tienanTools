##########################PM2.5
<<<<<<< HEAD
setwd("../pm2.5/")
=======
library("limma")
setwd("D:/R/pm2.5/")
>>>>>>> ffb67f991571f503b06742b37b2530d36e21bb61

lipidMetabolism = read.table("lipidMetabolism.txt")
lipidMetabolism$V1 = tolower(lipidMetabolism$V1)
lipidMetabolism = read.table("Autophagy.txt")

lipidMetabolism = c(as.character(unlist(lipidMetabolism)),"LPIN")



<<<<<<< HEAD
=======

>>>>>>> d98a0e7d3e5bc5e067c6b5f1aadb28376eda804e
genes_1975_pm_con = read.table("1975_gene_exp.diff",header = T,sep = "\t")
genes_A549_pm_con = read.table("A549_gene_exp.diff",header = T,sep = "\t")

genes_A549_pm_con[genes_A549_pm_con$gene%in%lipidMetabolism,]

diff_H1975 = genes_1975_pm_con [genes_1975_pm_con$q_value<0.1,]
diff_A549 = genes_A549_pm_con[genes_A549_pm_con$q_value<0.1,]

nrow(diff_H1975)
nrow(diff_A549)

diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.05,]$gene)



diff_gene = lipidMetabolism
colnames(diff_gene) = "geneName"
# gene_name

file_list = dir(pattern = "*.fpkm*")
# setwd("../DL/")
# tmp_list = dir(pattern = "*.fpkm*")
# file_list=c(tmp_list[grep(pattern = "con",x=tmp_list)],file_list)

?gsub
diff_gene$geneName = tolower(diff_gene$geneName)
remove(tmp_file)
tmp_file = data.frame(1:6)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp$gene_id=tolower(dat_tmp$gene_id)
  merge_tmp = merge(dat_tmp,diff_gene,by.x = "gene_id",by.y  = "geneName")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}

ncol(tmp_file)
nrow(tmp_file)
merge_tmp$gene_id
head(tmp_file)


#data governance
#rownames(tmp_file) = 
colnames(tmp_file) = c("No",gsub("_genes.fpkm_tracking", "", file_list))

tmp_file = cbind(merge_tmp$gene_id,tmp_file)
colnames(tmp_file)[1]="geneName"
colnames(tmp_file) = c("geneName","No",gsub(".fpkm_tracking", "", file_list))
colnames(tmp_file) = c("geneName","No",gsub("_genes", "", colnames(tmp_file)[3:14]))
tmp_file_1 = tmp_file[-1,]
rownames(tmp_file_1 ) = tmp_file_1$geneName
colnames(tmp_file_1)

write.csv(tmp_file,"autophagy_pm2.5.txt")


tmp_file_1 = na.omit(tmp_file[test_1<0.1,])
#tmp_file[test_1<0.1|test_2<0.1,]
# tmp_file_1 = na.omit(tmp_file[test_1<0.1|test_2<0.1,])
# 
# tmp_file_1 = na.omit(tmp_file[test_1<0.1,])

colnames(tmp_file_1)



tmp_file_1[(sum(tmp_file_1[3:5])-sum(tmp_file_1[6:8]))*(sum(tmp_file_1[9:11])-sum(tmp_file_1[12:14]))>0,]

head(tmp_file_1)
colnames(tmp_file_1) = gsub(".fpkm_tracking", "", colnames(tmp_file_1))

rownames(tmp_file_1) = tmp_file_1$`merge_tmp$gene_id`

rownames(tmp_file_1)

library(pheatmap)

pheatmap(tmp_file[,c(3:14)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
pheatmap(tmp_file[,c(2:13)])

pheatmap(tmp_file[,c(2:13)],cluster_rows = F)


tiff(filename = "Figure_PM_3a.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file_1[,c(3:8)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
#         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()


?tiff
tiff(filename = "Figure_PM_4b.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file_1[,c(9:14)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()

tiff(filename = "Figure-1c.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file[,c(2:13)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()





##########################DL analysis
source("https://bioconductor.org/biocLite.R")
BiocManager::install("edgeR")
install.packages("gdata")
install.packages("heatmaply")
install.packages("ggplot2")
install.packages("genefilter")
BiocManager::install("methylumi")

BiocManager::install("vsn")
BiocManager::install("preprocessCore")
BiocManager::install("gridExtra")
BiocManager::install("EBSeq")
biocLite("blockmodeling")
BiocManager::install("biomaRt")
library()
library("EBSeq")

data(GeneMat)



# BiocManager::install("tscvh")
# BiocInstaller::biocLite("tscvh")

library(tscvh)


library("edgeR")
library("gdata")
library("heatmaply")
library("ggplot2")
library("genefilter")
library("methylumi")


####################################################DL
getwd()
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")


genes_1975_dl_con[genes_1975_dl_con$gene=="SLPI",]

genes_A549_dl_con [genes_1975_dl_con$gene=="VCAN",]

#cutoff < 0.1 to get diff exp gene data
diff_H1975 = genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]
diff_A549 = genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]

nrow(diff_H1975)
nrow(diff_A549)

# intersect of gene name 
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                      genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)
# or Union
diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05,]$gene,
                     genes_A549_dl_con[genes_A549_dl_con$q_value<0.05,]$gene)


diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05&
                                      abs(genes_1975_dl_con$log2.fold_change.)>1,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.05&
                                      abs(genes_A549_dl_con$log2.fold_change.)>1,]$gene)


diff_gene

# extract the foldchange of diff exp gene
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.)

# extract the data with the same direction 
diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,4])>0,]

# gene name order
gene_name = as.data.frame(sort(tolower(diff_gene_filer_1[,1])))
colnames(gene_name)="gene_name"

#sort(gene_name)

# import exp of each samples
file_list = dir(pattern = "*.fpkm*")

#extract the data of target gene
?data.frame
tmp_file = data.frame(1:29)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp$gene_id=tolower(dat_tmp$gene_id)
  merge_tmp = merge(dat_tmp,geene_name,by.x = "gene_id",by.y  = "gene_name")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}
#######################################
file_list = dir(pattern = "*.fpkm*")
?data.frame
tmp_file = data.frame(1)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  merge_tmp = dat_tmp[dat_tmp$gene_id=="VCAN",]
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}
colnames(tmp_file) = c("No",gsub("_genes.fpkm_tracking", "", file_list))

########################################

ncol(tmp_file)

#data governance
rownames(tmp_file) = gene_name$gene_name
colnames(tmp_file) = c("No",gsub("_genes.fpkm_tracking", "", file_list))


# plot the heatmap 
library(pheatmap)

pheatmap(tmp_file[,c(2:13)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
pheatmap(tmp_file[,c(2:13)])

pheatmap(tmp_file[,c(2:13)],cluster_rows = F)


tiff(filename = "Figure-1a.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file[,c(2:7)],
#         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()

tiff(filename = "Figure-1b.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file[,c(8:13)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()

tiff(filename = "Figure-1c.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(tmp_file[,c(2:13)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)
dev.off()



# data of TCGA 

library("limma")
?lmFit
####TCGA RNA data 
#install.packages("devtools")
library("devtools")
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")

#download clinical data 
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)

#clinical data tidy
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = data.frame(clinical_LUAD$bcr_patient_barcode,
                             clinical_LUAD$tumor_stage,
                             clinical_LUAD$days_to_last_follow_up,
                             clinical_LUAD$days_to_death)
head(clinical_LUAD_m1 )

#clinical_LUAD


#TCGAbiolinks:::getGDCprojects()$project_id
# 
# 
# query = GDCquery(project = "TCGA-LUAD", data.category = "Raw Sequencing Data", experimental.strategy = "RNA-Seq")
# 
# 
# query <- GDCquery(project = "TCGA-LUAD",data.category = "Transcriptome Profiling")
# 
# 
# 
# query.exp <- GDCquery(project= "TCGA-LUAD", 
#                       
#                       data.category = "Gene expression", 
#                       data.type = "Gene expression quantification",
#                       platform = "Illumina HiSeq", experimental.strategy = "RNA-Seq",
#                       legacy = TRUE)
# 
# GDCdownload(query, method = "api")
# expdat <- GDCprepare(query = query.exp,
#                      save = TRUE, save.filename = "exp.rda")
# 
# 
# #Question: TCGAbiolinks TCGA-BRCA RNA-seq clinical data
# #https://www.biostars.org/p/320274/
# CancerProject <- "TCGA-LUAD"
# 
# query <- GDCquery(project = CancerProject,
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification", 
#                   workflow.type = "HTSeq - Counts")
# 
# samplesDown <- getResults(query,cols=c("cases"))
# 
# dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
#                                   typesample = "TP")
# 
# queryDown <- GDCquery(project = CancerProject, 
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification", 
#                       workflow.type = "HTSeq - Counts", 
#                       barcode = dataSmTP)
# 
# GDCdownload(query = queryDown,directory = "LUAD_RESULTS/TCGA/htseq_data/")                    
# 
# dataPrep <- GDCprepare(query = queryDown, 
#                        save = TRUE, 
#                        directory =  "LUAD_RESULTS/TCGA/htseq_data/",
#                        save.filename = "htseq_counts.rda", summarizedExperiment = TRUE)
# 
# rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
# 
# check = values(exp.hg38)
# check$original_ensembl_gene_id
# #Question: TCGA-biolink: GRCH38 RNA-seq, normalized counts matrix, Converting ensemble gene name into common gene name
#https://support.bioconductor.org/p/101276/
#https://gist.githubusercontent.com/tiagochst/fab346c19fa97d62a4bfb943024d1566/raw/0e7846e7b838c8c8403e4d915fc735001d8f1d06/getExp.R
BiocManager::install("TCGAbiolinks")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
######transcript data downlaod  hg 19  or  hg38
query.exp.hg19 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Gene expression", 
                           data.type = "Gene expression quantification",
                           legacy = TRUE)
exp.hg19 <- GDCprepare(query = query.exp.hg19)
exp.values <- assay(query.exp.hg19)

#hg 38 RNA se
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
GDCdownload(query.exp.hg38,files.per.chunk = 1)

exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
colnames(exp.hg38.values)
colnames(exp.hg38)

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


###########################dif analysis
# Query platform Illumina HiSeq with a list of barcode 
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


dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA=="SLPI",]

dif_gene = tolower(dataDEGsFiltLevel$mRNA)


target_dl_mutation = toupper(c("flnc","vcan","sh3tc1","c1s","xdh","smarca1","col12a1","ptprm",
  "trank1","atad5","ept1","hrnr") )  

library(maftools)
?subsetMaf

diff_gene

target_dl_deg_mut = dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%diff_gene, ]$mRNA

target_dl_deg_mut[target_dl_deg_mut == "SLPI"]


target_dl_deg_mut = dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%target_dl_mutation, ]$mRNA

laml_LUAD <- read.maf(maf="../TCGA.LUAD.varscan.acb6852e-dd48-4ca5-80f2-3d1a2c7d7ceb.DR-10.0.somatic.maf.gz")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top = 100)
mut_genes = subsetMaf(maf = laml)

mut_genes[mut_genes$Hugo_Symbol%in%target_dl_deg_mut,]$Hugo_Symbol 

dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%mut_genes[mut_genes$Hugo_Symbol%in%target_dl_deg_mut,]$Hugo_Symbol ,]



fold_1975 = genes_1975_dl_con[genes_1975_dl_con$gene%in%mut_genes[mut_genes$Hugo_Symbol%in%target_dl_deg_mut,]$Hugo_Symbol,]$log2.fold_change.
fold_A549 = genes_A549_dl_con [genes_1975_dl_con$gene%in%mut_genes[mut_genes$Hugo_Symbol%in%target_dl_deg_mut,]$Hugo_Symbol,]$log2.fold_change.

genes_1975_dl_con[genes_1975_dl_con$gene=="SLPI",]
genes_A549_dl_con[genes_1975_dl_con$gene=="SLPI",]

fold_1975*fold_A549>0

H1975 = genes_1975_dl_con[genes_1975_dl_con$gene%in%mut_genes[mut_genes$Hugo_Symbol%in%target_dl_deg_mut,]$Hugo_Symbol,]

H1975 = H1975[fold_1975*fold_A549>0,]
#H1975$gene_id=="SLPI"

H1975[H1975$log2.fold_change.*
  dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%H1975$gene,]$logFC<0,]$gene



                                                                    
mut_genes[mut_genes$Hugo_Symbol%in%H1975[H1975[fold_1975*fold_A549>0,]$log2.fold_change.*
                                           dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%H1975$gene,]$logFC<0,]$gene,]


<<<<<<< HEAD
<<<<<<< HEAD
laml_LUAD
=======
mut_genes[mut_genes$Hugo_Symbol%in%H1975[H1975[fold_1975*fold_A549>0,]$log2.fold_change.*
                                           dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%H1975$gene,]$logFC<0,]$gene,]
>>>>>>> d98a0e7d3e5bc5e067c6b5f1aadb28376eda804e
=======
mut_genes[mut_genes$Hugo_Symbol%in%H1975[H1975[fold_1975*fold_A549>0,]$log2.fold_change.*
                                           dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA%in%H1975$gene,]$logFC<0,]$gene,]
>>>>>>> d98a0e7d3e5bc5e067c6b5f1aadb28376eda804e


sign =c()

?ebayes
########## simulation & illumination
# Simulate gene expression data for 100 probes and 6 microarrays
# Microarray are in two groups
# First two probes are differentially expressed in second group
# Std deviations vary between genes with prior df=4
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,1))
options(digits=3)

# Ordinary fit
fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit,coef=2)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)


?lmFit
set.seed(2016)
sigma2 <- 0.05 / rchisq(100, df=10) * 10
y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
y[1,4:6] <- y[1,4:6] + 1
fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit,coef=2)


length(y)
y


#  Moderated t-statistic
fit <- eBayes(fit)
topTable(fit,coef=2)
#########


colnames(exp.hg38.values) = patient_id



exp.hg38.values_targeted_gene_dif = exp.hg38.values_targeted_gene

colnames(exp.hg38.values_targeted_gene_dif ) =  

###########################

#
min(gene_name_exp_carcer)

#ncol(gene_name_exp_carcer)

#??strsplit2

#tmp_file

#DL statue simulation 
DL_state = apply(tmp_file[,c(5:7,11:13)],1,mean)
normal_state = apply(tmp_file[,c(2:4,8:10)],1,mean)
diff_DL_nor = DL_state - normal_state 

#1. DL increase; 0. DL decease

DL_sign = c()
DL_sign = ifelse(diff_DL_nor<0,0,1)
DL_sign_sort = DL_sign[order(names(DL_sign ))]

gene_name_exp_carcer_sort = 
  gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),]
rownames(gene_name_exp_carcer_sort)
names(DL_sign_sort)
?ifelse

# DL statue simulation calculation  >75%  <25%  model 1 
gene_name_exp_carcer_sign=gene_name_exp_carcer_sort


for (i in 1:length(DL_sign_sort)){

  if (DL_sign[i]==1)# DL increase
    {
    gene_name_exp_carcer_sign[i,] = 
      gene_name_exp_carcer_sort[i,] - 
      fivenum(gene_name_exp_carcer_sort[i,])[4] 
    gene_name_exp_carcer_sign[i,] = 
      ifelse(gene_name_exp_carcer_sign[i,]<0,0,1)
  }else# DL decrease
  {
    gene_name_exp_carcer_sign[i,] = 
      gene_name_exp_carcer_sort[i,] - 
      fivenum(gene_name_exp_carcer_sort[i,])[2] 
    gene_name_exp_carcer_sign[i,] = ifelse(gene_name_exp_carcer_sign[i,]<0,1,0)
    }
#  gene_name_exp_carcer_sign[,i] 
}
gene_name_exp_carcer_sign_sum = apply( gene_name_exp_carcer_sign,2,sum)

table(gene_name_exp_carcer_sign_sum)


#gene_name_exp_carcer_sign

#ncol(gene_name_exp_carcer_sign)

# DL_sign_res=c()
# for (i in 1:ncol(gene_name_exp_carcer)){
#   DL_sign_res[i] = sum(abs(gene_name_exp_carcer_sign[,i] - DL_sign)) 
# }

#table(DL_sign_res)


# name tidy 
names = names(gene_name_exp_carcer_sign_sum)
name = c()
for (i in 1:length(gene_name_exp_carcer_sign_sum)){
  tmp = unlist(strsplit(names[i],split = "-"))
  name[i]=paste(tmp[1],tmp[2],tmp[3],sep = "-")
}
DL_statu = as.data.frame(cbind(name,gene_name_exp_carcer_sign_sum))

#
#DL_statu =as.data.frame(bind(colnames(gene_name_exp_carcer),gene_name_exp_carcer_sign))
# DL_statu = rbind(,)
# colnames(DL_statu) = c("bar_code","DL_level")
# DL_statu$bar_code=as.character(DL_statu$bar_code)
# DL_statu$bar_code <- gsub("_", "-", DL_statu$bar_code)
#head(DL_statu)
# head(gene_name_exp_carcer_sort[,c(1:5)])  head(gene_name_exp_carcer_sign)


# pheatmap(cbind(DL_sign,gene_name_exp_carcer_sign),
#          clustering_distance_cols  = "euclidean",
# #         scale="row",
#          show_colnames =   F,cluster_rows = F,cluster_cols=F)
# 
# ?pheatmap
# 
# library(pheatmap)
# tiff(filename = "Figure-2.tif",
#      width = 1800, height = 3000, units = "px", pointsize = 12,
#      compression = "lzw",
#      bg = "white", res = 400, family = "", restoreConsole = TRUE
# )
# pheatmap(apply(gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),],2,powerT),
#          clustering_distance_cols  = "euclidean",
#          scale="row",
#          show_colnames =   F,cluster_rows = F
# )
# dev.off()
# 
# 
# powerT = function(a){
#   return(a^(0.9))
# }
# 
# powerT(c(5,1))
# 
# apply(gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),],2,log)[3,]
# 
# ((gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),])[3,])^(0.01)
# 
# 
# pheatmap
# 
# 
# (diff_DL_nor + normal_state)/2
# 
# (gene_name_exp - (diff_DL_nor + normal_state)/2)[c(1:10),c(1:10)]
# 
# 
# 
# gene_name_exp[,1] - apply(gene_name_exp,1,median)
# 
# gene_name_exp = gene_name_exp[order(rownames(gene_name_exp)),]
# 
# ?pheatmap
# 
# 
# gene_name_exp = gene_name_exp[order(rownames(gene_name_exp)),]
#临床数据下载
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)

#clinical_LUAD$
#临床数据整理
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = DataFrame(clinical_LUAD$submitter_id,
                             clinical_LUAD$tumor_stage,
                             clinical_LUAD$days_to_last_follow_up,
                             clinical_LUAD$days_to_death)

survial_day=c()
for (i in 1:nrow(clinical_LUAD_m1)){
  survial_day[i] = ifelse(is.na(clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i]),clinical_LUAD$days_to_death[i],clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i])
}
survial_state = c()
for (i in 1:nrow(clinical_LUAD_m1)){
  survial_state[i] = ifelse(is.na(clinical_LUAD_m1$clinical_LUAD.days_to_last_follow_up[i]),1,0)
}

clinical_LUAD_m1$survial_day = survial_day
clinical_LUAD_m1$survial_state = survial_state

#head(clinical_LUAD_m1)
#DL statu & clinical data merge 
DL_statu$name
clinical_LUAD_m1$clinical_LUAD.bcr_patient_barcode
clin_DL = merge(DL_statu,clinical_LUAD_m1,by.x = "name",by.y="clinical_LUAD.bcr_patient_barcode")
#head(clin_DL)
#cbind(clin_DL$DL_level,clin_DL$survial_day,clin_DL$clinical_LUAD.tumor_stage)

plot(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum),clin_DL$survial_day)

stage = as.character(clin_DL$clinical_LUAD.tumor_stage)
stage_simple = c()

?grepl
for(i in 1:length(stage)){
  if(grepl('[i,ii][a,b]', stage[i]))
    {stage_simple[i]=1}
  else
    {stage_simple[i]=2}
}

         
    

clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>10,1,0)

boxplot(clin_DL$survial_day~clin_DL$group)

t.test(clin_DL$survial_day~clin_DL$group)

library(survival)
?survdiff
survival::survdiff(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL)



?ks.test()

?boxplot

# 
# plot(as.numeric(clin_DL$DL_level),clin_DL$survial_day)
# 
# DL_statu$bar_code
# clinical_LUAD_m1$clinical_LUAD.submitter_id
# 
# sink("tmp")
# sink()
# 
# ?sink
# 
# ?merge
# 
# ?ifelse
# 
# 
# ?strsplit2
# 

#http://www.bio-info-trainee.com/1980.html heatmap





# query <- GDCquery(project = "TCGA-GBM",
#                   legacy = TRUE,
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq",
#                   file.type = "normalized_results",
#                   experimental.strategy = "RNA-Seq")
# GDCdownload(query, method = "api")
# 
# data_1 <- GDCprepare(query,add.gistic2.mut = c("PTEN","FOXJ1"))



clin.gbm <- GDCquery_clinic("TCGA-LUAD", "clinical")
TCGAanalyze_survival(clinical_LUAD ,
                     "gender",
                     main = "TCGA Set\n GBM",height = 10, width=10)
