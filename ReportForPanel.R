getwd()
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
gene_name = as.data.frame(sort(tolower(diff_gene_filer_1[,1])))
colnames(gene_name)="gene_name"


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

rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
head(exp.hg38.values)
#write.csv(exp.hg38.values,file = "stad_exp_hg38_FPKM.csv")
# extract the targeted gene
rownames(exp.hg38.values) = tolower(rownames(exp.hg38.values))
# gene collection
exp.hg38.values_targeted_gene = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
# patient_id tidy
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
  patient_id[i] = tmp
}
colnames(exp.hg38.values_targeted_gene) = patient_id

#gene_name_exp = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
#rownames(gene_name_exp)
gene_name_exp_carcer = exp.hg38.values_targeted_gene[,sign == "01A"]
##################DL condition 
DL_state = apply(tmp_file[,c(5:7,11:13)],1,mean)
normal_state = apply(tmp_file[,c(2:4,8:10)],1,mean)
diff_DL_nor =as.data.frame(DL_state - normal_state)
rownames(diff_DL_nor)=tmp_file$gene_id

#1. DL increase; 0. DL decease

DL_sign = c()
DL_sign = ifelse(diff_DL_nor<0,0,1)
DL_sign_sort = DL_sign[order(rownames(DL_sign))]
gene_name_exp_carcer_sort = 
  gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),]
rownames(gene_name_exp_carcer_sort)

#?ifelse
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
gene_name_exp_carcer_sign_sum = apply(gene_name_exp_carcer_sign,2,sum)
table(gene_name_exp_carcer_sign_sum)

names = names(gene_name_exp_carcer_sign_sum)
name = c()
for (i in 1:length(gene_name_exp_carcer_sign_sum)){
  tmp = unlist(strsplit(names[i],split = "_"))
  name[i]=paste(tmp[1],tmp[2],tmp[3],sep = "-")
}
DL_statu = as.data.frame(cbind(name,gene_name_exp_carcer_sign_sum))

clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)


linical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = data.frame(clinical_LUAD$submitter_id,
                              clinical_LUAD$tumor_stage,
                              clinical_LUAD$days_to_last_follow_up,
                              clinical_LUAD$age_at_diagnosis,
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

clin_DL = merge(DL_statu,clinical_LUAD_m1,by.x = "name",by.y="clinical_LUAD.submitter_id")
head(clin_DL)
#cbind(clin_DL$DL_level,clin_DL$survial_day,clin_DL$clinical_LUAD.tumor_stage)

plot(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum),clin_DL$survial_day)



stage = as.character(clin_DL$clinical_LUAD.tumor_stage)
stage_simple = c()
?grepl
for(i in 1:length(stage)){
  if(grepl('iii|iv', stage[i]))
  {stage_simple[i]=2}
  else
  {stage_simple[i]=1}
}
clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>15,1,0)

boxplot(clin_DL$survial_day~clin_DL$group)

t.test(clin_DL$survial_day~clin_DL$group)

library(survival)
library(ggplot2)
require("survival")
library(survminer)

survival::survdiff(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL)
fit <- coxph(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL) 
summary(fit)




fit <- coxph(Surv(survial_day, survial_state)~
               as.numeric(gene_name_exp_carcer_sign_sum),data=clin_DL) 
summary(fit)

##########################optmal cufoff value
resP=c()
j=1
for (i in 4:16){
  clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>j,1,0)
  fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL)
  summary(fit)
  resP[j]= fit$chisq
  j=j+1
}

clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>8,1,0)
fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL)
fit
fit<- survfit(Surv(survial_day, survial_state)~group,data=clin_DL)

ggsurvplot(fit, data = clin_DL)



#################################################other dataset
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE42127", GSEMatrix =TRUE, AnnotGPL=TRUE )


exprSet <- exprs(gset[[1]])

pData <- pData(gset[[1]])

fdata<-fData(gset[[1]])

fdata_target = fdata[fdata[,3]%in%diff_gene_filer_1[,1],c(1,3)]

exprSet 



sample <- pData$geo_accession
group = ifelse(grepl(pattern = "non",pData[,1]),1,0)
group <- rep(c(1,0),times=c(12,11))
design <- model.matrix(~group)
rownames(design)=colnames(exprSet)














































#/media/tienan/00006784000048231/R/Shsmu_tienanTools
getwd()
setwd("../Shsmu_tienanTools/")
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")

# extract the foldchange of diff exp gene
diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                      genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)

# extract the data with the same direction 
diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                       as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$q_value,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                       genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$q_value
)

diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,5])>0,]

# gene name order
gene_name = as.data.frame(sort(diff_gene))
colnames(diff_gene_filer_1 )=c("gene_name","gene_id","H1975_foldChange","p-value","A549_foldChange","p-value")
diff_gene_filer_1_intersection = diff_gene_filer_1
#write.csv(x = diff_gene_filer_1_intersection ,file = "diff_gene_filer_1_intersection_DL.csv")

#TCGA validation
#installation readme is in the DL_analysis.R
library("TCGAbiolinks")
#hg19
query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  legacy = TRUE)

GDCdownload(query)



################################################hg38
#hg 38 RNA se
query <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")


LUADRnaseqSE <- GDCprepare(query)


################################################

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
LUADRnaseqSE <- GDCprepare(query)
GDCdownload(query)

LUADMatrix <- assay(LUADRnaseqSE) 

rownames(LUADMatrix) <- values(LUADRnaseqSE)$external_gene_name

dataNorm <-LUADMatrix

####################################

dataNorm <- TCGAanalyze_Normalization(tabDF = LUADMatrix, geneInfo =  geneInfo) #hg19 is requried this step


targetedGeneExp = dataNorm[rownames(dataNorm)%in%diff_gene_filer_1_intersection[,1],]

targetedGeneExp = targetedGeneExp[order(rownames(targetedGeneExp)),]



targetedGeneExp  = targetedGeneExp [,grepl("01A",colnames(targetedGeneExp))]


#######Patients ID##############################
patient_id = colnames(targetedGeneExp)
for (i in 1:length(colnames(targetedGeneExp))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(targetedGeneExp) = patient_id

##calculating DL condition 

DL_gene_up = diff_gene_filer_1_intersection[diff_gene_filer_1_intersection[,3]>0,1]
DL_gene_down = diff_gene_filer_1_intersection[diff_gene_filer_1_intersection[,3]<0,1]

#calculating median value 

targetedGeneExp[rownames(targetedGeneExp)%in%DL_gene_down,] = -targetedGeneExp[rownames(targetedGeneExp)%in%DL_gene_down,]
quantile = apply(targetedGeneExp,1,fivenum)
quantile = quantile[,order(colnames(quantile))]


targetedGeneExpDL = targetedGeneExp
for (i in 1:ncol(targetedGeneExp)){
  targetedGeneExpDL[,i] = targetedGeneExp[,i] - quantile[3,]
}

targetedGeneExpDL[targetedGeneExpDL>0]=1
targetedGeneExpDL[targetedGeneExpDL<0]=0
PatientDLCon = apply(targetedGeneExpDL, 2, sum)
PatientDLCon[order(names(PatientDLCon))]
length(PatientDLCon)
name = names(PatientDLCon) 
PatientDLCon = as.data.frame(PatientDLCon)
PatientDLCon$ID = name 

######################################Clincal Data

clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = data.frame(clinical_LUAD$submitter_id,
                              clinical_LUAD$tumor_stage,
                              clinical_LUAD$days_to_last_follow_up,
                              clinical_LUAD$age_at_diagnosis,
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
clinical_LUAD_m1= clinical_LUAD_m1[,c(1,2,6,7)]


colnames(clinical_LUAD_m1) = clinical_names



nrow(clinical_LUAD_m1)
head(clinical_LUAD_m1)

GeneExpCli = merge(clinical_LUAD_m1, PatientDLCon, by.x = "Tumor_Sample_Barcode", by.y = "ID")



stage_simple = c()
?grepl
for(i in 1:length(GeneExpCli$FAB_classification)){
  if(grepl('iii|iv', GeneExpCli$FAB_classification[i]))
  {stage_simple[i]=2}
  else
  {stage_simple[i]=1}
}
GeneExpCli$stage_simple = stage_simple
head(GeneExpCli)

GeneExpCli$group = ifelse(as.numeric(GeneExpCli$PatientDLCon)>14,1,0)



boxplot(GeneExpCli$days_to_last_followup~GeneExpCli$group )
plot(GeneExpCli$PatientDLCon,GeneExpCli$days_to_last_followup)


t.test(GeneExpCli$days_to_last_followup~GeneExpCli$group)


fit <- coxph(Surv(days_to_last_followup,  Overall_Survival_Status)~PatientDLCon+stage_simple,data=GeneExpCli) 
summary(fit)

fit <- coxph(Surv(days_to_last_followup,  Overall_Survival_Status)~PatientDLCon,data=GeneExpCli) 
summary(fit)


GeneExpCli_1 = GeneExpCli[GeneExpCli$stage_simple==1,] 

fit <- coxph(Surv(days_to_last_followup,  Overall_Survival_Status)~PatientDLCon,data=GeneExpCli_1) 
summary(fit)


GeneExpCli_2 = GeneExpCli[GeneExpCli$stage_simple==2,] 

fit <- coxph(Surv(days_to_last_followup,  Overall_Survival_Status)~PatientDLCon,data=GeneExpCli_2) 
summary(fit)












