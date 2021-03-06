getwd()
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
gene_name = as.data.frame(sort((diff_gene)))
colnames(diff_gene_filer_1 )=c("gene_name","gene_id","H1975_foldChange","p-value","A549_foldChange","p-value")
diff_gene_filer_1_intersection = diff_gene_filer_1
gene_name = as.data.frame(sort((diff_gene_filer_1[,1])))
colnames(gene_name)="gene_name"
write.csv(x = diff_gene_filer_1_intersection ,file = "diff_gene_filer_1_intersection_DL.csv")

#################################################diff analysis###########################################
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
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
                                  qnt.cut =  0.01)

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
library()

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
# intersect
setwd("D:/R/DL/")
setwd("/home/tienan/R/DL/")

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
# import exp of each samples
file_list = dir(pattern = "*.fpkm*")
dir()
#extract the data of target gene
?data.frame
tmp_file=c()
tmp_file = data.frame(1:23)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp$gene_id=dat_tmp$gene_id
  merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}
tmp_file$gene_id = merge_tmp$gene_id
colnames(tmp_file)=c("id",file_list,"gene_id")
# data of TCGA 
library("limma")
####TCGA RNA data 
#install.packages("devtools")
library("devtools")
#devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")
<<<<<<< HEAD

=======
>>>>>>> 787dbc1d85f552f3290f6a41c701f3eb9306ba3b
#BiocManager::install("TCGAbiolinks")
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
#rownames(exp.hg38.values) = tolower(rownames(exp.hg38.values))
# gene collection
exp.hg38.values_targeted_gene = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
# patient_id tidy
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(exp.hg38.values_targeted_gene) = patient_id

#gene_name_exp = exp.hg38.values[rownames(exp.hg38.values)%in%gene_name$gene_name,]
#rownames(gene_name_exp)
gene_name_exp_carcer = exp.hg38.values_targeted_gene[,sign == "01A"] # grepl is also can use
head(gene_name_exp_carcer)
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
#table(gene_name_exp_carcer_sign_sum)
gene_name_exp_carcer_sign_sum  = rbind(gene_name_exp_carcer_sign,gene_name_exp_carcer_sign_sum )

###########################################DL heatmap
DL_gene_exp = rbind(gene_name_exp_carcer_sort,ifelse(gene_name_exp_carcer_sign_sum>5,100,0))

DL_gene_exp_orderDLS = DL_gene_exp[,order(DL_gene_exp[24,])]

DL_gene_exp_orderDLS[c(1:23),] = apply(DL_gene_exp_orderDLS[c(1:23),],1,scale)


?pheatmap
install.packages("pheatmap")
library(pheatmap)
tiff(filename = "Figure-2.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400, family = "", restoreConsole = TRUE
)
pheatmap(DL_gene_exp_orderDLS,
         clustering_distance_cols  = "euclidean",
         show_colnames =   F,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F
)

pheatmap(DL_gene_exp_orderDLS[c(1,24),],
         clustering_distance_cols  = "euclidean",
         scale="row",
         
         show_colnames =   F,
         cluster_cols = F,
         cluster_rows = F
)


####

names = colnames(gene_name_exp_carcer_sign_sum)
# name = c()
# for (i in 1:length(gene_name_exp_carcer_sign_sum)){
#   tmp = unlist(strsplit(names[i],split = "_"))
#   name[i]=paste(tmp[1],tmp[2],tmp[3],sep = "-")
# }
DL_statu = as.data.frame(cbind(names,t(gene_name_exp_carcer_sign_sum)))

clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)
#clinical_LUAD$
#临床数据整理
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

#head(clinical_LUAD_m1)
#DL statu & clinical data merge 
DL_statu$names
clinical_LUAD_m1$clinical_LUAD.submitter_id
# DL_statu$name=rownames(DL_statu)
clin_DL = merge(DL_statu,clinical_LUAD_m1,by.x = "names",by.y="clinical_LUAD.submitter_id")
head(clin_DL)
#cbind(clin_DL$DL_level,clin_DL$survial_day,clin_DL$clinical_LUAD.tumor_stage)

plot(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum),clin_DL$survial_day)
##########cancer stage modify 
stage = as.character(clin_DL$clinical_LUAD.tumor_stage)
stage_simple = c()
?grepl
for(i in 1:length(stage)){
  if(base::grepl('iii|iv', stage[i]))
  {stage_simple[i]=2}
  else
  {stage_simple[i]=1}
}
clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>8,1,0)

clin_DL$gene_name_exp_carcer_sign_sum=as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)

boxplot(clin_DL$survial_day~clin_DL$group)

t.test(clin_DL$survial_day~clin_DL$group)

library(survival)
library(ggplot2)
require("survival")
library(survminer)

clin_DL_stage_1 = clin_DL[stage_simple==1,]
clin_DL_stage_2 = clin_DL[stage_simple==2,]


############
#analyse the genes HR 
############
head(clin_DL)
fit <- summary(coxph(Surv(survial_day, survial_state)~ genes[,i],data=clin_DL_stage_1 ) )



head(clin_DL)
genes = clin_DL_stage_1[,c(2:25)]
tmpHR=c()
tmpP=c()
i=24
for (i in 1:ncol(genes)){
  fit <- summary(coxph(Surv(survial_day, survial_state)~ genes[,i],data=clin_DL_stage_1 ) )
  tmpHR[i]=fit$coefficients[1,2]
  tmpP[i]=fit$coefficients[1,5]
}
tmpHR
tmpP
resEval_stage_1 = cbind(colnames(genes),tmpHR,tmpP)


head(clin_DL)
genes = clin_DL_stage_2[,c(2:25)]
tmpHR=c()
tmpP=c()
i=  
for (i in 1:ncol(genes)){
  fit <- summary(coxph(Surv(survial_day, survial_state)~ genes[,i],data=clin_DL_stage_2 ) )
  tmpHR[i]=fit$coefficients[1,2]
  tmpP[i]=fit$coefficients[1,5]
}
tmpHR
tmpP
resEval_stage_2 = cbind(colnames(genes),tmpHR,tmpP)

panelName = intersect(resEval_stage_1[resEval[,2]<1,1],resEval_stage_1[resEval[,2]<1,1])

tmp = as.data.frame(clin_DL[,colnames(clin_DL)%in%panelName])

tmp=tmp[,-8]
mean(base::unlist(tmp))
apply(tmp,1,sum)
HR_nagetive_sum = c()
for (i in 1:nrow(tmp))
  HR_nagetive_sum[i] = as.numeric(as.matrix(tmp)[i,]) %>% sum
clin_DL$HR_nagetive_sum=HR_nagetive_sum


?survdiff
??survival
nrow(clin_DL)
survival::survdiff(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL)
fit <- coxph(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL) 
summary(fit)




fit <- coxph(Surv(survial_day, survial_state)~
               as.numeric(gene_name_exp_carcer_sign_sum),data=clin_DL) 
summary(fit)

clin_DL$gene_name_exp_carcer_sign_sum = as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)

fit <- coxph(Surv(survial_day, survial_state)~gene_name_exp_carcer_sign_sum,data=clin_DL) 

fit <- coxph(Surv(survial_day, survial_state)~HR_nagetive_sum,data=clin_DL) 
summary(fit)

fit <- coxph(Surv(survial_day, survial_state)~gene_name_exp_carcer_sign_sum+stage_simple,data=clin_DL) 


clin_DL_stage_1 = clin_DL[stage_simple==1,]

fit <- coxph(Surv(survial_day, survial_state)~gene_name_exp_carcer_sign_sum,data=clin_DL_stage_1) 
fit <- coxph(Surv(survial_day, survial_state)~HR_nagetive_sum,data=clin_DL_stage_1) 

summary(fit)


clin_DL_stage_2 = clin_DL[stage_simple==2,]
fit <- coxph(Surv(survial_day, survial_state)~gene_name_exp_carcer_sign_sum,data=clin_DL_stage_2) 
fit <- coxph(Surv(survial_day, survial_state)~HR_nagetive_sum,data=clin_DL_stage_2) 

summary(fit)

# 
#                            coef exp(coef) se(coef)      z Pr(>|z|)  
#gene_name_exp_carcer_sign_sum -0.03164   0.96886  0.01585 -1.995    0.046 *
#
summary(fit)

resP=c()
j=1
for (i in 4:16){
  clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>j,1,0)
  fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL)
  summary(fit)
  resP[j]= fit$chisq
  j=j+1
}






resP  # 8 is the max value 
fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL)
fit
summary(fit)

clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>8,1,0)

fit <- coxph(Surv(survial_day, survial_state)~group,data=clin_DL) 
summary(fit)


tiff(filename = "DL_sur.tiff",
     width = 2480, height = 1880, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 300)
fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL[clin_DL$survial_day<5*365,])
fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL[clin_DL$survial_day<3*365,])
fit
splots <- list()
splots[[1]] = ggsurvplot(fit, data = clin_DL,
           font.main = c(16, "bold", "darkblue"),
           font.x = c(20, "bold.italic", "darkred"),
           font.y = c(20, "bold.italic", "darkred"),
           font.tickslab = c(20, "plain", "Black"))
splots[[2]] <- ggsurvplot(fit, data = clin_DL, risk.table = TRUE, ggtheme = theme_grey())# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
dev.off()

fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL[clin_DL$survial_day<3*365,])
fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL[clin_DL$survial_day<3*365,])
fit
splots <- list()
splots[[1]] = ggsurvplot(fit, data = clin_DL,
                         font.main = c(16, "bold", "darkblue"),
                         font.x = c(20, "bold.italic", "darkred"),
                         font.y = c(20, "bold.italic", "darkred"),
                         font.tickslab = c(20, "plain", "Black"))
splots[[2]] <- ggsurvplot(fit, data = clin_DL, risk.table = TRUE, ggtheme = theme_grey())# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)

fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL[clin_DL$survial_day<1*365,])
fit = survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL[clin_DL$survial_day<1*365,])
fit
splots <- list()
splots[[1]] = ggsurvplot(fit, data = clin_DL,
                         font.main = c(16, "bold", "darkblue"),
                         font.x = c(20, "bold.italic", "darkred"),
                         font.y = c(20, "bold.italic", "darkred"),
                         font.tickslab = c(20, "plain", "Black"))
splots[[2]] <- ggsurvplot(fit, data = clin_DL, risk.table = TRUE, ggtheme = theme_grey())# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)




survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL)

ggsurvplot(fit, data = clin_DL)
?ggsurvplot

fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL_stage_1)
survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL_stage_1)
ggsurvplot(fit, data = clin_DL)
summary(fit)
??ggsurvplot

clin_DL_stage_1 = clin_DL[stage_simple==1,]
fit <- coxph(Surv(survial_day, survial_state)~group,data=clin_DL_stage_1) 
summary(fit)
fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL_stage_1)
ggsurvplot(fit, data = clin_DL_stage_1)
summary(fit)
survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL_stage_1)
ggsurvplot(fit, data = clin_DL)


clin_DL_stage_2 = clin_DL[stage_simple==2,]
fit <- coxph(Surv(survial_day, survial_state)~group,data=clin_DL_stage_2) 
summary(fit)
fit<- survfit(Surv(survial_day, survial_state)~group, data=clin_DL_stage_2)
ggsurvplot(fit, data = clin_DL_stage_2)
summary(fit)
survival::survdiff(Surv(survial_day, survial_state)~group,data=clin_DL_stage_2)
ggsurvplot(fit, data = clin_DL)







####################################################################
file_list = dir(pattern = "*.fpkm*")

#extract the data of target gene
?data.frame
tmp_file = NULL
dat_tmp = read.table(file_list[i],header = T,sep = "\t")
tmp_file =  as.data.frame(dat_tmp$FPKM)
for (i in 2:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  tmp_file =cbind(tmp_file,as.data.frame(dat_tmp$FPKM))
}
