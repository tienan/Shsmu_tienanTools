# DL manuscript code
# Figure 1. Significant genetic profile by DL 


# installing/loading the package:
if(!require(installr))
  { install.packages("installr")
  require(installr)} #load / install+load installr
install.pandoc()



library("edgeR")
library("gdata")
library("heatmaply")
library("ggplot2")
library("genefilter")
library("methylumi")


############################################################
#Generating the profile of DEGs 
#############################################################
getwd()
setwd("/home/tienan/R/Shsmu_tienanTools/")
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
# import exp of each samples
file_list = dir(pattern = "*.fpkm*")
#dir()
#extract the data of target gene
?data.frame
tmp_file=c()
tmp_file = data.frame(1:23)
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp$gene_id=tolower(dat_tmp$gene_id)
  merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
  tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
}
tmp_file$gene_id = merge_tmp$gene_id

colnames(tmp_file) = c("No",base::gsub("_genes.fpkm_tracking", "", file_list),"geneName")





library(pheatmap)
rownames(tmp_file_dat) = tmp_file_dat[,1]
pheatmap(tmp_file_dat[,c(2:7)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         cluster_cols = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

pheatmap(tmp_file_dat[,c(8:13)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

pheatmap(tmp_file_dat[,c(2:13)],
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         cluster_cols = F,
         scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)

write.csv(file = "../Figures/UnionDLGenesName.csv",x = tmp_file_dat[,1])

write.csv(file = "../Figures/IntersectDLGenesName.csv",x = tmp_file_dat[,1])
##########################################Non-samoke Female data###############################
#Based on the sequencing of A549 cell lines, PDIA3P1 is the target of research
###########################################TCGA set#############################################


#####################################################################
#TCGA dataset of LUAD to generate the KEGG and GSEA pthway (geneset)#
#####################################################################

BiocManager::install("TCGAbiolinks")
install.packages("RSQLite")
install.packages("mime")
install.packages("curl")
install.packages("TCGAbiolinks")
install.packages("cmprsk")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
install.packages("devtools")
library("devtools")
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
library("TCGAbiolinks")

query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)

rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
LUADRnaseqSEvalues <- assay(LUADRnaseqSE)
ncol(LUADRnaseqSEvalues)
head(exp.hg38.values)

exp.hg38.values = read_csv("LUADMatrix.csv")

exp.hg38.values_1 = exp.hg38.values


query.exp.hg38 <- GDCquery(project = "TCGA-THCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
query.exp.hg38 <- GDCprepare(query.exp.hg38)
THCARnaseqSE = assay(query.exp.hg38)
rownames(THCARnaseqSE) <- values(query.exp.hg38 )$external_gene_name
THCARnaseqSEValue <- assay(THCARnaseqSE)
head(THCARnaseqSE)
ncol(THCARnaseqSE)
write.csv(x = THCARnaseqSE,file = "THCARnaseqV")


query.exp.hg19 <- GDCquery(project = "TCGA-THCA", 
                           data.category = "Gene expression", 
                           data.type = "Gene expression quantification",
                           legacy = TRUE)
GDCdownload(query.exp.hg38)
exp.hg19 <- GDCprepare(query = query.exp.hg19)
exp.values <- assay(query.exp.hg19)

exp.hg38.values$X1

exp.hg38.values_targeted_gene = exp.hg38.values[exp.hg38.values$X1%in%toupper(gene_name$gene_name),]

gene_name_exp_carcer = exp.hg38.values_targeted_gene[,grepl(pattern = "01A",x = colnames(exp.hg38.values_targeted_gene))]

rownames(gene_name_exp_carcer) =  exp.hg38.values_targeted_gene$X1

gene_name_exp_carcer = as.data.frame(gene_name_exp_carcer)



#####################################################################
#Clinical information                                               #
#####################################################################
library(data.table)
library(dplyr)
install.packages("regexPipes")
library(regexPipes)

clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>% 
  regexPipes::grep("TCGA",value=T) %>% 
  sort %>% 
  plyr::alply(1,GDCquery_clinic, .progress = "text") %>% 
  rbindlist
?alply

clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

clinical_THCA <- GDCquery_clinic(project = "TCGA-THCA", type = "clinical")
head(clinical_THCA)
write.csv(clinical_THCA,file = "clinical_THCA.csv")

clinical_LUAD = clinical_LUAD[!is.na(clinical_LUAD[,2]),]
head(clinical_LUAD)


clinical_md= clinical_LUAD [,c("submitter_id","primary_diagnosis","tumor_stage","age_at_diagnosis",
                  "days_to_last_follow_up","days_to_death","gender")]

#TCGA-05-4244 

clinical_md$age_at_diagnosis=clinical_md$age_at_diagnosis/365
clinical_md$days_to_last_follow_up=clinical_md$days_to_last_follow_up/30
clinical_md$days_to_death = clinical_md$days_to_death/30
clinical_md$outcome = as.factor(if_else(is.na(clinical_md$days_to_last_follow_up),1,0))


tab.noby <- tableby(~primary_diagnosis+tumor_stage+age_at_diagnosis+
                      outcome+days_to_last_follow_up+days_to_death+gender,data=clinical_md)
summary(tab.noby)

DL_state = apply(tmp_file[,c(5:7,11:13)],1,mean)
normal_state = apply(tmp_file[,c(2:4,8:10)],1,mean)
diff_DL_nor =as.data.frame(DL_state - normal_state)
rownames(diff_DL_nor)=tmp_file$geneName

#1. DL increase; 0. DL decease

DL_sign = c()
diff_DL_nor = ifelse(diff_DL_nor<0,0,1)
diff_DL_nor= diff_DL_nor[order(rownames(diff_DL_nor)),]
DL_sign = diff_DL_nor

gene_name_exp_carcer_sort = 
  gene_name_exp_carcer[order(rownames(gene_name_exp_carcer)),]

gene_name_exp_carcer_sign=gene_name_exp_carcer_sort


for (i in 1:length(DL_sign_sort)){
  
  if (DL_sign[i]==1)# DL increase
  {
    gene_name_exp_carcer_sign[i,] = 
      gene_name_exp_carcer_sort[i,] - 
      as.numeric(fivenum(gene_name_exp_carcer_sort[i,])[3]) 
    gene_name_exp_carcer_sign[i,] = 
      ifelse(gene_name_exp_carcer_sign[i,]<0,0,1)
  }else# DL decrease
  {
    gene_name_exp_carcer_sign[i,] = 
      gene_name_exp_carcer_sort[i,] - 
      as.numeric(fivenum(gene_name_exp_carcer_sort[i,])[3]) 
    gene_name_exp_carcer_sign[i,] = ifelse(gene_name_exp_carcer_sign[i,]<0,1,0)
  }
  #  gene_name_exp_carcer_sign[,i] 
}
gene_name_exp_carcer_sign_sum = apply( gene_name_exp_carcer_sign,2,sum)

names = names(gene_name_exp_carcer_sign_sum)
name = c()
for (i in 1:length(gene_name_exp_carcer_sign_sum)){
  tmp = unlist(strsplit(names[i],split = "-"))
  name[i]=paste(tmp[1],tmp[2],tmp[3],sep = "-")
}
DL_statu = rbind(gene_name_exp_carcer_sign_sum,gene_name_exp_carcer_sign)

t(DL_statu)

DL_statu = cbind(name,t(DL_statu)) 



head(clinical_md)


clin_DL = merge(DL_statu,clinical_md,by.x = "name",by.y="submitter_id")

for (i in 1:nrow(clin_DL))
clin_DL$surMonth[i] = ifelse(clin_DL$outcome[i]==1,clin_DL$days_to_death[i],clin_DL$days_to_last_follow_up[i])

clin_DL$surMonth=as.numeric(clin_DL$surMonth)
clin_DL$DL = as.numeric(clin_DL$DL)

clin_DL$outcome = as.numeric(clin_DL$outcome)



head(clin_DL)

library(survival)
library(ggplot2)
require("survival")
BiocManager::install("survminer")
BiocManager::install("ggpubr")
install.packages("survminer")
library(survminer)
BiocManager::install("ggsurvplot")
#library(ggsurvplot)
library(ggpubr)
?survdiff
??survival
head(clin_DL)

stage = as.character(clin_DL$tumor_stage)
stage_simple = c()

for(i in 1:length(stage)){
  if(base::grepl('iii|iv', stage[i]))
  {stage_simple[i]=2}
  else
  {stage_simple[i]=1}
}



clin_DL[,2] = as.numeric(clin_DL[,2])

clin_DL_stage_1 = clin_DL[stage_simple==1,]
clin_DL_stage_2 = clin_DL[stage_simple==2,]

tmp=c()
tmpP=c()
j=1
i=3
for(i in 3:25){
  fit <- coxph(Surv(surMonth,outcome) ~  clin_DL_stage_1[,i] ,data=clin_DL[stage_simple==1,]) 
  a = summary(fit)
  tmp[j] = a$coefficients[2] 
  tmpP[j]=a$coefficients[5]
  j=j+1
}


tmp

tmpP

tmp[tmpP<0.05]

tmp[tmpP<0.05]

clin_DL_stage_1_marker = clin_DL_stage_1[,3:25]


as.numeric(as.matrix(clin_DL_stage_1_marker[i,tmpP<0.05]))*tmp[tmpP<0.05]



as.numeric(as.character(clin_DL_stage_1_marker[i,tmpP<0.05]))

i=1
for (i in 1:ncol(clin_DL_stage_1)){
  clin_DL_stage_1[i,2] =  sum(as.numeric(as.matrix(clin_DL_stage_1_marker[i,tmpP<0.05]))*tmp[tmpP<0.05])
}


for (i in 1:ncol(clin_DL_stage_1)){
  clin_DL_stage_1[i,2] =  sum(as.numeric(clin_DL_stage_1_marker[i,tmpP<0.05]))
}


fit <- coxph(Surv(surMonth,outcome) ~  clin_DL_stage_1[,2] ,data=clin_DL[stage_simple==1,]) 
summary(fit)

survival::survdiff(Surv(clin_DL$surMonth,clin_DL$outcome) ~ DL  ,data=clin_DL)



fit <- coxph(Surv(survial_day, survial_state)~
               as.numeric(gene_name_exp_carcer_sign_sum),data=clin_DL) 
summary(fit)

nrow(clin_DL)

fit<- survfit(Surv(surMonth,outcome) ~DL, data=clin_DL)

ggsurvplot(fit, data = clin_DL_60)
summary(fit)
??ggsurvplot


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







#########################################Demography
install.packages("tidyverse","gapminder")
BiocManager::install("gapminder")
library(tidyverse)
library(gapminder)
data(gapminder)

install.packages("qwraps2")
library(qwraps2)
options(qwraps2_markup = "markdown")
gapminder <- as.data.frame(gapminder)
summary_statistics <-
  list(
    "Life Expectancy" =
      list(
        "mean (sd)" = ~qwraps2::mean_sd(lifeExp, na_rm = TRUE),
        "median (Q1, Q3)" = ~qwraps2::median_iqr(lifeExp, na_rm = TRUE),
        "min" = ~min(lifeExp, na.rm = TRUE),
        "max" = ~max(lifeExp, na.rm = TRUE),
        "Missing" = ~sum(is.na(lifeExp))
      ),
    "Population" =
      list(
        "mean (sd)" = ~qwraps2::mean_sd(pop, na_rm = TRUE),
        "median (Q1, Q3)" = ~qwraps2::median_iqr(pop, na_rm = TRUE),
        "min" = ~min(pop, na.rm = TRUE),
        "max" = ~max(pop, na.rm = TRUE),
        "Missing" = ~sum(is.na(pop))
      ),
    "GDP per Capita" =
      list(
        "High GDP per Capita" = ~qwraps2::n_perc(na.omit(gdpPercap) %in% "high"),
        "Low GDP per Capita" = ~qwraps2::n_perc(na.omit(gdpPercap) %in% "low"),
        "Missing" = ~sum(is.na(gdpPercap))
      )
  )

s = summary_table(gapminder, summary_statistics)

our_summaries <-
  list("Miles Per Gallon" = 
         list("min"  = ~ min(.data$mpg),
              "mean" = ~ mean(.data$mpg),
              "mean &plusmn; sd" = ~ qwraps2::mean_sd(.data$mpg),
              "max"  = ~ max(.data$mpg)),
       "Weight" = 
         list("median" = ~ median(.data$wt)),
       "Cylinders" = 
         list("4 cyl: n (%)" = ~ qwraps2::n_perc0(.data$cyl == 4),
              "6 cyl: n (%)" = ~ qwraps2::n_perc0(.data$cyl == 6),
              "8 cyl: n (%)" = ~ qwraps2::n_perc0(.data$cyl == 8)))

whole_table <- summary_table(mtcars, our_summaries)
whole_table

tab.noby <- tableby(~ bmi + sex + age, data=mockstudy)
summary(tab.noby)


pings.

table_four <- tableby(~year + continent + lifeExp + gdpPercap + pop, data = gapminder)
summary(table_four)

table_four <- tableby(~year + continent + lifeExp + gdpPercap + pop, data = gapminder)
summary(table_four)

BiocManager::install("arsenal")
library(arsenal)
require(knitr)
require(survival)
data(mockstudy) ##load data
dim(mockstudy)  ##look at how many subjects and variables are in the dataset 
summary(tab1, text=TRUE)
tab1 <- tableby(arm ~ sex + age, data=mockstudy)
summary(tab1)
as.data.frame(tab1)

tab.noby <- tableby(~ bmi + sex + age, data=mockstudy)
summary(tab.noby)

summary(tab.test, pfootnote=TRUE)

install.packages("rmarkdown")
library(rmarkdown)
render("1-example.Rmd")

#pdia3p1 Sign#########################################################################

targetGene="PDIA3P1"

##selecting the cancer sample

exp.hg38.values_caner = exp.hg38.values[,grep("01A",colnames(exp.hg38.values))]

?ifelse
# high:1 low:0
sign = ifelse(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]
              >median(exp.hg38.values_caner[rownames(exp.hg38.values_caner) == targetGene,]),1,0)

length(sign)#519

#explanation: the first sign is the ref, namely the denominator####################################
t1 = edgeR::DGEList(exp.hg38.values_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
write.csv(tableDEA,file = "PDIA3P1TCGA.csv")
getwd()

####################################CHECK
?boxplot
dataCheck = exp.hg38.values_caner[rownames(exp.hg38.values_caner) == "NPY",]
boxplot(log(dataCheck) ~ sign,col = "lightgray")


###########################################GSEA results
BiocManager::install("fgsea")
library(fgsea)
#load("../human_H_v5p2.rdata")
load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = row.names(tableDEA)
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
TCGAGsea = fgseaRes[order(padj, -abs(NES)), ] ################## result

topUp <- TCGAGsea %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- TCGAGsea %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              TCGAGsea, 
              gseaParam = 0.5)
################## result 
library(clusterProfiler)


head(fgseaRes[order(padj, -abs(NES)), ], n=10)
plotEnrichment(pathwaysH[["BENPORATH_ES_WITH_H3K27ME3"]], ranks)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

TCGApathway  =  kk@result[kk@result$pvalue<0.1,2]
#####################################################################
#non-female dataset generate the KEGG and GSEA pthway (geneset)
#####################################################################

getwd()
setwd("../non_female_smoker_pm2.5/")
dir()

file_list = dir(pattern = "^[1-9]")
i=1
tmp = read.table(file_list[i],header = T,sep = "\t")
dat = tmp[,c(1,10)]
head(dat)
for (i in 2:length(file_list)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
  dat =cbind(dat,tmp[,10])
}
colnames(dat) = c("Genes",gsub("_genes.fpkm_tracking", "", file_list))


targetGene="PDIA3P1"

##selecting the cancer sample

dat_caner = dat[,grep("A1",colnames(dat))]

dat_caner = as.matrix(dat_caner)

head(dat_caner)

sign = ifelse(dat_caner[dat[,1] == targetGene,]
              >median(dat_caner[dat[,1] == targetGene,]),1,0)

#explanation: the "second" sign is the ref, namely the denominator
t1 = edgeR::DGEList(dat_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)

?boxplot
dataCheck = dat_caner[28868,]
boxplot(log(dataCheck) ~ sign,col = "lightgray")

tableDEA$geneName = dat[rownames(tableDEA),1]

write.csv(tableDEA,file = "PDIA3P1non_female_smoker_pm2.5.csv")




load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = dat[row.names(tableDEA),1]
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
NonFemaleNonSmakerGsea = fgseaRes[order(padj, -abs(NES)), ] ################## result
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

NonFemaleNonSmakerpathway  =  kk@result[kk@result$pvalue<0.05,2]
################## result 

#####################################################

library(GEOquery)

gset <- getGEO("GSE86958", GSEMatrix =TRUE, AnnotGPL=TRUE )

getwd()
setwd("../GSE86958//")
dir()




file_list = dir(pattern = "_Tumor_fpkm.txt")
i=1
tmp = read.table(file_list[i],header = T,sep = "\t")
dat = tmp[,c(2,6)]
head(dat)
for (i in 1:length(file_list)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
#  dat =cbind(dat,tmp[,6])
  print(nrow(tmp))  
}

for (i in c(2:4)){
  tmp = read.table(file_list[i],header = T,sep = "\t")
  dat =cbind(dat,tmp[,6])
}

colnames(dat) = c("Genes",gsub("_fpkm.txt", "", file_list[c(1:4)]))



tmp = read.table(file_list[5],header = T,sep = "\t")
dat_1 = tmp[,c(2,6)]
tmp = read.table(file_list[6],header = T,sep = "\t")
dat_1 =cbind(dat_1,tmp[,6])
colnames(dat_1) = c("Genes",gsub("_fpkm.txt", "", file_list[c(5,6)]))

targetGene="PDIA3P"

dat_caner = as.matrix(dat[,c(2:5)])

sign = ifelse(dat_caner[dat[,1] == targetGene,]
              >median(dat_caner[dat[,1] == targetGene,]),1,0)

#explanation: the first sign is the ref, namely the denominator
t1 = edgeR::DGEList(dat_caner,group = as.factor(sign))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$PValue <= 0.01, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)


?boxplot
dataCheck = dat_caner[12005,]
boxplot(log(dataCheck) ~ sign,col = "lightgray")

tableDEA$geneName = dat[rownames(tableDEA),1]

write.csv(tableDEA,file = "PDIA3P1GSE86958.csv")


load("../human_c2_v5p2.rdata") # relative conprehensive

#input logFC vaule and entrezgene
rank = tableDEA$logFC
row.names(tableDEA)

mygenes = dat[row.names(tableDEA),1]
mapping <- getBM(
  attributes = c('entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

names(rank) = mapping[,1]

ranks = rank[na.omit(names(rank))]

barplot(sort(ranks, decreasing = T))
pathwaysH <- Hs.c2
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)

################## result 
GSE86958Gsea = fgseaRes[order(padj, -abs(NES)), ] ################## result
topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              ranks, 
              fgseaRes, 
              gseaParam = 0.5)

kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')

GSE86958pathway  =  kk@result[kk@result$pvalue<0.05,2]


################## result 


intersect(intersect(NonFemaleNonSmakerpathway,GSE86958pathway),TCGApathway)

intersect(GSE86958pathway,TCGApathway)

intersect(intersect(NonFemaleNonSmakerGsea[,1],TCGAGsea[,1]),GSE86958Gsea[,1])

pathwaysH 

############################################################KEGG


BiocManager::install("clusterProfiler")
library(clusterProfiler)
search_kegg_organism('hsa', by='kegg_code')
?enrichKEGG


head(kk, n=10)
data(geneList, package='DOSE')
browseKEGG(kk, 'hsa05034')
BiocManager::install("pathview")

library(pathview)

ranks


pathview(gene.data = ranks, 
         pathway.id = "hsa05034", 
         species = "hsa", 
         limit = list(gene=5, cpd=1))

data(examplePathways)
data(exampleRanks)
#######################################
library(biomaRt)

hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

hsmart

# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- c("TNF", "IL6", "IL1B", "IL10", "CRP", "TGFB1", "CXCL8")
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values = mygenes,
  mart = hsmart
)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mygenes <- c("TNF", "IL6", "IL1B", "IL10", "CRP", "TGFB1", "CXCL8")
getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene"), 
      filters = "hgnc_symbol",
      values = mygenes,
      mart = ensembl) # 前面选好的的数据库

affyids=c("202763_at","209310_s_at","207500_at")
getBM(attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol','chromosome_name','start_position','end_position', 'band'),
      filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)




library(GEOquery)

GSM118720 <- getGEO('GSM118720')
GSM118721 <- getGEO('GSM118721')

GSM118720 <- getGEO(filename=system.file("extdata/GSM118720.soft",package="GeneExpressionSignature"))
#control gene-expression profiles
GSM118721 <- getGEO(filename=system.file("extdata/GSM118721.soft",package="GeneExpressionSignature"))

control <- as.matrix(as.numeric(Table(GSM118721)[,2]))
head(control)
treatment <- as.matrix(as.numeric(Table(GSM118720)[,2]))
ranked_list <-getRLs(control,treatment)
data(exampleSet)
show(exampleSet)






library(dplyr)


topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(10, wt=-padj)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(10, wt=-padj)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-ES)
plotGseaTable(pathwaysH[topPathways$pathway], 
              rank, 
              fgseaRes, 
              gseaParam = 0.5)

BiocManager::install("goseq")
library(goseq)
supportedOrganisms() %>% filter(str_detect(Genome, "mm"))

genes =as.integer(tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR))
tableDEA[tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR),]
names(genes)=rownames(tableDEA[tableDEA$FDR < 0.01 & !is.na(tableDEA$FDR),])


mapping <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
  filters = 'hgnc_symbol',
  values =names(genes),
  mart = hsmart
)
names(genes) = mapping[,1]

genes = genes[na.omit(names(genes))]



pwf <- nullp(genes, "hg19", "ensGene")
?nullp

getlength(names(genes),'hg19','ensGene')

BiocManager::install("org.Hs.eg.db")

goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:CC"))
library(ggplot2)
goResults %>% 
  top_n(10, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, 
             y=term, 
             colour=over_represented_pvalue, 
             size=numDEInCat)) +
  geom_point() +
  expand_limits(x=0) +
  labs(x="Hits (%)", y="GO term", colour="p value", size="Count")


