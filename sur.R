getwd()
setwd("../DL/")
# read different exp gene data 
genes_1975_dl_con = read.table("1975_con_dl/gene_exp.diff",header = T,sep = "\t")
genes_A549_dl_con = read.table("A549_con_dl/gene_exp.diff",header = T,sep = "\t")
diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.05,]$gene,
                  genes_A549_dl_con[genes_A549_dl_con$q_value<0.05,]$gene)

#union
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


###########################


query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
GDCdownload(query.exp.hg38)
GDCdownload(query.exp.hg38,files.per.chunk = 1)
exp.hg38 <- GDCprepare(query =query.exp.hg38 )

exp.hg38.values <- assay(exp.hg38)
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name



colnames(exp.hg38.values)
colnames(exp.hg38)
rownames(exp.hg38.values)
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)

clinical_LUAD_m1 = data.frame(clinical_LUAD$submitter_id,
                              clinical_LUAD$tumor_stage,
                              clinical_LUAD$days_to_last_follow_up,
                              clinical_LUAD$age_at_diagnosis,
                              clinical_LUAD$days_to_death)
head(clinical_LUAD_m1)
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



clinical_LUAD_m1$clinical_LUAD.submitter_id

t_exp.hg38.values=t(exp.hg38.values)


t_exp.hg38.values_ca =as.data.frame(
  t_exp.hg38.values[grepl(rownames(t_exp.hg38.values),pattern = "01A"),])
library(limma)
tmp = strsplit2(rownames(t_exp.hg38.values_ca),split = "-")
t_exp.hg38.values_ca$pName = paste(tmp[,1],tmp[,2],tmp[,3],sep = "-")


clin_gene = merge(clinical_LUAD_m1,t_exp.hg38.values_ca,by.x = "clinical_LUAD.submitter_id","pName")


stage = as.character(clin_gene$clinical_LUAD.tumor_stage)
stage_simple = c()

for(i in 1:length(stage)){
  if(grepl('iii|iv', stage[i]))
  {stage_simple[i]=2}
  else
  {stage_simple[i]=1}
}
clin_gene$stage = stage_simple

colnames(clin_gene)[1:10]
??survfit


clin_gene$clinical_LUAD.age_at_diagnosis=as.numeric(clin_gene$clinical_LUAD.age_at_diagnosis)/365

library(survival)

clin_gene[1,56719+4]

HR_result = data.frame(1,2)
geneNameHR = NULL

for(i in 1:(ncol(t_exp.hg38.values_ca))){
  
  sfit <- coxph(Surv(survial_day, survial_state)
                ~as.numeric(stage)+as.numeric(clinical_LUAD.age_at_diagnosis)+clin_gene[,7+i], data=clin_gene)
  res=summary(sfit)
  geneNameHR[i] = clin_gene[,7+i]
  HR_result[i,2] = res$coefficients[3,5]
  HR_result[i,1] = res$coefficients[3,2]
}
dat = cbind(clin_gene[,7+1:i],HR_result)
colnames(dat)=c("geneName","HR","Pvalue")
write.csv(dat,file = "LUAD_sur_hg19_HR.csv")

gei = ncol(t_exp.hg38.values_ca)




HR_result[HR_result[,2]<0.01,]

HR_result[HR_result[,3]%in%gene_name$gene_name,]
HR_result[,4] = apply(t_exp.hg38.values_ca, 2, max)
write.csv(HR_result[HR_result[,3]%in%gene_name$gene_name|HR_result[HR_result[,2]<0.01,],],file = "sur_dl.csv")
