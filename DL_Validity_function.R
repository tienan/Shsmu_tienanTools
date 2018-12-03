DL_Vadility <- function(CancerName) {
  print("read DL genes")
###################################
  # intersect
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
  colnames(tmp_file)=c("id",file_list,"gene_id")
#######################################################
  print("read RNA genes")
######################################################
  library("limma")
  ####TCGA RNA data 
  #install.packages("devtools")
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
  print("DL condition")
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
  
  clinical_LUAD <- GDCquery_clinic(project = CancerName, type = "clinical")
  colnames(clinical_LUAD)
  #clinical_LUAD$
  #临床数据整理
  print("clinical_data read")
  clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")
  
  clinical_LUAD_m1 = data.frame(clinical_LUAD$submitter_id,
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
  clinical_LUAD_m1$clinical_LUAD.submitter_id
  clin_DL = merge(DL_statu,clinical_LUAD_m1,by.x = "name",by.y="clinical_LUAD.submitter_id")
  #head(clin_DL)
  #cbind(clin_DL$DL_level,clin_DL$survial_day,clin_DL$clinical_LUAD.tumor_stage)
  
  #plot(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum),clin_DL$survial_day)
  ##########cancer stage modify 
  stage = as.character(clin_DL$clinical_LUAD.tumor_stage)
  stage_simple = c()
  ?grepl
  for(i in 1:length(stage)){
    if(grepl('iii|iv', stage[i]))
    {stage_simple[i]=2}
    else
    {stage_simple[i]=1}
  }
  clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>10,1,0)
  
  #boxplot(clin_DL$survial_day~clin_DL$group)
  
  print(t.test(clin_DL$survial_day~clin_DL$group))
  
  library(survival)
  library(ggplot2)
  require("survival")
  library(survminer)
  surR = survival::survdiff(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL)
#  fit <- coxph(Surv(survial_day, survial_state)~group+stage_simple,data=clin_DL) 
  fit<- survfit(Surv(survial_day, survial_state)~group+stage_simple, data=clin_DL)
  tiff(filename = paste(CancerName,".tif"),
       width = 2480, height = 2480, units = "px", pointsize = 12,
       compression = "lzw", 
       bg = "white", res = 300
       )
  ggsurvplot(fit, data = clin_DL)
  dev.off()
# summary(fit)
  return(1-pchisq(surR$chisq,3))
}
