#Transcript of RNA Signature


main = function(){
  ######################Setting the work direction
  getwd()
  setwd("~/DL/1975_con_dl")
  ######################Reading different exp gene data 
  genes_1975_dl_con = read.table("gene_exp.diff",header = T,sep = "\t")
  genes_A549_dl_con = read.table("gene_exp.diff",header = T,sep = "\t")
  diff_gene = union(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                        genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)#
  diff_gene = intersect(genes_1975_dl_con[genes_1975_dl_con$q_value<0.1,]$gene,
                        genes_A549_dl_con[genes_A549_dl_con$q_value<0.1,]$gene)# extract the foldchange of diff exp gene
  diff_gene_fold = cbind(as.character(genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$gene),
                         as.character(genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$gene),
                         genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                         genes_1975_dl_con[(genes_1975_dl_con$gene)%in%diff_gene,]$q_value,
                         genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$log2.fold_change.,
                         genes_A549_dl_con[(genes_A549_dl_con$gene)%in%diff_gene,]$q_value
  )
  
  #Extracting the data with the same direction 
  diff_gene_filer_1 = diff_gene_fold[as.numeric(diff_gene_fold[,3])*as.numeric(diff_gene_fold[,5])>0,]
  colnames(diff_gene_filer_1 )=c("gene_name","gene_id","H1975_foldChange","p-value","A549_foldChange","p-value")
  gene_name=diff_gene_filer_1[,1]
  gene_name = data.frame(gene_name)
  #Outputing file
  write.csv(x = diff_gene_filer_1 ,file = "diff_gene_filer_1_intersection_DL.csv")
  write.csv(x = diff_gene_filer_1 ,file = "diff_gene_filer_1_Union_DL.csv")
  
  
  #####################Using the TCGA dataset
  
  
  library(SummarizedExperiment)
  library(TCGAbiolinks)
  library(limma)
  query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - FPKM")
  LUADRnaseqSE <- GDCdownload(query.exp.hg38)
  LUADRnaseqSE <- GDCprepare(query.exp.hg38)
  #GDCdownload(query.exp.hg38,files.per.chunk = 1)
  
  rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
  exp.hg38.values <- assay(LUADRnaseqSE)
  
  head(exp.hg38.values)
  #write.csv(exp.hg38.values,file = "stad_exp_hg38_FPKM.csv")
  # extract the targeted gene
  #rownames(exp.hg38.values) = tolower(rownames(exp.hg38.values))
  # gene collection
  exp.hg38.values_targeted_gene = exp.hg38.values[rownames(exp.hg38.values)%in%diff_gene_filer_1[,1],]
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
  head(gene_name_exp_carcer) #Generating the dataset of Patients with DL-affacted genes
  
  
  ##################Calculating the DL profile (The gene name and value)
  getwd()
  setwd("~/DL")
  file_list = dir(pattern = "*.fpkm*")
  
  tmp_file = data.frame(1:601)
  for (i in 1:length(file_list)){
    dat_tmp = read.table(file_list[i],header = T,sep = "\t")
    merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
    tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
  }

  
  colnames(tmp_file) = c("geneId",gsub("_genes.fpkm_tracking", "", file_list))
  tmp_file$gene_id = merge_tmp$gene_id

  
  ##################DL condition 
  source("DLCalculation.R")
  gene_name_exp_carcer_sign_sum = DLCalculation(gene_name_exp_carcer,tmp_file)
  names = colnames(gene_name_exp_carcer_sign_sum)
  DL_statu = as.data.frame(cbind(names,t(gene_name_exp_carcer_sign_sum)))
  
  ##################Clinical Data
  
  clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
  colnames(clinical_LUAD)
  clinical_LUAD$age_at_diagnosis/365
  clinical_LUAD$cigarettes_per_day
  clinical_LUAD$alcohol_intensity#no data
  clinical_LUAD$bmi#no data
  
  
  
  
  
  clinical_LUAD_m1 = data.frame(clinical_LUAD$submitter_id,
                                clinical_LUAD$tumor_stage,
                                clinical_LUAD$days_to_last_follow_up,
                                clinical_LUAD$age_at_diagnosis,
                                clinical_LUAD$days_to_death,
                                clinical_LUAD$age_at_diagnosis/365,
                                clinical_LUAD$cigarettes_per_day
                                )
  
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
  
  clin_DL = merge(DL_statu,clinical_LUAD_m1,by.x = "names",by.y="clinical_LUAD.submitter_id")
  
  
  stage = as.character(clin_DL$clinical_LUAD.tumor_stage)
  stage_simple = c()
  for(i in 1:length(stage)){
    if(base::grepl('iii|iv', stage[i]))
    {stage_simple[i]=2}
    else
    {stage_simple[i]=1}
  }
  
  
  
  library(survival)
  library(ggplot2)
  require("survival")
  library(survminer)
  
  colnames(clin_DL)
  
  
  table_1 = data.frame()
  
  fit <- coxph(Surv(survial_day, survial_state)~ stage_simple ,data=clin_DL)
  sumfit = summary(fit)
  table_1 = sumfit$coefficients
  fit <- coxph(Surv(survial_day, survial_state)~ clinical_LUAD.age_at_diagnosis ,data=clin_DL)
  sumfit = summary(fit)
  table_1 = rbind( table_1 ,sumfit$coefficients)
  fit <- coxph(Surv(survial_day, survial_state)~clinical_LUAD.cigarettes_per_day ,data=clin_DL)
  sumfit = summary(fit) 
  table_1 = rbind( table_1 ,sumfit$coefficients) 
  fit <- coxph(Surv(survial_day, survial_state)~clinical_LUAD.cigarettes_per_day+clinical_LUAD.age_at_diagnosis+stage_simple,data=clin_DL)
  sumfit = summary(fit) 
  table_0 = sumfit$coefficients
  
  
  table_2=data.frame()
  fit <- coxph(Surv(survial_day, survial_state)~ stage_simple + clin_DL[,2],data=clin_DL)
  sumfit = summary(fit)
  table_2 = sumfit$coefficients[2,]
  
  for (i in 3:24){
    fit <- coxph(Surv(survial_day, survial_state)~ stage_simple + clin_DL[,i],data=clin_DL)
    sumfit = summary(fit)
    table_2 = rbind(table_2,sumfit$coefficients[2,])
  }
  
  table_2=as.data.frame(table_2)
  
  rownames(table_2)=colnames(clin_DL)[2:24]
  
  
  table_3=data.frame()
  fit <- coxph(Surv(survial_day, survial_state)~ clin_DL[,2],data=clin_DL)
  sumfit = summary(fit)
  table_3 = sumfit$coefficients[1,]
  
  for (i in 3:24){
    fit <- coxph(Surv(survial_day, survial_state)~ clin_DL[,i],data=clin_DL)
    sumfit = summary(fit)
    table_3 = rbind(table_3,sumfit$coefficients[1,])
  }
  table_3=as.data.frame(table_3)
  
  rownames(table_3)=colnames(clin_DL)[2:24]
  
  
  
  
  
  
  
  
  plot(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum),clin_DL$survial_day)

  clin_DL$group = ifelse(as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)>12,1,0)
  
  clin_DL$gene_name_exp_carcer_sign_sum=as.numeric(clin_DL$gene_name_exp_carcer_sign_sum)
  
  boxplot(clin_DL$survial_day~clin_DL$group)
  
  t.test(clin_DL$survial_day~clin_DL$group)
  

  clin_DL_stage_1 = clin_DL[stage_simple==1,]
  clin_DL_stage_2 = clin_DL[stage_simple==2,]
  
  fit <- coxph(Surv(survial_day, survial_state)~
                 as.numeric(gene_name_exp_carcer_sign_sum),data=clin_DL ) 
  summary(fit)
  
  fit <- coxph(Surv(survial_day, survial_state)~
                 as.numeric(gene_name_exp_carcer_sign_sum),data=clin_DL_stage_1 ) 
  summary(fit)
  
  fit <- coxph(Surv(survial_day, survial_state)~
                 as.numeric(gene_name_exp_carcer_sign_sum),data= clin_DL_stage_2) 
  summary(fit)
  

  require("survival")  
  data_3_year = clin_DL[clin_DL$survial_day<3*365,]
  data_2_year = clin_DL[clin_DL$survial_day<2*365,]
  
  fit<- survfit(Surv(survial_day, survial_state)~data_3_year[,4], data=data_3_year)
  ggsurvplot(fit, data =   data_3_year )
  
  splots <- list()
  splots[[1]] = ggsurvplot(fit, data = data_3_year,
                           font.main = c(16, "bold", "darkblue"),
                           font.x = c(20, "bold.italic", "darkred"),
                           font.y = c(20, "bold.italic", "darkred"),
                           font.tickslab = c(20, "plain", "Black"))
  splots[[2]] <- ggsurvplot(fit, data = data_3_year, risk.table = TRUE, ggtheme = theme_grey())# Arrange multiple ggsurvplots and print the output
  arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
  
  
 survival::survdiff(Surv(survial_day, survial_state)~data_2_year[,4],data=data_2_year)
  
  fit<- survfit(Surv(survial_day, survial_state)~data_2_year[,4], data=data_2_year)
  ggsurvplot(fit, data =   data_3_year )
  
  splots <- list()
  splots[[1]] = ggsurvplot(fit, data = data_2_year,
                           font.main = c(16, "bold", "darkblue"),
                           font.x = c(20, "bold.italic", "darkred"),
                           font.y = c(20, "bold.italic", "darkred"),
                           font.tickslab = c(20, "plain", "Black"))
  splots[[2]] <- ggsurvplot(fit, data = data_2_year, risk.table = TRUE, ggtheme = theme_grey())# Arrange multiple ggsurvplots and print the output
  arrange_ggsurvplots(splots, print = TRUE,  ncol = 2, nrow = 1, risk.table.height = 0.4)
  
  ##############################GEO dataset
  install.packages("BiocManager")
  BiocManager::install("GEOquery")
  library(GEOquery)
  gset <- getGEO("GSE42127", GSEMatrix =TRUE, AnnotGPL=TRUE )
  
  
  exprSet <- exprs(gset[[1]])
  
  exprSet$id = rownames(exprSet)
  
  pData <- pData(gset[[1]])
  
  fdata<-fData(gset[[1]])
  
  fdata_target = fdata[fdata[,3]%in%diff_gene_filer_1[,1],c(1,3)]
  
  exprSet_target = as.data.frame(exprSet[rownames(exprSet)%in%fdata_target[,1],])
  exprSet_target$ID = rownames(exprSet_target)
  
  targetGene = merge(fdata_target,exprSet_target,by.x="ID",by.y="ID")
  targetGeneM = targetGene[,-1]
 
  targetGeneM = (as.data.frame(t(targetGeneM)))
  colnames(targetGeneM) = targetGene[,2]
  targetGeneM = targetGeneM[-1,]
  targetGeneM$ID = rownames(targetGeneM) 
  colnames(targetGeneM)
  targetGeneM[,34]
  
  source("DLCalculation.R")
  gene_name_exp_carcer_sign_sum = DLCalculation(gene_name_exp_carcer,tmp_file)
  
  
  
  
  colnames(pData)
  head(pData)
 
  
  clinData=data.frame(c(c(1:nrow(pData))))
  clinData$stage = pData$`final.pat.stage:ch1`
  clinData$histology = pData$`histology:ch1`
  clinData$month = as.numeric(pData$`overall survival months:ch1`)
  clinData$status = pData$`survival status:ch1`
  
  clinData$status = ifelse(clinData$status=="A",0,1)
  clinData$ID = rownames(pData)

  stage = as.character(clinData$stage)
  stage_simple = c()
  for(i in 1:length(stage)){
    if(base::grepl('III|IV', stage[i]))
    {stage_simple[i]=2}
    else
    {stage_simple[i]=1}
  }
  
  
  fit <- coxph(Surv(month, status)~ stage_simple ,data=clinData)
  sumfit = summary(fit)
  
  
  clinDataDL = merge(targetGeneM,clinData,by.x="ID",by.y="ID")
  
  
  fit <- coxph(Surv(month, status)~ stage_simple+clinData[,1] ,data=  clinDataDL )
  sumfit = summary(fit)
  
  table_4=data.frame()
  for (i in 2:34){
    fit <- coxph(Surv(month, status) ~ stage_simple +  clinDataDL[,i], data=  clinDataDL )
    sumfit = summary(fit)
    table_4 = rbind(table_4,sumfit$coefficients[1,])
  }
  table_4=as.data.frame(table_4)
  
 table_4$genes=colnames(clinDataDL)[1:34]
  
  
  
  
  sample <- pData$geo_accession
  group = ifelse(grepl(pattern = "non",pData[,1]),1,0)
  group <- rep(c(1,0),times=c(12,11))
  design <- model.matrix(~group)
  rownames(design)=colnames(exprSet)
  
  
}





test = function(){}

DLCalculation = function(data,DLSignature){
  #data
  DL_state = apply(DLSignature[,c(5:7,11:13)],1,mean)
  normal_state = apply(DLSignature[,c(2:4,8:10)],1,mean)
  diff_DL_nor =as.data.frame(DL_state - normal_state)
  rownames(diff_DL_nor)=DLSignature$gene_id
  DL_sign = c()
  DL_sign = ifelse(diff_DL_nor<0,0,1)
  DL_sign_sort = DL_sign[order(rownames(DL_sign))]
  gene_name_exp_carcer_sort = 
    data[order(rownames(gene_name_exp_carcer)),]
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
  return( gene_name_exp_carcer_sign_sum )
}
