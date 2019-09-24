#Transcript of RNA Signature

main = function(){
  ######################Setting the work direction
  getwd()
  setwd("..//DL/")
  ######################Reading different exp gene data 
  genes_1975_dl_con = read.table("1975_gene_exp.diff",header = T,sep = "\t")
  genes_A549_dl_con = read.table("A549_gene_exp.diff",header = T,sep = "\t")
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
  
  # gene name order
  gene_name = as.data.frame(sort((diff_gene)))
  colnames(diff_gene_filer_1 )=c("gene_name","gene_id","H1975_foldChange","p-value","A549_foldChange","p-value")
  diff_gene_filer_1_intersection = diff_gene_filer_1
  gene_name = as.data.frame(sort((diff_gene_filer_1[,1])))
  colnames(gene_name)="gene_name"
  #Outputing file
  write.csv(x = diff_gene_filer_1_intersection ,file = "diff_gene_filer_1_intersection_DL.csv")
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
  head(gene_name_exp_carcer) #Generating the dataset of Patients with DL-affacted genes
  
  
  ##################Calculating the DL profile (The gene name and value)
  getwd()
  
  file_list = dir(pattern = "*.fpkm*")
  
  tmp_file = data.frame(1:601)
  for (i in 1:length(file_list)){
    dat_tmp = read.table(file_list[i],header = T,sep = "\t")
    merge_tmp = merge(dat_tmp,gene_name,by.x = "gene_id",by.y  = "gene_name")
    tmp_file =cbind(tmp_file,as.data.frame(merge_tmp$FPKM))
  }

  colnames(tmp_file) = c("No",gsub("_genes.fpkm_tracking", "", file_list))
  tmp_file$gene_id=merge_tmp$gene_id
  head(tmp_file)
  
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
  
  
}

