#############################LUAD DEGs
#hg38
main = function(){
  source("/home/tienan/R/Shsmu_tienanTools/TGGA_analysis.R")
  library(limma)
  names = read.table("/home/tienan/R/Shsmu_tienanTools/TCGA-Cancer_1.txt")
  for(i in 1:length(names))DEG(names[i])
  DEG("LUAD")
  }


DEG <- function(name) {
#name = "LUAD" name  name="THCA"
  setwd("/media/tienan/00006784000048231/R/DEGs/")
  
  library(SummarizedExperiment)
  library(TCGAbiolinks)
  library(limma)
  
  query.exp.hg38 <- GDCquery(project = paste("TCGA",name,sep = "-"), #TCGA-LUAD
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "HTSeq - FPKM")
  GDCdownload(query.exp.hg38,files.per.chunk = 50,method = "client")
  LUADRnaseq <- GDCprepare(query.exp.hg38)
  LUADMatix <- assay(LUADRnaseq)
  rownames(LUADMatix) <- values(LUADRnaseq)$external_gene_name
  write.csv(LUADMatix,file = paste(name,"Matrix.csv"))
  
  
  sign=NULL
  patient_id = colnames(LUADMatix)
  for (i in 1:length(patient_id)){
    tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
    sign[i] = tmp[4]
    tmp = paste(tmp[1],tmp[2],tmp[3],sep = "_")
    patient_id[i] = tmp
  }
  
  
  gene_name_exp_dif = LUADMatix[,sign == "01A"|sign == "11A"]
  rownames(gene_name_exp_dif)
  
  
  t1 = edgeR::DGEList(gene_name_exp_dif,group = as.factor(sign[sign == "01A"|sign == "11A"]))
  t2 = edgeR::estimateCommonDisp(t1)
  t3 = edgeR::exactTest(t2)
  t3$table$logFC=-t3$table$logFC
  dat = cbind(rownames(gene_name_exp_dif),t3$table)
  dat[dat$`rownames(gene_name_exp_dif)`=="PDIA3P1",]
  
  nrow(dat[dat[,4]<0.05,])
  
  nrow(dat[dat[,4]<0.05,])/nrow(dat)
  setwd("/home/tienan/R/DL/")
  write.csv(dat,file = paste(name,"DEG_hg38_edgeR.csv",sep = "_"))
  
  ####################################################################################
  main_1 = function(){
  LUADMatrix = LUADMatix
  
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
                              fdr.cut = 0.1 ,
                              logFC.cut = 0.5,
                              method = "glmLRT")
  
  # DEGs table with expression values in normal and tumor samples
  dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                            dataFilt[,samplesTP],dataFilt[,samplesNT])
  
  
  nrow(dataDEGsFiltLevel)
  nrow(dataDEGsFiltLevel)/nrow(dataNorm)
  
    
  setwd("/home/tienan/R/DL/")  
  write.csv(dataDEGsFiltLevel,file = paste(name,"DEG_hg38_TCGAbiolinks.csv",sep = "_"))

}
#############################################################


########################################A549 H1975

#H1975 con vs dl
setwd("/home/tienan/R/DL/")
dat = read.csv("salmon_exp_1st_seq.csv")
colnames(dat)
gene_H1975_DLvsCon = dat[,c(3:8)]
group=NULL
group[1:3]="Con"
group[4:6]="DL"
t1 = edgeR::DGEList(gene_H1975_DLvsCon,group = group)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$GeneName,t3$table)

write.csv(resDat,file = "H1975ConvsDL.csv")

write.csv(resDat[resDat$PValue<0.05,],file = "H1975ConvsDL_sig.csv")


#A549 con vs dl
setwd("/home/tienan/R/DL/")
dat = read.csv("salmon_exp_1st_seq.csv")
colnames(dat)
gene_A549_DLvsCon = dat[,c(24:29)]
group=NULL
group[1:3]="Con"
group[4:6]="DL"
t1 = edgeR::DGEList(gene_A549_DLvsCon,group = group)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)

resDat = cbind(dat$GeneName,t3$table)
write.csv(resDat,file = "A549ConvsDL.csv")

write.csv(resDat[resDat$PValue<0.05,],file = "A549ConvsDL_sig.csv")



#######################################


########################################BASE2B  CT is the wild-type 
setwd("/home/tienan/R/DL/")
dat = read.csv("/home/tienan/R/DL/base2bexp.csv")
colnames(dat)

# BASE2B CC P10   Con vs DL

gene_Base2b_ConvsDL_p10 = dat[,c(37:44,5:12)] # change 
group=NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(gene_Base2b_ConvsDL_p10,group = group) # change 
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_ConvsDL_p10_CC.csv") #change 


# BASE2B CC P30   Con vs DL

gene_Base2b_ConvsDL_p30 = dat[,c(45:52,13:20)] # change 
group=NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(gene_Base2b_ConvsDL_p30,group = group) # change 
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_ConvsDL_p30_CC .csv") #change 

# BASE2B CT P10   Con vs DL

gene_Base2b_ConvsDL_p30 = dat[,c(117:124,77:84)] # change 
group=NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(gene_Base2b_ConvsDL_p30,group = group) # change 
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_ConvsDL_p10_CT.csv") #change 

write.csv(resDat[resDat$PValue<0.05,],file = "gene_Base2b_ConvsDL_p10_CT_sig.csv")


# BASE2B CT P30   Con vs DL

gene_Base2b_ConvsDL_p30 = dat[,c(125:132,85:92)] # change 
group=NULL
group[1:8]="Con"
group[9:16]="DL"
t1 = edgeR::DGEList(gene_Base2b_ConvsDL_p30,group = group) # change 
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_ConvsDL_p30_CT.csv") #change 

write.csv(resDat[resDat$PValue<0.05,],file = "gene_Base2b_ConvsDL_p30_CT_sig.csv")

####################################CC    
# P10 vs P0 CC 

tmp = dat[,c(69:76,37:44)] # change
group=NULL
group[1:8]="P0"
group[9:16]="P10"
t1 = edgeR::DGEList(tmp,group = group) # c
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P10vspP0_CC.csv") #change 

# P30 vs P0 CC 

tmp = dat[,c(69:76,45:52)] # change
group=NULL
group[1:8]="P0"
group[9:16]="P30"
t1 = edgeR::DGEList(tmp,group = group) # change
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P30vspP0_CC.csv") #change 


# P30 vs P10 CC 

tmp = dat[,c(37:44,45:52)] # change
group=NULL
group[1:8]="P10"
group[9:16]="P30"
t1 = edgeR::DGEList(tmp,group = group) # change
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P10vspP30_CC.csv") #change 

#########################CT

# P10 vs P0 CT 

tmp = dat[,c(109:116,117:124)]
group=NULL
group[1:8]="P0"
group[9:16]="P10"
t1 = edgeR::DGEList(gene_Base2b_ConvsDL_p30,group = group) # change 
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P10vsP0_CT.csv") #change 

# P30 vs P0 CT 

tmp = dat[,c(109:116,125:132)] # change
group=NULL
group[1:8]="Con"
group[9:16]="P30"
t1 = edgeR::DGEList(tmp,group = group) # change
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P30vsP0_CT.csv") #change 


# P30 vs P10 CT 

tmp = dat[,c(117:124,125:132)] # change
group=NULL
group[1:8]="P10"
group[9:16]="P30"
t1 = edgeR::DGEList(tmp,group = group) # change
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
resDat = cbind(dat$Name,t3$table)
write.csv(resDat,file = "gene_Base2b_P10vspP30_CT.csv") #change 

################################################################### analysis


setwd("/media/tienan/00006784000048231/R/DEGs/")

LUAD_DEG_hg38_edgeR = read.csv("LUAD_DEG_hg38_edgeR.csv")

LUAD_TCGAbiolinks = read.csv("LUAD_TCGAbiolinks.csv")

intersect(LUAD_DEG_hg38_edgeR$rownames.gene_name_exp_dif.,LUAD_TCGAbiolinks$mRNA)

setwd("/home/tienan/R/DL/")
list=dir()
list_1 =list[grepl(list,pattern = "gene_Base2b*")]
list_2 =list[grepl(list,pattern = "*.csv")]


H1975_DL_DEGs = read.csv("H1975ConvsDL.csv")
head(H1975_DL_DEGs)
nrow(H1975_DL_DEGs[H1975_DL_DEGs$PValue<0.05,])
nrow(H1975_DL_DEGs[H1975_DL_DEGs$PValue<0.05&H1975_DL_DEGs$logFC>0,])
nrow(H1975_DL_DEGs[H1975_DL_DEGs$PValue<0.05&H1975_DL_DEGs$logFC<0,])

A549ConvsDL = read.csv("A549ConvsDL.csv")
head(A549ConvsDL)
nrow(A549ConvsDL[A549ConvsDL$PValue<0.05,])
nrow(A549ConvsDL[A549ConvsDL$PValue<0.05&A549ConvsDL$logFC>0,])
nrow(A549ConvsDL[A549ConvsDL$PValue<0.05&A549ConvsDL$logFC<0,])


tmp = read.csv("gene_Base2b_ConvsDL_p10_CT.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_ConvsDL_p30_CT.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_ConvsDL_p10_CC.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_ConvsDL_p30_CC .csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_P10vsP0_CT.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_P30vsP0_CT.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])


tmp = read.csv("gene_Base2b_P10vspP30_CT.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_P10vspP0_CC.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])

tmp = read.csv("gene_Base2b_P30vspP0_CC.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])


tmp = read.csv("gene_Base2b_P10vspP30_CC.csv")
head(tmp)
nrow(tmp[tmp$PValue<0.05,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC>0,])
nrow(tmp[tmp$PValue<0.05&tmp$logFC<0,])


####################################################analysis#################################################

LUAD_1 = read.csv("LUAU_DEG_hg38_edgeR.csv")

LUAD_1 = LUAD_1[LUAD_1$PValue<0.05,]

LUAD_2 = read.csv("/media/tienan/0000678400004823//R/DEGs/LUAD_TCGAbiolinks.csv")

H1975_DL_DEGs = read.csv("H1975ConvsDL.csv")

A549_DL_DEGs = read.csv("A549ConvsDL.csv")

write.csv(cbind(H1975_DL_DEGs,A549_DL_DEGs),file = "A549_H1975_DL.csv")

#
gene_DL = intersect(H1975_DL_DEGs[H1975_DL_DEGs$PValue<0.1,]$dat.GeneName,A549_DL_DEGs[A549_DL_DEGs$PValue<0.1,]$dat.GeneName)


dat = read.csv("salmon_exp_1st_seq.csv")

gene_A549_H1975_DLvsCon = dat[,c(2,3:8,24:29)]


gene_A549_H1975_DLvsCon_target = gene_A549_H1975_DLvsCon[gene_A549_H1975_DLvsCon$GeneName%in%gene_DL,][,c(1:13)]


colnames(gene_A549_H1975_DLvsCon_target)=
  c("geneName","H1975_Control_1","H1975_Control_2","H1975_Control_3","H1975_DL_1","H1975_DL_2",
    "X1975_DL_3","A549_Control_1","A549_Control_2","A549_Control_3","A549_DL_1","A549_DL_2",
    "A549_DL_3")

#rownames(gene_A549_H1975_DLvsCon_target) = gene_A549_H1975_DLvsCon_target[,1]

setwd("/home/tienan/R/DL/")

#install.packages("pheatmap")
library(pheatmap)
tiff(filename = "Figure-H1975.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400)

pheatmap(gene_A549_H1975_DLvsCon_target[,c(2:7)],
         clustering_distance_cols  = "euclidean",
         show_colnames =   T,
         show_rownames = F,
         scale = "row",
         cluster_cols = T,
         cluster_rows = F,
         fontsize_row = 3, 
         fontsize_col = 12
)

dev.off()

tiff(filename = "Figure-A549.tif",
     width = 1800, height = 3000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400)

pheatmap(gene_A549_H1975_DLvsCon_target[,c(8:13)],
         clustering_distance_cols  = "euclidean",
         show_colnames =   T,
         show_rownames = F,
         scale = "row",
         cluster_cols = T,
         cluster_rows = F,
         fontsize_row = 3, 
         fontsize_col = 12
)

dev.off()


###################################
gene_DL = intersect(H1975_DL_DEGs[H1975_DL_DEGs$PValue<0.1,]$dat.GeneName,A549_DL_DEGs[A549_DL_DEGs$PValue<0.1,]$dat.GeneName)

intersect(diff_gene_filer_1[,1],gene_DL)

cbind(H1975_DL_DEGs,A549_DL_DEGs)[cbind(H1975_DL_DEGs,A549_DL_DEGs)[,2]%in%diff_gene_filer_1[,1],]


mean(cbind(H1975_DL_DEGs,A549_DL_DEGs)[,4])
mean(cbind(H1975_DL_DEGs,A549_DL_DEGs)[,9])
########################################

H1975_subset = H1975_DL_DEGs[H1975_DL_DEGs$dat.GeneName%in%gene_DL,]
nrow(H1975_subset)
H1975_subset = H1975_subset[order(H1975_subset$dat.GeneName),]

A549_subset = A549_DL_DEGs[A549_DL_DEGs$dat.GeneName%in%gene_DL,]
nrow(A549_subset)
A549_subset = A549_subset[order(A549_subset$dat.GeneName),]

DL_genes = H1975_subset[H1975_subset$logFC*A549_subset$logFC>0,]
nrow(DL_genes)
setwd("/home/tienan/R/DL/")
#output data##############################
write.csv(DL_genes,file = "H1975A549_DL_diff_0.1.csv")



###############################LUAD_1
LUAD_1_DL_gene = intersect(LUAD_1$rownames.gene_name_exp_dif.,DL_genes$dat.GeneName)

LUAD_1_subset = LUAD_1[LUAD_1$rownames.gene_name_exp_dif.%in%LUAD_1_DL_gene,]
nrow(LUAD_1_subset)
LUAD_1_subset_1= LUAD_1_subset[LUAD_1_subset$rownames.gene_name_exp_dif.!="CLIC1",]
nrow(LUAD_1_subset_1)

DL_genes_subset = DL_genes[DL_genes$dat.GeneName%in%LUAD_1_DL_gene,]
nrow(DL_genes_subset)
DL_genes_subset = DL_genes_subset[order(DL_genes_subset$dat.GeneName),]
DL_genes_subset_1 = DL_genes_subset[DL_genes_subset$dat.GeneName!="CLIC1",]
nrow(DL_genes_subset_1)


LUAD_DL_geneset =  DL_genes_subset_1[LUAD_1_subset_1$logFC*DL_genes_subset_1$logFC<0,]
nrow(LUAD_DL_geneset)



write.csv(LUAD_DL_geneset,file = "H1975A549_DL_diff_cancer.csv")



head(LUAD_sur_hg19_HR_DLset_modify_sig)

LUAD_sur_hg19_HR_DLset_modify_sig$logHR = log(LUAD_sur_hg19_HR_DLset_modify_sig$HR)


###############################LUAD_2

LUAD_2_DL_gene = intersect(LUAD_2$mRNA,DL_genes$dat.GeneName)

# LUAD_1_subset = LUAD_1[LUAD_1$rownames.gene_name_exp_dif.%in%LUAD_1_DL_gene,]
# nrow(LUAD_1_subset)
# LUAD_1_subset_1= LUAD_1_subset[LUAD_1_subset$rownames.gene_name_exp_dif.!="CLIC1",]
# nrow(LUAD_1_subset_1)
# 
# DL_genes_subset = DL_genes[DL_genes$dat.GeneName%in%LUAD_1_DL_gene,]
# nrow(DL_genes_subset)
# DL_genes_subset = DL_genes_subset[order(DL_genes_subset$dat.GeneName),]
# DL_genes_subset_1 = DL_genes_subset[DL_genes_subset$dat.GeneName!="CLIC1",]
# nrow(DL_genes_subset_1)
# 
# 
# LUAD_DL_geneset =  LUAD_1_subset[LUAD_1_subset_1$logFC*DL_genes_subset_1$logFC<0,]
# nrow(LUAD_DL_geneset)
# 
# 
# write.csv(LUAD_DL_geneset,file = "LUAD_DL_geneset.csv")


intersect(LUAD_2_DL_gene,LUAD_DL_geneset$rownames.gene_name_exp_dif.)


LUAD_DL_geneset$dat.GeneName
LUAD_DL_geneset$logFC
write.csv(LUAD_DL_geneset,file = "LUAD_DL_geneset.csv")
nrow(LUAD_DL_geneset[LUAD_DL_geneset$logFC>0,])
nrow(LUAD_DL_geneset[LUAD_DL_geneset$logFC<0,])


}


###########################clinical data analysis

clinical =function(){
#  write.csv(LUADMatix,file = paste(name,"Matrix.csv",sep = ""))
  getwd()
  setwd("/media/tienan/00006784000048231/R/DEGs/")
  LUADMatix = read.csv("LUADMatrix.csv")
  head(LUADMatix)
  LUAD_sur_hg19_HR = read.csv("/home/tienan/R/DL/LUAD_sur_hg19_HR.csv")
  head(LUAD_sur_hg19_HR)
  
  LUAD_sur_hg19_HR_DLset =  LUAD_sur_hg19_HR[LUAD_sur_hg19_HR$geneName%in%LUAD_DL_geneset$dat.GeneName,]
  nrow(LUAD_sur_hg19_HR_DLset)
  LUAD_sur_hg19_HR_DLset = LUAD_sur_hg19_HR_DLset[order(LUAD_sur_hg19_HR_DLset$geneName),]
  LUAD_sur_hg19_HR_DLset_modify = LUAD_sur_hg19_HR_DLset[(LUAD_sur_hg19_HR_DLset$HR-1)*LUAD_DL_geneset$logFC<0,]
  LUAD_sur_hg19_HR_DLset_modify_sig = LUAD_sur_hg19_HR_DLset_modify[LUAD_sur_hg19_HR_DLset_modify$Pvalue<0.05,]
  nrow(LUAD_sur_hg19_HR_DLset_modify_sig)
  write.csv(LUAD_sur_hg19_HR_DLset_modify_sig,file = "LUAD_sur_hg19_HR_DLset_modify_sig.csv")
  
  clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
  colnames(clinical_LUAD)
 
  
  
  
}

##########################################################prevention analysis


#####CAR

#LUAD_1 = read.csv("LUAU_DEG_hg38_edgeR.csv")

#LUAD_2 = read.csv("/media/tienan/0000678400004823//R/DEGs/LUAD_TCGAbiolinks.csv")

LUAD_sur_hg19_HR_DLset_modify_sig = 
  read.csv("/media/tienan/0000678400004823/R/DEGs/LUAD_sur_hg19_HR_DLset_modify_sig.csv")

########## DL_R

gene_Base2b_ConvsDL_p30_CT = read.csv("gene_Base2b_ConvsDL_p30_CT.csv")

gene_Base2b_ConvsDL_p30_CT

nrow(gene_Base2b_ConvsDL_p30_CT)

gene_Base2b_ConvsDL_p10_CT = read.csv("gene_Base2b_ConvsDL_p10_CT.csv")


gene_Base2b_ConvsDL = 
  intersect(gene_Base2b_ConvsDL_p30_CT[gene_Base2b_ConvsDL_p30_CT$PValue<0.05,]$dat.Name,
                    gene_Base2b_ConvsDL_p30_CT[gene_Base2b_ConvsDL_p10_CT$PValue<0.05,]$dat.Name)

gene_Base2b_ConvsDL_same_direction = 
  gene_Base2b_ConvsDL_p30_CT[gene_Base2b_ConvsDL_p30_CT$logFC*gene_Base2b_ConvsDL_p10_CT$logFC>0,]


c("ACAT2","DHCR24","EBP",  "SCD","SQLE")


gene_Base2b_ConvsDL_targets = 
  gene_Base2b_ConvsDL_p30_CT[
    gene_Base2b_ConvsDL_p30_CT$dat.Name%in%intersect(gene_Base2b_ConvsDL,gene_Base2b_ConvsDL_same_direction$dat.Name),]


nrow(gene_Base2b_ConvsDL_targets)

write.csv(gene_Base2b_ConvsDL_targets,file = "gene_Base2b_ConvsDL_targets_p10p30.csv")


LUAD_1_target = LUAD_1[LUAD_1$rownames.gene_name_exp_dif.%in%gene_Base2b_ConvsDL_targets$dat.Name,]
nrow(LUAD_1_target)
gene_Base2b_ConvsDL_targets_1 = gene_Base2b_ConvsDL_targets[
  gene_Base2b_ConvsDL_targets$dat.Name%in%LUAD_1_target$rownames.gene_name_exp_dif.,]
nrow(gene_Base2b_ConvsDL_targets_1)
nrow(gene_Base2b_ConvsDL_targets[
  gene_Base2b_ConvsDL_targets$dat.Name%in%LUAD_1_target$rownames.gene_name_exp_dif.,])

LUAD_1_target[order(LUAD_1_target$rownames.gene_name_exp_dif.),]
gene_Base2b_ConvsDL_targets_1[order(gene_Base2b_ConvsDL_targets_1$dat.Name),]
?merge
m1 = merge(LUAD_1_target,gene_Base2b_ConvsDL_targets_1,by.x ="rownames.gene_name_exp_dif.",by.y = "dat.Name")

write.csv(m1[m1$logFC.x*m1$logFC.y<0,],"gene_Base2b_ConvsDL_targets_p10p30_cancer.csv")


################################################generation###################################


# P10 vs P0 CT 

gene_Base2b_P10vsP0_CT = read.csv("gene_Base2b_P10vsP0_CT.csv")


# P30 vs P0 CT 

gene_Base2b_P30vsP0_CT = read.csv("gene_Base2b_P30vsP0_CT.csv")

# P30 vs P10 CT 

gene_Base2b_P10vspP30_CT = read.csv("gene_Base2b_P10vspP30_CT.csv")
?intersect

generation = 
intersect(
intersect(
gene_Base2b_P10vsP0_CT[gene_Base2b_P10vsP0_CT$PValue<0.05,]$dat.Name,
gene_Base2b_P30vsP0_CT[gene_Base2b_P30vsP0_CT$PValue<0.05,]$dat.Name),
gene_Base2b_P10vspP30_CT[gene_Base2b_P10vspP30_CT$PValue<0.05,]$dat.Name
)

setdiff(m1[m1$logFC.x*m1$logFC.y<0,]$rownames.gene_name_exp_dif.,generation)



