#YugangWang, 
getwd()
dir()

#setwd("C:/Users/tienan/Documents/R")
dat = read.table("../yugangWangdata.txt",sep = "\t",header = T)
head(dat)

colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(4:5,10:11)],group = as.factor(sign[c(4:5,10:11)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(4:5,10:11)],as.numeric(dat[20609,c(4:5,10:11)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(4:5)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(4:5)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

ControlDiseaseModel = dat[rownames(tableDEA),c(4:5,10:11),]

colnames(ControlDiseaseModel)=c("Normal_1","Normal_2","DiseaseModel_1","DiseaseModel_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  ControlDiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
  )
dev.off()
write.csv(file = "Figure1aData",x = cbind(ControlDiseaseModel,dat[rownames(ControlDiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(ControlDiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "ControlDiseaseModel")






colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(10:11,2:3)],group = as.factor(sign[c(10:11,2:3)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(10:11,2:3)],as.numeric(dat[3861,c(10:11,2:3)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(2:3)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(2:3)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

CellDiseaseModel = dat[rownames(tableDEA),c(10:11,2:3)]

colnames(CellDiseaseModel)=c("DiseaseModel_1","DiseaseModel_2","Cell_1","Cell_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  CellDiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()
write.csv(file = "Figure1bData",x = cbind(CellDiseaseModel,dat[rownames(CellDiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(CellDiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "CellDiseaseModel")


#####################################################

colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(10:11,2:3)],group = as.factor(sign[c(10:11,2:3)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(10:11,2:3)],as.numeric(dat[3861,c(10:11,2:3)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(2:3)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(2:3)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

CellDiseaseModel = dat[rownames(tableDEA),c(10:11,2:3)]

colnames(CellDiseaseModel)=c("DiseaseModel_1","DiseaseModel_2","Cell_1","Cell_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  CellDiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()
write.csv(file = "Figure1bData",x = cbind(CellDiseaseModel,dat[rownames(CellDiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(CellDiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "CellDiseaseModel")

#############################################################################
colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(10:11,6:7)],group = as.factor(sign[c(10:11,6:7)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(10:11,6:7)],as.numeric(dat[51971,c(10:11,6:7)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(6:7)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(6:7)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

EPFDDiseaseModel = dat[rownames(tableDEA),c(10:11,6:7)]

colnames(EPFDDiseaseModel)=c("DiseaseModel_1","DiseaseModel_2","EPFD_1","EPFD_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  EPFDDiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()
write.csv(file = "Figure1cData",x = cbind(EPFDDiseaseModel,dat[rownames(EPFDDiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(EPFDDiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "EPFDDiseaseModel")


######################################################################

colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(10:11,8:9)],group = as.factor(sign[c(10:11,8:9)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(10:11,8:9)],as.numeric(dat[10582,c(10:11,8:9)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(8:9)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(8:9)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

EXODiseaseModel = dat[rownames(tableDEA),c(10:11,8:9)]

colnames(EXODiseaseModel)=c("DiseaseModel_1","DiseaseModel_2","EXO_1","EXO_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  EXODiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()
write.csv(file = "Figure1dData",x = cbind(EXODiseaseModel,dat[rownames(EXODiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(EXODiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "EXODiseaseModel")



######################################################################

colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

t1 = edgeR::DGEList(dat[,c(10:11,12:13)],group = as.factor(sign[c(10:11,12:13)]))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]

#Check
plot(sign[c(10:11,12:13)],as.numeric(dat[10582,c(10:11,12:13)]))

sign_2=c()
for (i in 1:nrow(tableDEA)){
  if ( tableDEA$logFC[i]<0){
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(12:13)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  } else {
    if(min(as.numeric(dat[rownames(tableDEA)[i],c(12:13)]))- 
       max(as.numeric(dat[rownames(tableDEA)[i],c(10:11)]))>0)
      sign_2[i] = 1
    else sign_2[i] = 0
  }
}

tableDEA=tableDEA[sign_2==1,]

PFDDiseaseModel = dat[rownames(tableDEA),c(10:11,12:13)]

colnames(PFDDiseaseModel)=c("DiseaseModel_1","DiseaseModel_2","FPD_1","FPD_2")


biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  PFDDiseaseModel,
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 15,
  show_rownames = F,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()
write.csv(file = "Figure1eData",x = cbind(PFDDiseaseModel,dat[rownames(PFDDiseaseModel),1]))

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[rownames(PFDDiseaseModel),1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)


?enrichMKEGG
kk <- enrichKEGG(gene = mapping [,3], organism = 'mmu')

keGG =  kk@result[kk@result$pvalue<0.05,]

write.csv(x = keGG,file = "PFDDiseaseModel")


