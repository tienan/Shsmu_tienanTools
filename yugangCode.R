#YugangWang, 
library(edgeR)

getwd()
dir()

#setwd("C:/Users/tienan/Documents/R")

dat = read.table("/media/tienan/0000678400004823/R/yugangWangdata.txt",sep = "\t",header = T)
head(dat)

colnames(dat) = c("seq_id",4,4,1,1,6,6,3,3,2,2,5,5)
sign=colnames(dat)
#"seq_id" "CELL1"  "CELL2"  "CN1"    "CN2"    "EPFD1"  "EPFD2"  "EXO1"  
#"EXO2"   "M1"     "M2"     "PFD1"   "PFD2" 

library(edgeR)
dat_1 = dat[,-1]
rownames(dat_1) = dat[,1]
group = gsub(pattern = "[0-9]",replacement = "",colnames(dat_1))
head(dat)
dat$Symbol<-mapIds(org.Mm.eg.db,keys = as.character(dat[,1]),keytype="ENSEMBL",column="SYMBOL",multiVals = "first")

library(dplyr)




y = DGEList(dat_1,group = group,genes = dat_1)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
rownames(y)
?mapIds
y$genes$Symbol<-mapIds(org.Mm.eg.db,rownames(y),keytype="ENSEMBL",column="SYMBOL",multiVals = "first")

rownames(y) <- mapIds(org.Mm.eg.db, keys = rownames(y), keytype = "ENSEMBL", column="ENTREZID",multiVals = "first")



require(hgu95av2.db)
columns(hgu95av2.db)
select(hgu95av2.db, keys=rownames(y), columns = c("SYMBOL","ENTREZID"))


head(y$genes)
y<-y[!is.na(y$genes$Symbol), ]
dim(y)
rownames(y)=y$genes$Symbol
y$genes$Symbol<-mapIds(org.Hs.eg.db,y$genes$Symbol,keytype="ENTREZID",column="SYMBOL")
#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#library(org.Mm.eg.db)

#y$genes$Symbol = mapIds(org.Hs.eg.db,rownames(y),keytype = "ENTREZID",column = "SYMBOL")
#group=factor(name)
#table(group)
#y = y[!is.na(y$genes$Symbol ),]
#dim(y)
#keep 
y=calcNormFactors(y)
y$samples

pch<- c(0,1,2,15)
colors<- rep(c("blue","red"),2)
?plotMDS
plotMDS(y,col=colors[group],pch=pch[group],labels=rownames(y$samples))
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)

head(dat_1)

dat_2 = dat_1[,c(-3,-9)]
head(dat_2)

y = DGEList(dat_2,group = group[c(-3,-9)],genes = dat_2)
library(org.Hs.eg.db)
rownames(y)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="SYMBOL",column="ENTREZID")
head(y$genes)
y<-y[!is.na(y$genes$Symbol), ]
dim(y)
rownames(y)=y$genes$Symbol
y$genes$Symbol<-mapIds(org.Hs.eg.db,y$genes$Symbol,keytype="ENTREZID",column="SYMBOL")
head(y$genes)




design<-model.matrix(~0+group[c(-3,-9)])
colnames(design)<-levels(group[c(-3,-9)])
design
y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
H1975dlcon<-makeContrasts(H1975_dl-H1975_con,levels=design)
res<-glmQLFTest(fit,contrast=dlcon)
tmp1 = topTags(res,n = 500)
#View(tmp1$table)
is.de<-decideTestsDGE(res,p.value = 0.2)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
go<-goana(res,species="Hs",FDR = 0.1)
go[go$P.Up<0.05&go$P.Down<0.05,]
# #Find the key words of mianyi, zhidaixie, chundaixie
# #grepl(pattern = )
# library(gplots)
# logCPM.1 = logCPM[as.logical(is.de@.Data),]
# View(is.de)
# 
# rownames(logCPM.1)<-logCPM.1
# colnames(logCPM)=colnames(y$genes)[1:10]
# logCPM<-t(scale(t(logCPM)))
# col.pan<-colorpanel(100,"blue","white","red")
# heatmap.2(logCPM[as.logical(is.de@.Data),],col=col.pan,Rowv=TRUE,scale="none",
#           trace="none",dendrogram="both",cexRow=1,cexCol=1.4,density.info="none",
#           margin=c(10,9),lhei=c(2,10),lwid=c(2,6))
keg<-kegga(res,species="Hs")
keg[keg$P.Up<0.1|keg$P.Down<0.1,]
load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
cam<-camera(y, idx, design,contrast=H1975dlcon,inter.gene.cor=0.01)
cam[cam$FDR<0.05,]

















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


#
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


??enrichMKEGG
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



install.packages("VennDiagram")
library(VennDiagram)

A = rownames(ControlDiseaseModel) 
B = rownames(CellDiseaseModel)
C = rownames(EPFDDiseaseModel) 
D = rownames(EXODiseaseModel)
E = rownames(PFDDiseaseModel) 

venn.diagram(list(ControlDiseaseModel=A,
                  CellDiseaseModel=B,EPFDDiseaseModel = C, 
                  EXODiseaseModel=D,PFDDiseaseModel=E ),
             height = 1550, width = 1700,
             compression = "lzw",
             resolution = 300, imagetype = "tiff", 
             #alpha=c(0.5,0.5,0.5),
             #fill=c("red","yellow","blue"), 
             cat.fontface=4,fontfamily=3,
             #main="Intersection of WD40 genes identified by different methods",
             main.cex = 2, main.fontface = 2, main.fontfamily = 3,
             cat.pos = c(-20, 0, -70,50,-20),
             filename = "VennDiagram.tif")
?venn.diagram


a = venn.diagram(list(ControlDiseaseModel=A,
                  CellDiseaseModel=B,EPFDDiseaseModel = C, 
                  EXODiseaseModel=D,PFDDiseaseModel=E ),
             height = 1550, width = 1700,
             compression = "lzw",
             resolution = 300, imagetype = "tiff", 
             #alpha=c(0.5,0.5,0.5),
             #fill=c("red","yellow","blue"), 
             cat.fontface=4,fontfamily=3,
             #main="Intersection of WD40 genes identified by different methods",
             main.cex = 2, main.fontface = 2, main.fontfamily = 3,
             cat.pos = c(-20, 0, -70,50,-20),
             filename = "VennDiagram.tif")

intersectsSet = intersect(intersect(intersect(intersect(A,B),C),D),E)
dataInsect = dat[intersectsSet,-1]
colnames(dataInsect)= c("CELL1","CELL2","CN1","CN2","EPFD1","EPFD2","EXO1","EXO2","M1","M2","PFD1","PFD2")
rownames(dataInsect) = dat[intersectsSet,1]

write.csv(file = "dataInsect",x=dataInsect)

pheatmap(
  as.data.frame(dataInsect),
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 10, 
  fontsize_col = 15,
  show_colnames = T,
  show_rownames = T,
  angle_col = 45
)

library(biomaRt)

musmart <- useMart(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")


# Object of class 'Mart':
#   Using the ENSEMBL_MART_ENSEMBL BioMart database
#   Using the hsapiens_gene_ensembl dataset
mygenes <- dat[intersectsSet,1]
mapping <- getBM(
  attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'entrezgene_id', 'hgnc_symbol'), 
  filters = 'ensembl_gene_id',
  values = mygenes,
  mart = musmart
)
mapping[,4]
