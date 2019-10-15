#Transcript_of_RNA_Signature_2

######################Setting the work direction
getwd()
setwd("../DL/")
######################Reading FPKM value
geneNmae = read.table("uniqGeneName",header = F)
file_list = dir(pattern = "*.fpkm*")
dat_tmp = read.table(file_list[1],header = T,sep = "\t")
tmp=c()
tmp_name=dat_tmp$tracking_id
for (i in 1:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%geneNmae$V1,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  tmp_name = intersect(tmp_name,dat_tmp$tracking_id)
  tmp[i] = nrow(dat_tmp)
  #head(dat_tmp)
  #tmp_file = cbind(tmp_file,dat_tmp[,2])
}
length(tmp_name)

tmp_file = read.table(file_list[1],header = T,sep = "\t")
#head(tmp_file)
#nrow(tmp_file)
#colnames(tmp_file)
tmp_file = tmp_file[,c(1,10)]
tmp_file = tmp_file[tmp_file$tracking_id%in%tmp_name,]
tmp_file = tmp_file[order(tmp_file$tracking_id),]
head(tmp_file)
library(dplyr)
tmp_file = summarise(group_by(tmp_file,tracking_id),maxFPKM=max(FPKM))
nrow(tmp_file)

for (i in 2:length(file_list)){
  dat_tmp = read.table(file_list[i],header = T,sep = "\t")
  dat_tmp = dat_tmp[,c(1,10)]
  dat_tmp = dat_tmp[dat_tmp $tracking_id%in%tmp_name,]
  dat_tmp = dat_tmp[order(dat_tmp$tracking_id),]
  dat_tmp = summarise(group_by(dat_tmp,tracking_id),maxFPKM=max(FPKM))
  #nrow(dat_tmp)
  #head(dat_tmp)
  tmp_file = cbind(tmp_file,dat_tmp[,2])
}
colname = gsub(pattern = "_genes.fpkm_tracking",replacement = "",file_list) 
colname = gsub(pattern = "1975",replacement = "H1975",colname) 
tmp = strsplit2(x = colname,split = "_")
group = as.factor(paste(tmp[,1],tmp[,2],sep="_"))
colnames(tmp_file) = c("Gene",colname)
condition = ifelse(grepl(pattern = "con",colname),"Con","DL")
library(edgeR)
dat_1 = tmp_file[,-1]
rownames(dat_1) = tmp_file[,1]

library(edgeR)

y = DGEList(dat_1,group = group,genes = dat_1)
library(org.Hs.eg.db)
rownames(y)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="SYMBOL",column="ENTREZID")


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
######################################################
######################################################

dat_Third = read.csv("ExpressDLThird.csv")
head(dat_Third)
dat_3 = dat_Third[,-1]
rownames(dat_3) = dat_Third[,1]
group = as.factor(gsub(pattern = "[0-9]",replacement = "",colnames(dat_3) ))
y = DGEList(dat_3,group = group,genes = dat_3)
library(org.Hs.eg.db)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="ENTREZID",column="SYMBOL")
head(y$genes)
y<-y[!is.na(y$genes$Symbol), ]
dim(y)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
pch<- c(1:24)
colors<- rep(c("darkgreen","red","blue"),8)
plotMDS(y,col=colors[group],pch=pch[group])
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)


design<-model.matrix(~0+group)
colnames(design)<-levels(group)
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




######################################################

res$table[res$table$PValue<0.05,]

tmp1 = topTags(res)
tmp1$table
boxplot(tmp1$table$PValue)
fivenum(tmp1$table$PValue)

tr<-glmTreat(fit,contrast=dlcon,lfc=log2(1.5))
tr$table[tr$table$PValue<0.05,]
is.de<-decideTestsDGE(tr)
summary(is.de)
is.de<-decideTestsDGE(res,lfc = 0.8)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")

# library(Hmisc)
# capitalize(tr$genes$Symbol)
# tr$genes$Symbol = capitalize(tr$genes$Symbol)
# tr$genes$Symbol = toupper(tr$genes$Symbol)
# 
# mapIds(org.Hs.eg.db,rownames(y),keytype="SYMBOL",column="ENTREZID")
# 
#

go<-goana(res,species="Hs",FDR = 0.1)
go
summary(res)
topTags(res)
res$table$PValue = res$table$PValue/7
keg<-kegga(res,species="Hs")
topKEGG(keg,n=15,truncate=34)
library(GO.db)
Rkeys(org.Hs.egGO2ALLEGS)

term = select(GO.db,keys=Rkeys(org.Hs.egGO2ALLEGS),columns="TERM")
head(term)  
term[grepl(pattern = "tnf",term$TERM,ignore.case = T),]


load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
idx<-ids2indices(Hs.c2,id=rownames(y))
dlvscom<-makeContrasts(H1975_dl-H1975_con,levels=design)
cam<-camera(y, idx, design,contrast=dlvscom,inter.gene.cor=0.01)
head(cam,14)

?decideTestsDGE
is.de<-
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")

topTags(res)
res$genes$Symbol = mapIds(org.Hs.eg.db,rownames(res),keytype = "ENTREZID",column = "SYMBOL")
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
logCPM<-cpm(y,prior.count=2,log=TRUE)
rownames(logCPM)<-y$genes$Symbol
colnames(logCPM)<-paste(y$samples$group,1:3,sep="-")
tr<-glmTreat(fit,contrast=ASAL,lfc=log2(1.5))
o<-order(tr$table$PValue)
nrow(tr)
is.de<-decideTestsDGE(tr)
summary(is.de)
tr
logCPM<-logCPM[is.de!=0,]
tr[is.de!=0,]$table

##############################GEO dataset GSE31210
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE31210", GSEMatrix =TRUE, AnnotGPL=TRUE )



exprSet <- exprs(gset[[1]])

#exprSet$id = rownames(exprSet)

pData <- pData(gset[[1]])

fdata<-fData(gset[[1]])

fdata_target = fdata[fdata[,3]%in%diff_gene_filer_1[,1],c(1,3)]

exprSet_target = as.data.frame(exprSet[rownames(exprSet)%in%fdata_target[,1],])
exprSet_target$ID = rownames(exprSet_target)

targetGene = merge(fdata_target,exprSet_target,by.x="ID",by.y="ID")
targetGene
targetGeneM = targetGene[,c(-1,-2)]

targetGeneM = (as.data.frame(t(targetGeneM)))
targetGeneM[c(1:10),c(1:10)]

head(targetGeneM[,c(1:10)])
colnames(targetGeneM) = targetGene[,2]
targetGeneM$ID = rownames(targetGeneM) 

colnames(pData)
clinData=pData[,c(52,54,63)]
clinData$ID = rownames(pData)
as.data.frame(clinData)


clin_DL = as.data.frame(merge(targetGeneM,clinData,by.x="ID",by.y="ID"))
head(clin_DL[,c(1:10)])
colnames(clin_DL)
ncol(clin_DL)
group = ifelse(as.numeric(clin_DL[,1037])/365>3,0,1)
tmp = colnames(clin_DL)
tmp[1036]
exp = as.data.frame(t(clin_DL[,c(2:1036)]))
exp_1 = as.data.frame((exp))

exp_1  = exp[,!is.na(group) ]
group_1 = group[!is.na(group) ]

t1 = edgeR::DGEList(exp_1,group = as.factor(group_1))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$PValue <= 0.01, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)

##############################GEO dataset GSE37745
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)

gset <- getGEO("GSE37745", GSEMatrix =TRUE, AnnotGPL=TRUE )

exprSet <- exprs(gset[[1]])

#exprSet$id = rownames(exprSet)

pData <- pData(gset[[1]])

fdata<-fData(gset[[1]])

fdata_target = fdata[fdata[,3]%in%diff_gene_filer_1[,1],c(1,3)]

exprSet_target = as.data.frame(exprSet[rownames(exprSet)%in%fdata_target[,1],])
exprSet_target$ID = rownames(exprSet_target)

targetGene = merge(fdata_target,exprSet_target,by.x="ID",by.y="ID")
targetGene
targetGeneM = targetGene[,c(-1,-2)]

targetGeneM = (as.data.frame(t(targetGeneM)))
targetGeneM[c(1:10),c(1:10)]

head(targetGeneM[,c(1:10)])
colnames(targetGeneM) = targetGene[,2]
targetGeneM$ID = rownames(targetGeneM) 

colnames(pData)
clinData=pData[,c(45,43,47,50)]
clinData$ID = rownames(pData)
as.data.frame(clinData)
clinData = clinData[clinData$`histology:ch1`=="adeno",] 


clin_DL = as.data.frame(merge(targetGeneM,clinData,by.x="ID",by.y="ID"))
head(clin_DL[,c(1:10)])
colnames(clin_DL)
ncol(clin_DL)
group = ifelse(as.numeric(clin_DL$`days to determined death status:ch1`)/365>5,0,1)
tmp = colnames(clin_DL)
tmp[1036]
exp = as.data.frame(t(clin_DL[,c(2:1036)]))
exp_1  = exp[,!is.na(group) ]
group_1 = group[!is.na(group) ]

t1 = edgeR::DGEList(exp_1,group = as.factor(group_1))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$PValue <= 0.01, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)


