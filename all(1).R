
bio
getwd()
setwd("~/DL/")
setwd("../R/DL/")
dat_Third = read.csv("ExpressDLThird.csv")
head(dat_Third)
dat_3 = dat_Third[,-1]
rownames(dat_3) = dat_Third[,1]
group = as.factor(gsub(pattern = "[0-9]",replacement = "",colnames(dat_3) ))
library(edgeR)
y = DGEList(dat_3,group = group,genes = dat_3)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="ENTREZID",column="SYMBOL")
y<-y[!is.na(y$genes$Symbol), ]
library(dplyr)
y$genes<- y$genes%>%
  group_by(Symbol) %>% 
  summarise_all(max)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
pch<- c(1:24)
colors<- rep(c("darkgreen","red","blue"),8)
plotMDS(y,col=colors[group],pch=pch[group])
legend("topright",legend=levels(group),pch=pch,col=colors,ncol=2)


design<-model.matrix(~0+group)
colnames(design)<-levels(group)

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)

#limma
design <- model.matrix(~ group)
fit <- lmFit(dat_3, design)
fit <- eBayes(fit)
contrast.matrix<-makeContrasts(groupHL-groupHS,levels=design)
fitC <- contrasts.fit(fit, contrast.matrix)
fitC <- eBayes(fitC)
go.fisher <- goana(fitC, species="Hs")
k <- kegga(fitC, species.KEGG="hsa") # equivalent to previous
#

###################################### H1975
contrast.matrix<-makeContrasts(HL-HS,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
is.de<-decideTestsDGE(res)
summary(is.de)
tr <- glmTreat(fit, contrast=contrast.matrix, lfc=log2(1.5))
topTags(tr)
is.de<-decideTestsDGE(tr)
summary(is.de)
plotMD(tr,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
tr$table$Symbol=mapIds(org.Hs.eg.db,rownames(tr),keytype="ENTREZID",column="SYMBOL")
DegH1975=tr$table[as.logical(abs(is.de)),]
write.csv(DegH1975,"~/lqR/DL/第二次/fc=1.5/DegH1975.csv")

############################################ A549
contrast.matrix<-makeContrasts(AL-AS,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)
tr <- glmTreat(fit, contrast=contrast.matrix, lfc=log2(1.5))
topTags(tr)
is.de<-decideTestsDGE(tr)
summary(is.de)
plotMD(tr,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
tr$table$Symbol=mapIds(org.Hs.eg.db,rownames(tr),keytype="ENTREZID",column="SYMBOL")
DegA549=tr$table[as.logical(abs(is.de)),]
write.csv(DegA549,"~/lqR/DL/第二次/fc=1.5/DegA549.csv")

################################# 取交集
IntDegUp=intersect(DegA549[DegA549$logFC>0,]$Symbol,DegH1975[DegH1975$logFC>0,]$Symbol)
write.csv(IntDegUp,"~/lqR/DL/第二次/fc=1.5/IntDegUp.csv")
IntDegDown=intersect(DegA549[DegA549$logFC<0,]$Symbol,DegH1975[DegH1975$logFC<0,]$Symbol)
write.csv(IntDegDown,"~/lqR/DL/第二次/fc=1.5/IntDegDown.csv")

################################# 取并集
UniDegUp=union(DegA549[DegA549$logFC>0,]$Symbol,DegH1975[DegH1975$logFC>0,]$Symbol)
write.csv(UniDegUp,"~/lqR/DL/第二次/fc=1.5/UniDegUp.csv")
UniDegDown=union(DegA549[DegA549$logFC<0,]$Symbol,DegH1975[DegH1975$logFC<0,]$Symbol)
write.csv(UniDegDown,"~/lqR/DL/第二次/fc=1.5/UniDegDown.csv")

################################ 交集基因集GO\KEGG分析
library(clusterProfiler)
genelist=mapIds(org.Hs.eg.db,IntDegUp,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/GoUp.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/KeggUp.csv")

genelist=mapIds(org.Hs.eg.db,IntDegDown,keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/GoDown.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/KeggDown.csv")

################################ GSE31210
install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE31210", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
fdata<-fData(gset[[1]])
############# 上调3year
clinData=pData[,c(2,52)]
clinData=clinData[!is.na(clinData$`days before death/censor:ch1`),]
clinData$group=ifelse(as.numeric(clinData$`days before death/censor:ch1`)/365>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegUp,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.symbol
exprSet_target=t(exprSet_target)
exprSet_target=data.frame(exprSet_target)
exprSet_target$geo_accession=rownames(exprSet_target)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[829,]
head(ExpCli)

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,829,828),])),nrow=nrow(ExpCli[-c(1,829,828),]), dimnames = list(rownames(ExpCli[-c(1,829,828),]),colnames(ExpCli[-c(1,829,828),])))


################limma
group1=as.factor(colnames(ExpCli1))
group1=as.factor(grou)

design <- model.matrix(~ group1-1)
head(design)
fit <- lmFit(ExpCli1, design)
fit <- eBayes(fit)
contrast.matrix<-makeContrasts(group1alive - group1dead,levels=design)
fitC <- contrasts.fit(fit, contrast.matrix)
fitC <- eBayes(fitC)
tmp = topTable(fitC,number = nrow(fitC))
nrow(tmp[tmp$P.Value<0.05,])
#################################
##########################
t1 = DGEList(ExpCli1,group = grou ,genes = ExpCli1)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]
head(tableDEA)
tableDEA$geneName = dat[rownames(tableDEA),1]
is.de<-decideTestsDGE(t3)
summary(is.de)
go.fisher <- goana(t3, species="Hs")
k <- kegga(t3, species="Hs")
###################

# 
# grou = as.factor(colnames(ExpCli1))
# colnames(ExpCli1) = c(1:ncol(ExpCli1))
y = DGEList(ExpCli1,group = grou ,genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+grou)
colnames(design)<-levels(grou)




y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
nrow(tmp1$table[tmp1$table$logFC>0.5,]) 
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
Up3yGSE31210=res$table[as.logical(abs(is.de)),]
write.csv(Up3yGSE31210,"~/lqR/DL/第二次/fc=1.5/GSE31210/Up3yGSE31210.csv")


###### 下调3year
clinData=pData[,c(2,52)]
clinData=clinData[!is.na(clinData$`days before death/censor:ch1`),]
clinData$group=ifelse(as.numeric(clinData$`days before death/censor:ch1`)/365>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegDown,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.symbol
exprSet_target=t(exprSet_target)
exprSet_target=data.frame(exprSet_target)
exprSet_target$geo_accession=rownames(exprSet_target)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[802,]

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,801,802),])),nrow=nrow(ExpCli[-c(1,801,802),]), dimnames = list(rownames(ExpCli[-c(1,801,802),]),colnames(ExpCli[-c(1,801,802),])))
y = DGEList(ExpCli1,group = as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
Down3yGSE31210=res$table[as.logical(abs(is.de)),]
write.csv(Down3yGSE31210,"~/lqR/DL/第二次/fc=1.5/GSE31210/Down3yGSE31210.csv")


################################ GO\KEGG分析
library(clusterProfiler)
genelist=mapIds(org.Hs.eg.db,rownames(Up3yGSE31210[Up3yGSE31210$logFC>0,]),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/GSE31210/GoUp.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/GSE31210/KeggUp.csv")

genelist=mapIds(org.Hs.eg.db,rownames(Down3yGSE31210[Down3yGSE31210$logFC<0,]),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/GSE31210/GoDown.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/GSE31210/KeggDown.csv")




################################## GSE37745
library(GEOquery)
gset <- getGEO("GSE37745", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
fdata<-fData(gset[[1]])

###### 上调3year
clinData=pData[,c(2,43,47)]
clinData = clinData[clinData$`histology:ch1`=="adeno",] 
clinData$group=ifelse(as.numeric(clinData$`days to determined death status:ch1`)/365>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegUp,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.ID
exprSet_target=t(exprSet_target)
exprSet_target=exprSet_target[-1,]
exprSet_target=as.matrix(exprSet_target)
geo_accession=rownames(exprSet_target)
exprSet_target=cbind(exprSet_target,geo_accession)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[830,]

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,828:830),])),nrow=nrow(ExpCli[-c(1,828:830),]), dimnames = list(rownames(ExpCli[-c(1,828:830),]),colnames(ExpCli[-c(1,828:830),])))
y = DGEList(ExpCli1,group = as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)#########无结果



####### 下调3year
clinData=pData[,c(2,43,47)]
clinData = clinData[clinData$`histology:ch1`=="adeno",] 
clinData$group=ifelse(as.numeric(clinData$`days to determined death status:ch1`)/365>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegDown,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.ID
exprSet_target=t(exprSet_target)
exprSet_target=exprSet_target[-2,]
exprSet_target=as.matrix(exprSet_target)
geo_accession=rownames(exprSet_target)
exprSet_target=cbind(exprSet_target,geo_accession)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[803,]

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,801:803),])),nrow=nrow(ExpCli[-c(1,801:803),]), dimnames = list(rownames(ExpCli[-c(1,801:803),]),colnames(ExpCli[-c(1,801:803),])))
y = DGEList(ExpCli1,group = as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)#########无结果

################################### GSE42127
library(GEOquery)
gset <- getGEO("GSE42127", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
fdata<-fData(gset[[1]])

###### 上调3year
clinData=pData[,c(2,45,46)]
clinData = clinData[clinData$`histology:ch1`=="Adenocarcinoma",] 
clinData$group=ifelse(as.numeric(clinData$`overall survival months:ch1`)/12>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegUp,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.ID
exprSet_target=t(exprSet_target)
exprSet_target=exprSet_target[-2,]
exprSet_target=data.frame(exprSet_target)
exprSet_target$geo_accession=rownames(exprSet_target)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[795,]

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,793:795),])),nrow=nrow(ExpCli[-c(1,793:795),]), dimnames = list(rownames(ExpCli[-c(1,793:795),]),colnames(ExpCli[-c(1,793:795),])))
y = DGEList(ExpCli1,group = as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)#########无结果


####### 下调3year
clinData=pData[,c(2,45,46)]
clinData = clinData[clinData$`histology:ch1`=="Adenocarcinoma",] 
clinData$group=ifelse(as.numeric(clinData$`overall survival months:ch1`)/12>3,"alive","dead")
exprSet=data.frame(exprSet)
exprSet$ID=rownames(exprSet)
fdata_target = fdata[fdata[,3]%in%UniDegDown,c(1,3,4)]
exprSet_target =merge(exprSet,fdata_target,by.x="ID",by.y="ID")
exprSet_target<-exprSet_target[!is.na(exprSet_target$`Gene symbol`), ]
exprSet_target<- exprSet_target%>%
  group_by(exprSet_target$`Gene symbol`) %>% 
  summarise_all(max)
exprSet_target=data.frame(exprSet_target)
rownames(exprSet_target) =exprSet_target$Gene.ID
exprSet_target=t(exprSet_target)
exprSet_target=exprSet_target[-2,]
exprSet_target=data.frame(exprSet_target)
exprSet_target$geo_accession=rownames(exprSet_target)
ExpCli=merge(exprSet_target,clinData,by.x="geo_accession",by.y="geo_accession")
ExpCli=t(ExpCli)
colnames(ExpCli)=ExpCli[802,]

ExpCli1=matrix(as.numeric(unlist(ExpCli[-c(1,800:802),])),nrow=nrow(ExpCli[-c(1,800:802),]), dimnames = list(rownames(ExpCli[-c(1,800:802),]),colnames(ExpCli[-c(1,800:802),])))
y = DGEList(ExpCli1,group = as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)#########无结果



#################### TCGA
library(SummarizedExperiment)
library(TCGAbiolinks)
library(limma)
query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")
LUADRnaseqSE <- GDCdownload(query.exp.hg38)
LUADRnaseqSE <- GDCprepare(query.exp.hg38)
rownames(LUADRnaseqSE) <- values(LUADRnaseqSE)$external_gene_name
exp.hg38.values <- assay(LUADRnaseqSE)
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

############# 上调3year
clinical_LUAD_m1$group=ifelse(clinical_LUAD_m1$survial_day/365>3,"alive","dead")
clinical_LUAD_m1=clinical_LUAD_m1[!is.na(clinical_LUAD_m1$group), ]
TargetExp = exp.hg38.values[rownames(exp.hg38.values)%in%UniDegUp,]
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(TargetExp) = patient_id
ExpCancer = TargetExp[,sign == "01A"]

ExpCancer = t(ExpCancer)
submitter_id=rownames(ExpCancer)

ExpCancer=cbind(submitter_id,ExpCancer)

ExpClin =merge(ExpCancer,clinical_LUAD_m1,by.x="submitter_id",by.y="clinical_LUAD.submitter_id")
ExpClin = t(ExpClin)

colnames(ExpClin)=ExpClin[917,]

ExpCli1=matrix(as.numeric(unlist(ExpClin[-c(1,809:917),])),nrow=nrow(ExpClin[-c(1,809:917),]), dimnames = list(rownames(ExpClin[-c(1,809:917),]),colnames(ExpClin[-c(1,809:917),])))

y = DGEList(ExpCli1,group =as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
Up3yTCGA=res$table[as.logical(abs(is.de)),]
write.csv(Up3yTCGA,"~/lqR/DL/第二次/fc=1.5/TCGA/Up3yTCGA.csv")


############ 下调3年

TargetExp = exp.hg38.values[rownames(exp.hg38.values)%in%UniDegDown,]
sign =c()
patient_id = colnames(exp.hg38.values)
for (i in 1:length(colnames(exp.hg38.values))){
  tmp = (strsplit2(as.character(patient_id[i]),split = "-"))
  sign[i] = tmp[4]
  tmp = paste(tmp[1],tmp[2],tmp[3],sep = "-")
  patient_id[i] = tmp
}
colnames(TargetExp) = patient_id
ExpCancer = TargetExp[,sign == "01A"]

ExpCancer = t(ExpCancer)
submitter_id=rownames(ExpCancer)

ExpCancer=cbind(submitter_id,ExpCancer)

ExpClin =merge(ExpCancer,clinical_LUAD_m1,by.x="submitter_id",by.y="clinical_LUAD.submitter_id")
ExpClin = t(ExpClin)

colnames(ExpClin)=ExpClin[890,]

ExpCli1=matrix(as.numeric(unlist(ExpClin[-c(1,881:890),])),nrow=nrow(ExpClin[-c(1,881:890),]), dimnames = list(rownames(ExpClin[-c(1,881:890),]),colnames(ExpClin[-c(1,881:890),])))

y = DGEList(ExpCli1,group =as.factor(colnames(ExpCli1)),genes = ExpCli1)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
y=calcNormFactors(y)
design<-model.matrix(~0+as.factor(colnames(ExpCli1)))
colnames(design)<-levels(as.factor(colnames(ExpCli1)))

y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
contrast.matrix<-makeContrasts(alive-dead,levels=design)
res<-glmQLFTest(fit,contrast=contrast.matrix)
tmp1 = topTags(res,n = 500)
is.de<-decideTestsDGE(res)
summary(is.de)
plotMD(res,status=is.de,values=c(1,-1),col=c("red","blue"),legend="topright")
Down3yTCGA=res$table[as.logical(abs(is.de)),]
write.csv(Down3yTCGA,"~/lqR/DL/第二次/fc=1.5/TCGA/Down3yTCGA.csv")




################################ GO\KEGG分析
library(clusterProfiler)
genelist=mapIds(org.Hs.eg.db,rownames(Up3yTCGA[Up3yTCGA$logFC>0,]),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/TCGA/GoUp.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/TCGA/KeggUp.csv")

genelist=mapIds(org.Hs.eg.db,rownames(Down3yTCGA[Down3yTCGA$logFC<0,]),keytype="SYMBOL",column="ENTREZID")
go=enrichGO(gene=genelist,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID",pvalueCutoff = 0.05)
write.csv(go,"~/lqR/DL/第二次/fc=1.5/TCGA/GoDown.csv")
kegg=enrichKEGG(genelist,organism = "hsa",pvalueCutoff = 0.05)
write.csv(kegg,"~/lqR/DL/第二次/fc=1.5/TCGA.csv")

