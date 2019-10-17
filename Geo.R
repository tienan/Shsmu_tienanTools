#Geo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE68072", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
fdat = fData(gset[[1]])

write.csv2()
read.csv2()
pData <- pData(gset[[1]])
sample <- pData$geo_accession
group = ifelse(grepl(pattern = "non",pData[,1]),1,0)
group <- rep(c(1,0),times=c(12,11))
design <- model.matrix(~group)
rownames(design)=colnames(exprSet)
#colnames(design)=levels(factor(group))
design=as.matrix(cbind(Grp1=1,Grp2vs1=group))
#差异基因分析
fit <- lmFit(exprSet,design)
fit1 <- eBayes(fit)
topTable(fit1,coef = "Grp2vs1",number = 33292,p.value = 0.05,lfc = 0.7,adjust="BH")



fit2 <- treat(fit,lfc=0.1)
topTable(fit)
topTable(fit, number = 33292,p.value = 0.05,lfc = 1.5,coef="ATHvsNATH", adjust="BH")
topTreat(fit2,coef = 2,number = 100)
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)
?options

fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("SummarizedExperiment")
biocLite("TCGAbiolinks")
biocLite(limma)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)




group = ifelse(grepl(pattern = "non",pData[,1]),"1","0")

t1 = edgeR::DGEList(exprSet,group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]

tableDEA <- tableDEA[tableDEA$PValue <= 0.1, ]






colnames(exprSet)
pData[,1]


require(graphics)

t.test(1:10, y = c(7:20))      # P = .00001855
t.test(1:10, y = c(7:20, 200)) # P = .1245    -- NOT significant anymore

## Classical example: Student's sleep data
plot(extra ~ group, data = sleep)
## Traditional interface
with(sleep, t.test(extra[group == 1], extra[group == 2]))
## Formula interface
t = t.test(extra ~ group, data = sleep)
t$p.value
