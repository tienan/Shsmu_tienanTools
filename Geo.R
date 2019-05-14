#Geo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
library(GEOquery)
gset <- getGEO("GSE68072", GSEMatrix =TRUE, AnnotGPL=TRUE )
exprSet <- exprs(gset[[1]])
pData <- pData(gset[[1]])
sample <- pData$geo_accession
group <- rep(c(1,0),times=c(12,11))
design <- model.matrix(~group)
rownames(design)=colnames(exprSet)
colnames(design)=levels(factor(group))
#差异基因分析
fit <- lmFit(exprSet, design)
fit2 <- eBayes(fit)
topTable(fit2, number = 33297,p.value = 0.05,lfc = 1.5,coef="ATHvsNATH", adjust="BH")