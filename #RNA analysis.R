#RNA analysis

getwd()
setwd("../")

RNAexpDate = read.csv("RNAExpThirdC1comma.txt")
head(RNAexpDate)
NameDat = read.csv("trans2geneComma.txt",header = F)
head(NameDat)
head(NameDat[,2])

NameDat[981,2]



source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("SummarizedExperiment")
biocLite("TCGAbiolinks")
biocLite(limma)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(edgeR)

grepl()



group = ifelse(grepl(pattern = "T",colnames(RNAexpDate)[4:149]),"T","N")


ifelse(RNAexpDate[1,])

head(RNAexpDate[4:149])

RNAexpDateCalc = RNAexpDate[4:149]

rownames(RNAexpDateCalc)=RNAexpDate[,1]

t1 = edgeR::DGEList(RNAexpDateCalc,group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]


t3$table



DGE <- edgeR::DGEList(TOC, group = rep(c(Cond1type, 
                                         Cond2type), c(Cond1num, Cond2num)))
disp <- edgeR::estimateCommonDisp(DGE)
tested <- edgeR::exactTest(disp, pair = c(Cond1type, 
                                          Cond2type))
logFC_table <- tested$table
tableDEA <- edgeR::topTags(tested, n = nrow(tested$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= fdr.cut, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= logFC.cut, ]
tableDEA$name = rownames(tableDEA)

nrow(tableDEA[tableDEA$name%in%NameDat$V1,])

nrow(tableDEA)

colnames(NameDat)=c("name","gene")

merge(tableDEA, NameDat, by.x = "name", by,y = "name")


difGene = merge(tableDEA, NameDat, by.x = "name", by.y = "name")

head(difGene)

difGene[difGene$gene=="a",]$gene="LPAL2"


head(difGene)

library(dplyr)

?group_by

by_cyl <- mtcars %>% group_by(cyl)

by_cyl %>% filter(disp == max(disp))


by_gene <-difGene  %>% group_by(gene)

by_gene


write.csv(by_gene,file = "mRNADiff.csv")

write.csv(summarise(by_gene,logFC=max(logFC),logCPM=max(logCPM),FDR=min(FDR)),file = "gene.csv")


RNAvalue = RNAexpDateCalc[RNAexpDate$Name%in%by_gene$name,]

apply(RNAvalue, 1, mean)>1


biocLite("pheatmap")
library(pheatmap)
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
         RNAvalue,
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
         cluster_rows = T,
         cluster_cols = T,
         scale="row",
         fontsize_row = 0.1, 
         fontsize_col = 6)
dev.off()

sumpheat = summary(a)





############################################################# subgroup analysis
getwd()
# RNAexpDate = read.csv("RNAExpThirdC1comma.txt")
# head(RNAexpDate)
thyriodData = read.table("thyriodData.txt",header = T,sep = "\t")
head(thyriodData)
tail(thyriodData)
raw = thyriodData
rownames(raw)=raw[,1]
thyriodData=raw[,raw[nrow(raw),]!=2]
thyriodData=thyriodData[,-1]
t1 = edgeR::DGEList(thyriodData[c(1:nrow(thyriodData)-1),],group = ``)
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
logFC_table <- t3$table
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.05, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= 1, ]





tableDEA <- edgeR::topTags(tested, n = nrow(tested$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= fdr.cut, ]
tableDEA <- tableDEA[abs(tableDEA$logFC) >= logFC.cut, ]
tableDEA$name = rownames(tableDEA)


