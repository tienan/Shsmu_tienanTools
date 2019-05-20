#setwd("C:\Users\tienan\Documents\R\")
# zhangjingfu Tongji 
setwd("C:/Users/tienan/Documents/R/")
dat = read.table("Human_jingfuzhang_1.txt",header = T,sep = "\t")
head(dat)
dat = read.table("Human_jingfuzhang_1.txt",header = T,sep = "\t")
id = read.table("tongren_zhangjingfu_id_name.txt",header = F,sep = "\t")
id = c(as.character(id[,1]),"TP53BP1","NUDT16L1","KAT5","RNF168")
dat$Gene=toupper(dat$Gene_id)

dat_1 =dat[dat$Gene%in%id,c(2:ncol(dat))]
rownames(dat_1) = dat[dat$Gene%in%id,1]


library(pheatmap)
tiff(filename = "Human_tongren_jingfu.tif",
     width = 3800, height = 4000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400
)
#?pheatmap
pheatmap(log(dat_1+0.1),
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         cluster_cols = T,
         #scale="row",
         fontsize_row = 10, 
         fontsize_col = 10,
         show_colnames = F
         )
dev.off()



order()
?merge

dat_1 = merge(dat[dat$Gene_id%in%id$V1,],id,by.x = "Gene_id",by.y = "V1")

dat_1=dat_1[order(dat_1$V2),]

rownames(dat_1)=dat_1$Gene_id

library(pheatmap)
tiff(filename = "?.tif",
     width = 1800, height = 4000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400
)
pheatmap(dat_1[,c(2:9)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         scale="row",
         fontsize_row = 10, 
         fontsize_col = 15)
dev.off()


dat[toupper(dat$Gene_id)%in%c("SIRT1","SIRT6"),]


source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
library(limma)


t1 = edgeR::DGEList(dat[,2:9], group = as.factor(c("Con","Con","Radiation",
                                                  "Radiation","Radiation","Radiation","Radiation","Radiation")))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)

diff_genes = t3$table[t3$table$PValue<0.01&abs(t3$table$logFC)>1,]

diff_genes = t3$table[t3$table$PValue<0.05,]

diff_genes[rownames(diff_genes)%in%c("16392","16397"),]

t3$table[rownames(t3$table)%in%c("16392","16397"),]

?kmeans
?scale
cl = kmeans(t(scale(t(dat_1[,c(2:9)]))),6)
cl$cluster


pheatmap(t(scale(t(dat_1[,c(2:9)]))))

cbind(dat_1[cl$cluster=="1",]$Gene_id,as.character(dat_1[cl$cluster=="1",]$V2))


cbind(dat_1[cl$cluster%in%c("2","3","4","5"),]$Gene_id,
      as.character( dat_1[cl$cluster%in%c("2","3","4","5"),]$V2))


dat_1$V2

ifelse()


dat = read.table("../tongren-mouse-data.txt",sep = "\t",header = T)

head(dat)
conditon = c("ctrl","CCL4","T1","T2","T3","T4")
# #样本是1-PBS组 2-CCL4造模组 
# 3-巨噬细胞exosome组  T1
# 4-巨噬细胞exosome包药物组 T2
# 5-药物组 T3
# 6-巨噬细胞组  T4
# 3-6组别都是治疗ccl4造肝纤维化 
# 



degs = read.table("../tongren-mouse-DEGs.txt",sep = "\t",header = T)
head(degs)
levels(degs$condition)
#("2vs1" "3vs2" "4vs2" "5vs2" "6vs2")
# 2vs1-> disease relative genes 
a0 = degs[degs$condition=="2vs1",]$Symbol
# 2vs1vs3 -> disease relative genes affected by T1
a1 = intersect(degs[degs$condition=="2vs1",]$Symbol,degs[degs$condition=="3vs2",]$Symbol)
# 2vs1vs4 -> disease relative genes affected by T2
a2 = intersect(degs[degs$condition=="2vs1",]$Symbol,degs[degs$condition=="4vs2",]$Symbol)
# 2vs1vs5-> disease relative genes affected by T3
a3 = intersect(degs[degs$condition=="2vs1",]$Symbol,degs[degs$condition=="5vs2",]$Symbol)
# 2vs1vs6-> disease relative genes affected by T4
a4 = intersect(degs[degs$condition=="2vs1",]$Symbol,degs[degs$condition=="6vs2",]$Symbol)




target_gene = read.table("../tongren_mouse_target_genes.txt",header = F,sep = "\t")
as.character(target_gene$V1)

dat[dat$Symbol%in%as.character(target_gene$V1),c(2:7)]

dat[dat$Symbol%in%union(union(union(a1,a2),a3),a4),c(2:7)]

target = dat[dat$Symbol%in%a0,]
rownames(target)=="49069"
target_1 = target[!rownames(target)=="49069",]
rownames(target_1) = target_1[,1]

rownames(target_1) = target_1[,1]
colnames(target_1) = c("Symbol","Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")


library(pheatmap)


tiff(filename = "tongren-mouse-figure-1.tif",
     width = 1800, height = 4000, units = "px", pointsize = 12,
     compression = "lzw",
     bg = "white", res = 400
)
pheatmap(target_1[,c(2:7)],
         #         clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         scale="row",
         fontsize_row = 6, 
         fontsize_col = 15)
dev.off()



target_gene[1,1]

dat[dat$Symbol=="Ap1ar",]


print(as.data.frame(union(union(union(a1,a2),a3),a4))[,1])

union(union(union(a1,a2),a3),a4)

i1 = intersect(degs[degs$condition=="3vs2",]$Symbol,degs[degs$condition=="4vs2",]$Symbol)
i1 = degs[degs$condition=="4vs2",]$Symbol
i2 = intersect(degs[degs$condition=="5vs2",]$Symbol,i1)
intersect(degs[degs$condition=="6vs2",]$Symbol,i2)




