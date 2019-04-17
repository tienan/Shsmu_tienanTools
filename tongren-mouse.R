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




