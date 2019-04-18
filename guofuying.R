dat = read.table(file = "../R/guofuyin.txt",header = T,sep = "\t")
head(dat)
row.names(dat)=dat[,1]

library(pheatmap)

pheatmap(t(dat[,c(2,3)]),
         # clustering_distance_cols = "", 
         clustering_distance_rows = "euclidean",
         cluster_rows = T,
         #scale="row",
         fontsize_row = 15, 
         fontsize_col = 15)


pheatmap(tmp_file[,c(2:13)])

pheatmap(tmp_file[,c(2:13)],cluster_rows = F)

