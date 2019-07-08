#YugangWang, 
getwd()
dir()

setwd("C:/Users/tienan/Documents/R")
dat = read.table("yugangWangdata.txt",sep = "\t",header = T)
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
  dat[rownames(tableDEA),c(4:5,10:11),],
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 20,
  show_rownames = F,
  show_colnames = T
  )
dev.off()





