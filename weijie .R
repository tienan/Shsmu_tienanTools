#weijie 
dat = read.table("../weijie/lncRNA/Upload/04.Quantification/data/all/all.gene.FPKM.mRNA.xls",header = T,sep = "\t")
dat = read.table("../weijie/lncRNA/Upload/04.Quantification/data/all/all.gene.FPKM.known_lncRNA.xls",header = T,sep = "\t")

head(dat)
colnames(dat) = gsub(x = colnames(dat),pattern = ".FPKM",replacement = "")
dat[c(1:3),1]
group = c(1,1,2,2,1,2,2,1,1,2)

#1975 PC vs TC
t1 = edgeR::DGEList(dat[,c(2:11)],group = as.factor(group))
t2 = edgeR::estimateCommonDisp(t1)
t3 = edgeR::exactTest(t2)
tableDEA <- edgeR::topTags(t3, n = nrow(t3$table))$table
tableDEA <- tableDEA[tableDEA$FDR <= 0.1, ]
nrow(tableDEA)

#foldChange>0, TC higher expression 
tableDEA$id = dat[rownames(tableDEA),1]


datTarget = dat[dat$GeneID%in%tableDEA$id,]
dat[9638,]
#biocLite("pheatmap")
library(pheatmap)
?pheatmap
tiff(filename = "thy_fig_1.tif",
     width = 12000, 
     height = 8000, 
     compression = "lzw",
     bg = "white", res = 800
)
a = pheatmap(
  datTarget[,-1],
  # clustering_distance_cols = "", 
  clustering_distance_rows = "euclidean",units = "px", pointsize = 1,
  cluster_rows = T,
  cluster_cols = T,
  scale="row",
  fontsize_row = 0.1, 
  fontsize_col = 14,
  show_rownames = T,
  show_colnames = T,
  angle_col = 45
  
)
dev.off()






hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mygenes =dat[rownames(tableDEA),1]
#write.csv(file = "tmp.txt",x = mygenes,row.names = F,col.names = F,quote = F)
getwd()
#attributes( hsmart )

mapping <- getBM(
  attributes = c('entrezgene_id', 'hgnc_symbol'), 
  filters = 'entrezgene_id',
  values = mygenes,
  mart =ensembl 
)
mapping = as.data.frame(mapping)

tmp = merge( tableDEA  ,mapping,by.x = "id",by.y= "entrezgene_id" )
tmp = na.omit(tmp)
write.csv(tmp,"DEGsPC&TC.csv")
ranks = tmp$logFC
names(ranks)=rownames(tmp)
ranks =sort(ranks,decreasing = T)

?fgsea
fgseaRes <- fgsea(pathwaysH, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaRes
plot(fgseaRes)
write.csv(fgseaRes[,c(1:5)],file = "fgseaResPC&TC.csv")


library(clusterProfiler)
kk <- enrichKEGG(gene = names(ranks), organism = 'hsa')
pathway  =  kk@result[kk@result$pvalue<0.1,]

write.csv(pathway,file = "KGGGPC&TC.csv")

