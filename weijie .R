#weijie 
dat = read.table("/media/tienan/0000678400004823/R/weijie/lncRNA/Upload/04.Quantification/data/all/all.gene.FPKM.mRNA.xls",header = T,sep = "\t")
dat = read.table("../weijie/lncRNA/Upload/04.Quantification/data/all/all.gene.FPKM.known_lncRNA.xls",header = T,sep = "\t")
nrow(dat)
head(dat)
colnames(dat) = gsub(x = colnames(dat),pattern = ".FPKM",replacement = "")
dat[c(1:3),1]
group = c(1,1,2,2,1,2,2,1,1,2)
group = gsub(pattern = "[0-9]",replacement = "",colnames(dat))
group = as.factor(group[-1])

y = DGEList(dat[,-1],group = group,genes = dat[,-1])
library(org.Hs.eg.db)
y$genes$Symbol<-mapIds(org.Hs.eg.db,rownames(y),keytype="ENTREZID",column="SYMBOL")
head(y$genes)
y<-y[!is.na(y$genes$Symbol), ]
dim(y)
keep<-rowSums(cpm(y) > 0.5) >= 2
y<-y[keep, ,keep.lib.sizes=FALSE]
colors<- rep(c("darkgreen","red"),1)
plotMDS(y,col=colors[group],pch=pch[group])
legend("topleft",legend=levels(group),pch=pch,col=colors,ncol=2)


design<-model.matrix(~0+group)
colnames(design)<-levels(group)
design
y<-estimateDisp(y, design,robust=TRUE)
#install.packages("statmod")
fit<-glmQLFit(y, design,robust=TRUE)
#################################################PCvsTC
PCvsTC<-makeContrasts(PC-TC,levels=design)
res<-glmQLFTest(fit,contrast=PCvsTC)
res$table[abs(res$table$logFC)>1.5,]$PValue=1.750675e-05
is.de<-decideTestsDGE(res)
summary(is.de)
tmp = topTags(res,n = 3000)
tmp
degs = rownames(tmp)
go<-goana(res,species="Hs")

cyt.go = go[go$P.Up<0.05|go$P.Down<0.05,]
nrow(cyt.go)
library(GO.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
#############
con <-  file("PCvsTCGoGeneList.txt", open = "w")
for (i in 1:nrow(cyt.go)){
  #Rkeys(org.Mm.egGO2ALLEGS) = "GO:0032465"
  Rkeys(org.Hs.egGO2ALLEGS) = rownames(cyt.go)[i]
  cyt.go.genes = as.list(org.Hs.egGO2ALLEGS)
  # length(cyt.go.genes[[i]])
  # length(unique(cyt.go.genes[[i]]))
  tmpName = names(cyt.go.genes)
  goGeneList = intersect(degs,cyt.go.genes[[1]])
  goGeneListSym="NA"
  if(length(goGeneList)>0){
    goGeneListSym = mapIds(org.Hs.eg.db,keys =  goGeneList,keytype="ENTREZID",column="SYMBOL")
  }
  writeLines(text = tmpName, con = con )
  writeLines(text = goGeneListSym,sep = "\t", con = con )
  writeLines(text = "",con = con )
  rm(org.Hs.egGO2ALLEGS)
}
close(con)
library(KEGGREST)
keg<-kegga(res,species="Hs")
nrow(keg)
cyt.keg = keg[keg$P.Up<0.05|keg$P.Down<0.05,]
row.names(cyt.keg)
con <-  file("PCvsTCKeggGeneList.txt", open = "w")
for (i in 1:nrow(cyt.keg)){
  cyt.kegList = keggGet(row.names(cyt.keg)[i])
  tmpName = cyt.kegList[[1]]$NAME
  writeLines(text = tmpName, con = con )
  if(!is.null(cyt.kegList[[1]]$GENE)){
    kegGeneList = intersect(degs,cyt.kegList[[1]]$GENE[seq(1,length(cyt.kegList[[1]]$GENE),2)])
  }
  kegGeneListSym = mapIds(org.Hs.eg.db,keys =  kegGeneList,keytype="ENTREZID",column="SYMBOL")
  writeLines(text =   kegGeneListSym,sep = "\t", con = con )
  writeLines(text = "",con = con ) 
}
close(con)

load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
#load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata"))
idx<-ids2indices(Hs.c2,id=rownames(y))
cam<-camera(y, idx, design,contrast=PCvsTC,inter.gene.cor=0.4)
cam[cam$FDR<0.05,]
write.csv(x = tmp,file = "PCvsTC.csv")
write.csv(x = cyt.go,file = "PCvsTCGo.csv")
write.csv(x = cyt.keg ,file = "PCvsTCKegg.csv")
write.csv(x = cam[cam$FDR<0.05,],file = "PCvsTCGASE.csv")













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

