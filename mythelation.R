if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("MEDIPS")
library("BSgenome")
available.genomes()
BSgenome.Rnorvegicus.UCSC.rn5

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")


devtools::install_github("andigoni/meinter")




###
getwd()
#Venn

S_1 <- read.table("../R/RadiationMythylation/S-1.bed")
S_1$region = paste(S_1$V2,S_1$V3,sep = "_")
head(S_1)
summary(S_1$V5)
T_1 = read.table("../R/RadiationMythylation/T-1.bed")
T_1$region = paste(T_1$V2,T_1$V3,sep = "_")
head(T_1)
summary(T_1$V5)

intersect(S_1$region,T_1$region)


e28 <- read.csv("e28_4g.dsg2.csv")
m17_e28 <- read.csv("m17_e28_deg.csv",header = TRUE)
m17 <- as.vector(unlist(m17[1]))
e28 <- as.vector(unlist(e28[1]))
m17_e28 <- as.vector(unlist(m17_e28[1]))

#主要是获取对应的差异表达基因，转换成向量。

#veen.diagram中变量名称不能使用-,否则会报错。变量名要求是字母开头。
venn.plot <- venn.diagram(
  x = list(
    Mo17_DY13 = m17,       #对应三个元的，前面是图中显示字符，后面是实际变量
    E28_DY13 = e28,
    Mo17_E28 = m17_e28
  ),
  filename = "DY13.tiff",  #输出的文件名
  col = "transparent",
  fill = c("red", "blue", "green"),  #各个圈圈的颜色
  alpha = 0.5,  #透明度
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",  #此处默认是标签在圆内部
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)


########mythylationKit 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("methylKit")

BiocManager::install("methylKit")
