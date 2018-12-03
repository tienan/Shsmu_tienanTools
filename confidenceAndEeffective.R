dat = read.table("shaolei/ca.txt")
head(dat)
table(dat[,1])
table(dat[,2])
table(dat[,3])
table(dat[,4])

dat = read.table("shaolei/cs.txt")

table(dat[,1])


source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager")
