source("http://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("pachterlab/sleuth")
library('sleuth')
s2c = read.table("slueths2c.txt",header = T,sep = "\t")

s2c$sample=as.factor(s2c$sample)
s2c$condition=as.factor(s2c$condition)
s2c$path=as.character(s2c$path)
so <- sleuth_prep(s2c,~condition)
source("https://bioconductor.org/biocLite.R")