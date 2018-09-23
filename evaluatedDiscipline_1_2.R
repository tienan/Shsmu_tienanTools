#
dat = read.table("evaluatedDiscipline_1_2.txt",head=F,sep = "\t")
sdResult = apply(dat, 2, sd)
weight = sdResult/sum(sdResult)

sumResult = (as.matrix(dat) %*% as.matrix(weight))
for (i in 1:nrow(sumResult)){
  cat(sumResult[i])
  cat("\n")
}
