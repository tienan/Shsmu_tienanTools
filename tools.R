#tools

getwd()
setwd("C:\\Users\\tienan\\Documents\\R")
data = read.table("thyroidCancerList_1.txt",header = F)
data
?matrix
data.frame(as.character(data), nrow = 4, ncol = 73, byrow = TRUE)
