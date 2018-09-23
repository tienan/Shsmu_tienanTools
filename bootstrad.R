dat = read.table("tmp",sep = "\t") # read data 
dat[,3] = dat[,2]-dat[,1]
IR <- function( data, indices) {
  d <- data[indices,3] # allows boot to select sample 
  InterquatileRange = fivenum(d)[4]-fivenum(d)[2]
  return(InterquatileRange)
} 
results <- boot(dat, IR, R=1000)
results
boot.ci(results, type="norm")


dat = read.table("tmp",sep = "\t")
dat[,3] = abs(dat[,2]-dat[,1])
IR <- function( data, indices) {
  d <- data[indices,3] # allows boot to select sample 
  InterquatileRange = abs(fivenum(d)[4]-fivenum(d)[2])
  return(InterquatileRange)
} 
results <- boot(dat, IR, R=1000)
results
boot.ci(results, type="bca")




dat=read.table("tmp",sep = "\t")
P30 = function( data, indices) {
  a = abs(data[indices,2]-data[indices,1])/data[indices,1]
  b = length(a[a<=0.3])/length(a)
  b
}
results <- boot(dat, P30, R=1000)
results
boot.ci(results, type="norm")



P10 = function( data, indices) {
  a = abs(data[indices,2]-data[indices,1])/data[indices,1]
  b = length(a[a<=0.1])/length(a)
  b
}
results <- boot(dat, P10, R=1000)
results
boot.ci(results, type="norm")

rmse = function( data, indices) {
  a = sqrt(mean(abs(data[indices,2]-data[indices,1])^2))
  a
}
results <- boot(dat, rmse, R=1000)
results
boot.ci(results, type="norm")



