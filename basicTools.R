writeOut = function(dat,name){
  write.table(dat,file=name,quote = F,sep="\t",row.names = T,fileEncoding = "utf-8")
}

writeOutNoRowname = function(dat,name){
  write.table(dat,file=name,quote = F,sep="\t",row.names = F,fileEncoding = "utf-8")
}

t.test.mean=function(m1,sd1,n1,m2,sd2,n2){
#  sc = sum()
#  s = 
  t = abs((m1-m2)/sqrt(sd1^2/n1 + sd2^2/n2))
  p = 1-pt(t,n1+n2-2)
  print(t)
  print(p)
}

z.prop = function(x1,x2,n1,n2){
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  r = c(numerator,p.common,z.prop.ris)
  return(r)
}

#1-(1-dnorm(z.prop(44 * 0.26, 40 * 0.58,44,40)))^3
#r= prop.test(x = c(66*0.32, 66*0.48), n = c(66, 66), alternative = "two")#
#r1 = 1 - (1-r$p.value)^3
#r2 = qnorm(r1/2)

#正常值区间
#绘制散点图（序列）
#par(mgp=c(1.6,0.6,0),mar=c(3,3,2,1))



#name=""

main=function(){
  dat = read.table("tmp",header = T,sep = "\t")
  for (i in 1:length(colnames(dat))){
    name = colnames(dat)[i]
    point = na.omit(dat[1:(length(dat[,i])-2),i])
    lower = dat[length(dat[,i])-1,i]
    upper = dat[length(dat[,i]),i]
    plot_mean_diff(name,point,upper,lower)
  }
}

plot_mean_diff = function(name,point,upper,lower){
  Value=seq(min((min(point)-fivenum(point)[2]),lower)-0.5, max((max(point)+fivenum(point)[2]),upper)+0.5,length.out = 20)
  Sample = c(1:length(Value))
  line_type=c(5,5)
  line_width=seq(3,3)
  tiff(filename = paste(name,".tif",collapse = ""),
       width = 2500, height = 2000, units = "px", pointsize = 12,
       compression = "lzw", 
       bg = "white", res = 300)
  plot(Sample,Value,col="blue",pch=5,type="n",font.lab=1,cex.lab=1.5,main=name)
  abline(h=c(upper,lower),lty=line_type,col=colors()[125:126],lwd=line_width)
  points(point)
  dev.off()
}






