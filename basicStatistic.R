#two-side-t-test
t_test_n = function(d){
  r=c()
  m1 = d[1]
  m2 = d[4]
  sd1 = d[2]
  sd2 = d[5]
  n1 = d[3]
  n2 = d[6]
  m = m1 - m2
  s = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
  s = s*sqrt(1/n1+1/n2)
  t = m/s
  r[1] = abs(round(t,2))
  r[2] = round(1-(pt(r[1],n1+n2-2)),2)*2
  write.table(r,row.names = F,col.names = F, append = FALSE, quote = F,sep="\t", eol = "\t")
  
#  return(r)
} 

#incidence_chisq_test
chisq_test_n = function(d){
  d = as.vector(d)
#  d = c(d[1],d[2]-d[1],d[4],d[5]-d[4])
  
  d=matrix(d,nr=2)
  t_min=min(t_count(d))
  if(t_min>=5){
  t = chisq.test(d,correct = F)
  }else
  t = chisq.test(d,correct = T)
  r=c()
  r[1] = round(as.vector(t$statistic),2)
  r[2] = round(as.vector(t$p.value),2)
  write.table(r,row.names = F,col.names = F, append = FALSE, quote = F,sep="\t", eol = "\t")
  
}

# theory count
t_count = function(x){
  A = matrix(0,nrow(x),ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)) A[i,j] = sum(x[i,])*sum(x[,j])/sum(x)
  }
  A
}

#plot_bar
#plot_bar = function(dat,horiz){
#  a=matrix(c(0,1,1,3,3,2),byrow=T,nrow=3)
#  rownames(a)<-c('S3','S6','S7')
#  colnames(a)<-c('ins','del')
#  barplot(t(a), horiz=horiz)  
#  legend.text = c("ALL", "China")
#  
#}
  
# 