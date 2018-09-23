#TCGA
dat1=read.table("ZZY3.txt",header=T,sep='\t')
a=c()

j=1
for(i in 15111:18116){
  names(dat1)[i]="i"
  dat2=t.test(i~Condition,data=dat1)
  if(dat2$p.value<0.05){
    a[j]=i
    
    j=j+1
  }
}

dat1[,a]

read.table("PM_Heart.txt",header = T,sep = "\t")
