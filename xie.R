a1 = c(0:9)
a2 = c(0:9)
a3 = c(0:9)
a4 = c(0:9)
r=data.frame()
m=1
for (i in 1:10){
  for (j in 1:10){
    for (k in 1:10){
      for (l in 1:10){
        a=a1[i]*100+a2[j]*10+5
        b=100+a3[k]*10+a4[l]
        c = a*b
        c=as.character(c)
        c=unlist(strsplit(c,split = ""))
        if(length(c)==5&c[1]=="4"&c[3]=="7"&c[4]=="7"&c[5]=="5"){
          r[m,1]=a
          r[m,2]=b
          m=m+1
        }
      }
    }
  }
}


i=1
j=1
k=1
l=1