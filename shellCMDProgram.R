# salmon exp calculating
library(limma)
dat = read.table("../R/file_list",header = F)
cmd = c()
#salmon quant -i ../../hg19_38/GRCh38_latest_rna_salmon_index -l A -1 AB-1_1.fq.gz -2 AB-1_2.fq.gz -p 16 --validateMappings -o AB-1
cmd=c()
i=1
while (i < nrow(dat)){
  tmp = strsplit2(as.character(dat[i,]),split = "/") 
  tmp_1 = strsplit2(tmp[2],split = "_1")
  tmp_2 = paste("nohup salmon quant -i ../../hg19_38/GRCh38_latest_rna_salmon_index -l A -1 ", dat[i,]," -2 ",
              dat[i+1,]," -p 16 --validateMappings -o  ", tmp_1[1],"   &",sep = " ") 
  cmd = rbind(cmd,tmp_2)
  i = i+2
  
}
cmd
write.table(file = "DLSeq.sh",row.names = F,col.names = F,x = cmd,quote = F,fileEncoding = "UTF-8")
getwd()
