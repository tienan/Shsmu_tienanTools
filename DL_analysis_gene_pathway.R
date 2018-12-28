setwd("/home/tienan/R/DL/")
file_list = dir()
file_list_1 = file_list[grepl(file_list,pattern = "csv")]
file_list_2 = file_list_1[grepl(file_list_1,pattern = "CT|CC|base")]

all = read.csv(file_list_2[1])

res = read.csv(file_list_2[2],c(2:3,5))

for (i in 3:length(file_list_2)){
  tmp = read.csv(file_list_2[i])
  res = cbind(res,tmp[,c(3,5)])
}

CT_p10vsp0 = (res[res[,11]<0.05,c(2)])
CT_p10vsp30 = (res[res[,13]<0.05,c(2)])
CT_p30vsp0 = (res[res[,15]<0.05,2])

generate_diff = intersect(intersect(CT_p10vsp0,CT_p10vsp30),CT_p30vsp0)

p10_CT_con_vs_dl = (res[res[,17]<0.05,2])

p30_CT_con_vs_dl = (res[res[,19]<0.05,2])

dl = intersect(p10_CT_con_vs_dl,p30_CT_con_vs_dl)


setdiff(dl,generate_diff)

intersect(dl,generate_diff)


"PGK"%in%intersect(dl,generate_diff)
"PGK"%in%setdiff(dl,generate_diff)
    
target_gene_1 = res[res[,2]%in%setdiff(dl,generate_diff),2]

#DL affacted genes
target_gene_2 = target_gene_1[res[res[,2]%in%target_gene_1 ,16]*
                res[res[,2]%in%target_gene_1,18]>0]

LUAD_DEG_hg38 = read.csv(file = "LUAD_DEG_hg38.csv")
head(LUAD_DEG_hg3)
LUAD_DEG_hg38[grepl(LUAD_DEG_hg38$X,pattern = "PDIA"),]


LUAD_DEG_hg38[LUAD_DEG_hg38[LUAD_DEG_hg38$PValue<0.01,1]%in%target_gene_2,]

write.csv(LUAD_DEG_hg38[LUAD_DEG_hg38[LUAD_DEG_hg38$PValue<0.01,1]%in%target_gene_2,1],"DL_cancer_hg30.csv")

    
head(res)

seq(4,ncol(res),by = 2)
#p:

CC_p10vsp0_diff = res[]




head(tmp)

tmp = read.csv(file_list_2[9])

tmp[tmp$genes=="HMGCR",]

all[all$Name=="HMGCR",]

mean(as.numeric(all[all$Name=="HMGCR",(1:8)+4]))

mean(as.numeric(all[all$Name=="HMGCR",(33:40)+4]))


