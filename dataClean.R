library(randomizr)
data(HairEyeColor)
HairEyeColor <- data.frame(HairEyeColor)

# Transform so each row is a subject
# Columns describe subject's hair color, eye color, and gender
hec <- HairEyeColor[rep(1:nrow(HairEyeColor),
                        times = HairEyeColor$Freq), 1:3]

N <- nrow(hec)

Z <- block_ra(blocks = hec$Hair)
table(Z, hec$Hair)


blocks <- with(hec, paste(Hair, Sex, sep = "_"))
Z <- block_ra(blocks = blocks)
table(blocks, Z)



getwd()

name = read.table("../../Documents/R/thyroidCancerList.txt",header = F)


cmd_1 = "nohup salmon quant -p 16 -i /mnt/Vol01/hg19_38/hg38_index/ -l IU -1 "
cmd_2 = "-2"
cmd_3 = "-o"
cmd_4 = "&"
sink(file = "tmp.txt")

while (i<nrow(name)) {
  


tmp = unlist(strsplit(as.character(name[i,1]), "_AH"))
tmp[1]

print(paste(cmd_1,as.character(name[i,1]),cmd_2,as.character(name[i+1,1]),cmd_3, tmp[1], cmd_4,sep = " "))
i=i+2


}
sink()
