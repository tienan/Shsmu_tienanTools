# TCGA 
find -name "*isoforms.nor*" | xargs paste | awk '{{printf "%s",$1}{printf "\t"}for(i=2;i<=NF;i=i+2){printf "%s\t",$i}{printf "\n"}}' > ../isoform.exp 
find -name "*isoforms.nor*" | awk '{{printf "%s",substr($1,3,length($1))}{printf "\t"}}' > ../head.txt
cd ..
cat head.txt | awk '{{printf "ID"}{printf "\t"}{print $0}}' > head_1.txt
sed -i '1d' isoform.exp
cat head_1.txt isoform.exp > isoform_1.exp

len=`awk '{print NF}' head_1.txt`
i=1
while(($i <= $len))
  do cut -f"$i"  isoform_1.exp | paste -s
  let i+=1
done > isoform.exp.out
cd ..
cd ..
cat FILE_SAMPLE_MAP.txt | grep isoforms.nor | awk '{print $1"\t"substr($2,1,12)"\t"substr($2,14,2)}' > connector.txt

echo -e  "ID\tPatient\tCondition" > connector.txt.head
cat connector.txt.head connector.txt > connector_1.txt



#for i in {1..554}; do cut -f"$i" isoform_1.exp | paste -s; done | less

#var=`echo $result|awk '{print substr($result,16,3)}'`; 
# TCGA-18-3407
# 关联ID

#find -name "*isoforms.nor*" | xargs awk '{for(i=1;i<=NR;i=i+1){printf "%s\t%s",$1,$2}{printf "\n"}}' | less


####基因
