# Scripts for plotting PCA

We used plink to perform PCA analysis. Before analyzing, gatk v4.1.9 was applied to obtain vcf file of 77 P. polysora isolates. The rD was calculated by R package _poppr_

## Obtain combined.vcf of P. polysora population via gatk v4.1.9

```
samtools faidx hapA_all.fasta
java -jar /PATH/TO/picard-tools-2.18.27/picard.jar CreateSequenceDictionary R=hapA_all.fasta O=genome.dict
```
#list all isolate name in isolate.list file, do a cycle as follows
```
!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
        continue
fi
/PATH/gatk-4.1.9.0/gatk AddOrReplaceReadGroups -I ../bwa/${line}.rd.bam -O ${line}.group.bam -LB ${line} -PL Illumina -PU run -SM ${line}
samtools index ${line}.group.bam ${line}.group.bam.bai
/PATH/gatk-4.1.9.0/gatk --java-options -Xmx4G HaplotypeCaller -R hapA_all.fasta -I ${line}.group.bam -O ${line}.gatk.gvcf --emit-ref-confidence GVCF
done<isolate.list

#combine gvcf files
input_gvcf=(/mnt/data/liangjunmin/genome_resequencing/gatk_with_pgt/*.gvcf)
GVCF_IN=()
for s in "${input_gvcf[@]}"
do
GVCF_IN+=(-V $s)
done
gatk CombineGVCFs -O combine.gvcf "${GVCF_IN[@]}" -R hapA_all.fasta 

#convert combined.gvcf to combine.vcf
/PATH/gatk-4.1.9.0/gatk GenotypeGVCFs -R hapA_1aa.fasta -V combine.gvcf -O combine.vcf

#extract and filter snp.vcf
gatk IndexFeatureFile -I combine.vcf
gatk SelectVariants -V combine.vcf -select-type SNP -O snp.vcf
gatk VariantFiltration -V snp.vcf --filter-expression 'QD < 2.0 || QUAL < 30.0 || MQ < 40.0 ||FS > 60.0' --filter-name "filter" -O snp.labeled.vcf
awk '/^#/||$7=="PASS"' snp.labeled.vcf > snp.filtered.vcf
```

## Perform PCA analysis

move contaminated isolates if have
```
bcftools view -S clean_name.txt snp.filtered.vcf -Ov > cleaned.snp.vcf
```
Run PCA analysis via plink
```
vcftools --vcf cleaned.snp.vcf --plink --out snp.plink.vcf
/PATH/TO/plink --vcf snp.plink.vcf --recode --out SNP_OUT --allow-extra-chr --vcf-half-call m --noweb
/PATH/TO/plink --allow-extra-chr --file SNP_OUT --noweb --make-bed --out SNP_OUT
/PATH/TO/plink --allow-extra-chr --threads 30 -bfile SNP_OUT --pca 50 --out SNP_OUT
awk '{print $2"\t"$3"\t"$4}' SNP_OUT.eigenvec > SNP_OUT.pca.txt
```

## Plot PCA in R
```

#a <- read.table ("/PATH/pca/SNP_OUT.eigenvec",header=TRUE)
#plot(a$pc1,a$pc2,pch=c(rep(1,30),rep(2, 31),rep(3,19)),col=c(rep("blue",30),rep("red",31),rep("black",19)),main="pca", xlab="pc1",ylab="pc2")
#legend("bottomleft",c("south","north","central"), pch=c(rep(1),rep(2),rep(3)), col=c("blue","red","black"))
install.packages("FactoMineR")
install.packages("factoextra")
library("FactoMineR")
library("factoextra")
library("ggplot2")
setwd("E:/RStudio/pca/")
data <- read.table("SNP_OUT.pca.txt", header = T, sep="\t")
#ggplot(data, aes(x=PC1,y=PC2, color=group))+geom_point(shape=19,size=4)+scale_shape_manual(values=c(1:8))"
pdf(file="PCA.pdf")
ggplot(data=data, aes(x = PC1, y = PC2, fill = group)) + xlab("PC1 (3.10%)") + ylab("PC2 (3.0%)") + geom_point(size = 5,color = "black", shape=21) + theme(axis.title = element_text(size = 15),axis.text = element_text(size=15, color="black"), legend.title=element_blank(), legend.text=element_text(size=15,color="black"),panel.border = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(size = 0.8)) + ylim(-0.3,0.4)+xlim(-0.25,0.5)
dev.off()
#pca file column 1-3 = isolate name, pc1, pc2
```




