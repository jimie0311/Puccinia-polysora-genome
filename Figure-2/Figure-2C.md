# Ageing LTR insertions
We followed a pipeline https://github.com/SIWLab/Lab_Info/wiki/Ageing-LTR-insertions with few modification

## step 1 LTR_retriever installation
```
conda install LTR_retriever
```

## step 2 Identification of intact LTR-RT via LTR_harvest
```
gt suffixerator -db hapA_all.fasta -indexname hapA_all -tis -suf -lcp -des -ssp -sds -dna
gt harvest -index hapA_all -v -gff3 hapA_all.gff3 -out hapA_all.ltr.out > hapA_all.stdout
```

## step 3 Extract bed files for LTR_RTs
LTRharvest doesn't use the scaffold name you give it but outputs as seqX where Xis the 0-base number of your scaffold. 
Please remove the first few lines started with "#" and add "###" at the first line. see an example as below.

```
ls hapA_all.gff3

##gff-version 3
##sequence-region   seq0 1 67484389
##sequence-region   seq1 1 55950575
##sequence-region   seq2 1 55030165
##sequence-region   seq3 1 44624142
##sequence-region   seq4 1 50342496
##sequence-region   seq5 1 42457113
##sequence-region   seq6 1 55595191
##sequence-region   seq7 1 55557465
##sequence-region   seq8 1 54465677
#chr01
#chr02
#chr03
#chr04
#chr05
#chr06
#chr07
#chr08
#chr09
seq0	LTRharvest	repeat_region	4380	13346	.	?	.	ID=repeat_region1
seq0	LTRharvest	target_site_duplication	4380	4383	.	?	.	Parent=repeat_region1
seq0	LTRharvest	LTR_retrotransposon	4384	13342	.	?	.	ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=86.67;seq_number=0
seq0	LTRharvest	long_terminal_repeat	4384	4514	.	?	.	Parent=LTR_retrotransposon1
seq0	LTRharvest	long_terminal_repeat	13208	13342	.	?	.	Parent=LTR_retrotransposon1
seq0	LTRharvest	target_site_duplication	13343	13346	.	?	.	Parent=repeat_region1
###
seq0	LTRharvest	repeat_region	28572	31425	.	?	.	ID=repeat_region2
seq0	LTRharvest	target_site_duplication	28572	28575	.	?	.	Parent=repeat_region2
seq0	LTRharvest	LTR_retrotransposon	28576	31421	.	?	.	ID=LTR_retrotransposon2;Parent=repeat_region2;ltr_similarity=93.89;seq_number=0
seq0	LTRharvest	long_terminal_repeat	28576	28886	.	?	.	Parent=LTR_retrotransposon2
seq0	LTRharvest	long_terminal_repeat	31113	31421	.	?	.	Parent=LTR_retrotransposon2
seq0	LTRharvest	target_site_duplication	31422	31425	.

ls hapA_all.gff3.rename
###
chr01	LTRharvest	repeat_region	4380	13346	.	?	.	ID=repeat_region1
chr01	LTRharvest	target_site_duplication	4380	4383	.	?	.	Parent=repeat_region1
chr01	LTRharvest	LTR_retrotransposon	4384	13342	.	?	.	ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=86.67;seq_number=0
chr01	LTRharvest	long_terminal_repeat	4384	4514	.	?	.	Parent=LTR_retrotransposon1
chr01	LTRharvest	long_terminal_repeat	13208	13342	.	?	.	Parent=LTR_retrotransposon1
chr01	LTRharvest	target_site_duplication	13343	13346	.	?	.	Parent=repeat_region1
###
chr01	LTRharvest	repeat_region	28572	31425	.	?	.	ID=repeat_region2
chr01	LTRharvest	target_site_duplication	28572	28575	.	?	.	Parent=repeat_region2
chr01	LTRharvest	LTR_retrotransposon	28576	31421	.	?	.	ID=LTR_retrotransposon2;Parent=repeat_region2;ltr_similarity=93.89;seq_number=0
chr01	LTRharvest	long_terminal_repeat	28576	28886	.	?	.	Parent=LTR_retrotransposon2
chr01	LTRharvest	long_terminal_repeat	31113	31421	.	?	.	Parent=LTR_retrotransposon2
chr01	LTRharvest	target_site_duplication	31422	31425	.
```
## Extract LTR bed files via bedtools
```
mkdir -p ltr && cd ltr
cp ../hapA_all.gff3.rename ./
ulimit -n 100000  #Lift the limit on the number of files, otherwise an error will be reported
awk  '/###/{filename=NR".txt"}; /long_terminal_repeat/{gsub(/Parent=/,"");print $1,$4,$5,$9  >filename}' OFS="\t"  hapA_all.gff3.rename  #generate ltr bed files of each LTR-RT 
ls | grep txt  | sort -n > file.list   #generate list of bed files for each LTR-RT
cp file.list ../
cd ../
```

## Step 4 Extract sequences of LTR-RTs and perform allignment via muscle
wirte a cycle to a shell file named extract.sh and bash extract.sh
```
#!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
continue
fi
bedtools getfasta -fi ../hapA_all.fasta -bed ${line} > ${line}.fa
muscle -in ${line}.fa -out ${line}.afa
done<file.list
```
The outputs are hundreds or thousands of afa files and with 100000 afa files in each directory.

## Step 5 Insertion time estimation via Ape in R
```
library(ape)
library("xlsx")
setwd("E:/rstudio/ltr")

list <- list.files()
fas.F1 = read.FASTA(list[1])
mat1 = dist.dna(fas.F1,as.matrix = T, model = "K80")
merge.data = mat1[1,2]
mutate_rate <- 2e-8  # the 2.0×10-8 mutations per nucleotide site of a closely related fungus, Schizophyllum commune （Baranova et al. 2015）
time1 = merge.data/(2*mutate_rate)
v1.merge = c(list[1],merge.data, time1)

dir = paste("./",list,sep="")
n = length(dir)
for (i in 2:n)
  {
  fas.F.new = read.FASTA(list[i])
  mat.new = dist.dna(fas.F.new, as.matrix = T, model = "K80")
  time = mat.new[1,2]/(2*mutate_rate)
  v.new = c(list[i], mat.new[1,2], time)
  v1.merge = rbind(v1.merge, v.new)
}

write.table(v1.merge, file = "ltr_time.out")
```
## Step 6 LTR-RT annotation via REXdb
``
conda install tesorter  #https://github.com/zhangrengang/TEsorter#outputs
cat *.txt.fa > intact_ltr.fa
tesorter intact_ltr.fa  #the output is intact_ltr.fa.rexdb.cls.tsv
```
match insert time (mya) from ltr_time.out and LTR-type from intact_ltr.fa.rexdb.cls.tsv via LTR name and generate ltr.txt
```

## Step 7 Ploting of Figure 2C

```
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)

ltr <- read.table(file = "ltr.txt", header = T, sep = "\t") %>% mutate(Mya=Year/1e6) %>% mutate(LTRtype = LTR_superfamily)
ltr$LTRtype[which(ltr$LTRtype == 'Bel-Pao' )] <- 'Others'
ltr$LTRtype[which(ltr$LTRtype == 'MuDR_Mutator' )] <- 'Others'
ltr$LTRtype[which(ltr$LTRtype == 'PIF_Harbinger' )] <- 'Others'
ltr$LTRtype[which(ltr$LTRtype == 'Retrovirus' )] <- 'Others'
ltr$LTRtype[which(ltr$LTRtype == 'unknown' )] <- 'Others'
ltr_keep <- filter(ltr, LTRtype != "")

scaleFUN <- function(x) sprintf("%.2f", x) #set digits for x axis label
ggplot(ltr) + 
geom_line(ltr[ltr$LTRtype != "", ], mapping = aes(x = Mya, y=..count../nrow(ltr_keep), color = LTRtype), stat="bin", binwidth = 0.5) +
geom_line(ltr, mapping = aes(x = Mya, y=..count../nrow(ltr)),color = "black", stat="bin", binwidth = 0.5) +
labs(x = 'Insertion time (Mya)', y = 'Frequency') +
scale_color_manual(values = c( "#4ea5b9", "#f08f54", "#61733b")) +
theme_classic()+
scale_y_continuous(expand = c(0,0))+
scale_x_continuous(expand = c(0,0), limits=c(0,7), labels=scaleFUN) +
theme(legend.position = c(0.9, 0.8), legend.title = element_blank())

pdf(file = "LTR1.pdf", width = 5, height = 5)
last_plot()  
dev.off()
```
