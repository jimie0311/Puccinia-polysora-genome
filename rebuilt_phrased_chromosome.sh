#step 1 read hapA.fasta from SALSA by seqretsplit and wrote each scaffold to individual files
#Step 2 seqtk seq -r contig1.fa > contig1_r.fa #Based on hi-C map if inversion is necessary
#step 3 bedtools getfasta -fi tig00000617.fasta -bed bed.file -o tig00000617_1.fasta #based on hi-C map if breaking is necessary
#step 4 generate a readme file in which the contigs were reordered according to hi-C map. please see an example of readme like below (16 contigs are in scaffold1):
   tig00000617_2
   tig00002674
   tig00002981_r
   tig00000617_1_r
   tig00002822
   tig00002255
   tig00003203
   tig00002121_r
   tig00002643
   tig00002784_r
   tig00002744
   tig00001397_r
   tig00001374
   tig00002127
   tig00002362_r
   tig00003315
#step 5 create a gap.fasta with 100 N
#Step 6 rebuild the contigs together to a new chromosome
#putting together a chromsome from mutiple contig
#!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
         continue
fi
cat ${line}.fasta gap.fasta >> temp.fasta
done<readme
c="Chr1"
echo ${c}
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' temp.fasta > temp1.fasta
union -sequence temp1.fasta -outseq temp2.fasta && rm -f temp.fasta temp1.fasta
sed "s/>.*/>${c}/" temp2.fasta > ${c}.fasta && rm -f temp2.fasta
