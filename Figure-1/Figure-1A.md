# From HIFI raw data to chromosome assembly
The scripts introduce the pipeline for P. polysora-GD1913 assembly from HIFI raw data to chromosome level

## Step 1 Assemble HIFI raw data via Hicanu v2.1
```
samtools view hifi_raw.bam | awk '{OFS="\t"}; print ">"$1"\n"$10' > GD1913_hifi.fasta
/PATH/TO/canu-2.1.1/build/bin/canu -p GD1913 -d ./ genomeSize=1700000000 -pacbio-hifi /PATH//GD1913_hifi.fasta
```
The output file is GD1913.contigs.fasta

## Step 2 Move contaminated reads and low coverage reads
We followed a workflow described by Jana Sperschneider https://github.com/JanaSperschneider/GenomeAssemblyTools/tree/master/ContaminantScreening
The output file is GD1913.fasta

## Step 3 The mtDNA assembly
The contigs hit to mitochondrion were extracted to a separate FASTA file, In our study, 63 contigs in Table S1 (MER) were used
```
python ~/opt/biosoft/fasta-extract.py ../GD1913.contigs.fasta mit63_ID mit63.fasta
/data/Liangjunmin/opt/biosoft/canu-2.1.1/build/bin/canu -p mit63 -d ./ genomeSize= 70000 -pacbio-hifi mit63.fasta #the output file is mit63.contigs.fasta
samtools faidx mit63.contigs.fasta
makeblastdb -in mit63.contigs.fasta -dbtype nucl -parse_seqids -out mit_out
blastn -query mit63.contigs.fasta -db mit_out -outfmt 7 > mit_blast.out
less mit_blast.out

# BLASTN 2.2.31+
# Query: tig00000001 len=89056 reads=38 class=contig suggestRepeat=no suggestBubble=no suggestCircular=yes trim=16340-84170
# Database: mit_out
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 5 hits found
tig00000001     tig00000001     100.00  89056   0       0       1       89056   1       89056   0.0     1.645e+05
tig00000001     tig00000001     99.94   21229   0       9       67836   89056   1       21224   0.0     39118
tig00000001     tig00000001     99.94   21229   0       9       1       21224   67836   89056   0.0     39118
tig00000001     tig00000001     100.00  39      0       0       30222   30260   30186   30224   8e-13   73.1
tig00000001     tig00000001     100.00  39      0       0       30186   30224   30222   30260   8e-13   73.1
```
1-21224 and 67836-89056 showed that the head and end are matched to a ring. Thus the mtDNA is from 1 to 67835
```
samtools faidx mit63.contigs.fasta tig00000001:1-67835 > GD1913.mtDNA.fasta
```
## Step 4 Preliminary haplotype-phased via Haplomerger2
```
wget https://github.com/mapleforest/HaploMerger2
```
data preparation
-genome fasta file assembled from Canu, FALCON or other de novo softwares after moving contaminants. Here is GD1913.fasta 
-softmask repeats windowmasker
```
windowmasker -checkdup true -mk_counts -in GD1913.clean.fasta -infmt fasta -out GD1913.count -sformat obinary
```
-clean genome by moving illegal charactoer and compress genome 
```
gunzip -c GD1913.fasta.masked.gz | /PATH/TO/HaploMerger2_20180603/bin/fanaPolishing.pl --legalizing --maskShortPortion=1 --noleand:qingN --removeShortSeq=1 > GD1913.cleaned.fasta
```
-update path

All PATH in project_tempate are not the full PATH but ../bin， update as follows, pay attentation to make a copy before updating.
```
cp -r project_template project_template_bck
cd project_template 
grep'\.\.' hm*   #check what you will update
sed -i 's/\.\./\/PATH\/TO\/HaploMerger2_20180603/' hm*
mkdir -p HaploMerger && cd HaploMerger
cp /PATH/TO/HaploMerger2_20180603/project_template/hm* ./
cp /PATH/TO/HaploMerger2_20180603/project_template/*.ctl ./
```
- run HaploMerger
```
sh ./hm.batchA1.initiation_and_all_lastz GD1913.cleaned.fasta  
sh ./hm.batchA2.chainNet_and_netToMaf GD1913.cleaned.fasta
sh ./hm.batchA3.misjoin_processing
```
If the program stops after few seconds, there are something wrong on your PATH
There are two important prameter files (all_lastz.ctl and scoreMatrix.q) which control the behavior of alignment process. 
The file scoreMatrix.q contains a nucleotide substitution score matrix which is used by both LASTZ and chainNet. Because differntt diploid genomes have different heterozygosity rates and GC bias, it is highly recommended to infer a score matrix specific for our own genome.
We need to divide our genome into two parts. The part1.fa included the top few contigs which account for 5%-15% genome size and part2.fa includes all remaining contigs.
If your genome size > 1G, recommend less than 10% genome cotings in partA
```
samtools faidix GD1913.cleaned.fasta
seqkit sort --by-length --reverse GD1913.cleaned.fasta GD1913.cleaned.sorted.fasta
grep -n 'tig0004591' GD1913.cleaned.sorted.fasta # identify contigs accounting for 5% genome size and grep the line number NR of the last contig name
head -n A GD1913.cleaned.sorted.fasta > part1.fa #A = the last line number of the shortest contig in part1.fa
sed -n 'B,Cp'  GD1913.cleaned.sorted.fasta > part2.fa #B = A+1 C=total lines of genome.fasta file
lastz_D_Wrapper.pl --target=part1.fa.gz --query=part2.fa.gz --identity=90  # used lastz_D_Wrapper.pl in HaploMerger2 to calculate scoreMatrix.q
```
Because misjoins can affect the genome alignments, one round of hm.batchA1-3 may not identify all potential misjoins. 
You can use the output diploid assembly file as the input for a second round of hm.batchA1-3. Usually 2-3 rounds of hm.batchA1-3 will be sufficient.

The output files describe the co-linearity violation events:
- hm.scaffolds-information of original scaffolds -
- hm.nodes-information of the alignment blockes
- hm.portions-information of the partions of original scaffolds
- hm.new_scaffolds-the solved relations between two haplotypic scaffolds  #This is the most important output file which includes the all-vs-all alignment
- hm.assemblyz_errors- potential assembly errors in the original assembly
For more details please see https://github.com/mapleforest/HaploMerger2/releases/download/HaploMerger2_20180603/manual.v3.6.pdf

## Step 5 Duplicates exchange
Duplicate means the contig name appeared duplicately within each haplotype or between two haplotypes. These contigs are not sequence repeat but mis-break by HaploMerger2
The exchange criterion of duplicate contigs are as follows. Please see details in Table S3.

- Duplicates between hap1 and hap2
1. if part1 is far larger than part2 and part2 was assigned to a single scaffold with nothing else, I kept all sequences of this contig to the haplotype holding the longer part.
2. If the two parts have similar length in two haplotypes, such as tig00009080 (1.25M, 604KB and 643KB in each haplotype), I split this contig to two haplotypes
3. However, there are many cases. The part1 is not far longer or shorter than part 2, such as tig00060620 (8719bp in hap1 and 5961bp in hap2). I didn't split them and keep all sequences of this contig to the haplotype having the longer part.

- Duplicates within hap1 or hap2
1. Part1 is far longer than part2 and part 2 was assigned to a single scaffold, I kept all sequences to the hap with longer part;
2. Part1 is a little longer than part2 and part 1 was assigned to a single scaffold. such as tig00029892_pilon (6310bp in contig74 vs 15000bp in contig473) , tig00056055 and tig00097736, I also kept all sequences to the scaffold with the longer part. (This dispose can increase scaffold number, Shall we do like this?)
3. part1 is not far longer or shorter than part 2, such as tig00047080 (5584bp in contig80 vs 6961bp in contig128) and tig00060358, I didn't spilit them but kept all sequences to the contig with longer part.

These parts were done in Excel based on hm.new_scaffolds. Finally, generate two contig lists, hapA.ID and hapB.ID
```
python fasta-extract.py GD1913.cleaned.fasta hapA.ID hapA.contig.fasta
python fasta-extract.py GD1913.cleaned.fasta hapB.ID hapB.contig.fasta
```
## Step 6 Assemble contigs to scaffolds
The following scripts assembled contig files to scaffolds by using Hi-C data.
To move experimental artifacts, the Hi-C data went though the mapping workflow from the Arima Genomics pipeline https://github.com/ArimaGenomics/mapping_pipeline/blob/master/01_mapping_arima.sh 
```
mkdir -p raw_dir filter_dir tem_dir pair_dir rep_dir
bwa index hapA.contig.fasta
bwa mem -t 12 hapA.contig.fasta /PATH/HiC/Clean_data/LJM4_R1.fastq.gz | samtools view -@ 10 -Sb > ./raw_dir/LJM4_1.bam
bwa mem -t 12 hapA.contig.fasta /PATH/HiC/Clean_data/LJM4_R2.fastq.gz | samtools view -@ 10 -Sb > ./raw_dir/LJM4_2.bam
samtools view -h ./raw_dir/LJM4_1.bam | perl /PATH/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > ./filter_dir/LJM4_1.bam
samtools view -h ./raw_dir/LJM4_2.bam | perl /PATH/mapping_pipeline/filter_five_end.pl | samtools view -Sb - > ./filter_dir/LJM4_2.bam
perl /PATH/mapping_pipeline/two_read_bam_combiner.pl ./filter_dir/LJM4_1.bam ./filter_dir/LJM4_2.bam samtools 12 | samtools view -bS -t hapA.contig.fasta.fai - | sam    tools sort -@ 12 -o ./tem_dir/LJM4.bam
picard AddOrReplaceReadGroups -INPUT ./tem_dir/LJM4.bam -OUTPUT ./pair_dir/LJM4.bam -ID LJM4 -LB LJM4 -SM GD1913_2 -PL ILLUMINA -PU none
picard MarkDuplicates --INPUT ./pair_dir/LJM4.bam --OUTPUT ./rep_dir/LJM4_LABEL.bam --METRICS_FILE ./rep_dir/LJM4_LABEL.metrics --TMP_DIR ./tem_dir/ --ASSUME_SORTED TRUE --VALIDATION_ST    RINGENCY LENIENT --REMOVE_DUPLICATES TRUE
samtools index ./rep_dir/LJM4_LABEL.bam
perl /PATH/mapping_pipeline/get_stats.pl ./rep_dir/LJM4_LABEL.bam > ./rep_dir/LJM4_LABEL.bam.stats
echo "Finished Mapping Pipeline through Duplicate Removal"
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed
samtools faidx hapA.contig.fasta
python /PATH/SALSA/run_pipeline.py -a hapA.contig.fasta -l hapA.contigs.fasta.fai -b alignment.bed -e GATC -o scaffolds
```  
The output file is hapA/scaffolds/scaffolds_FINAL.fasta
Repeat above steps to obtain hapB/scaffolds/scaffolds_FINAL.fasta
  
## Step 7 Scaffolds rearrangement
Follow Figure-2B, and plot hic-map for each scaffold. Adjust scaffold according to hic map by using emboss
```
sudo apt-get install emboss
#Split each scaffold and wirte all contigs to individual files via seqretsplit
git clone https://github.com/lh3/seqtk.git
cd seqtk && make
seqtk seq -r contig1.fa > contig1_r.fa #Reverse complement FASTA file if hic map shows reverse signal
bedtools getfasta -fi tig00000617.fasta -bed bed.file -o tig00000617_1.fasta #split contigs to two parts based on hm.new_scaffolds if hic map showes mis-join 
```

## Step 8 Rebuild multiple contigs to a single fasta file

Generate a readme file in which the contigs were reordered according to hi-C map. please see an example of readme like below (16 contigs are in scaffold1):
```
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
```
Write the following cycle program in a .sh file and bash it.
```
#!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
        continue
fi
cat ${line}.fasta gap.fasta >> temp.fasta
done<contig.ID

c="Chr1"
echo ${c}
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' temp.fasta > temp1.fasta
union -sequence temp1.fasta -outseq temp2.fasta && rm -f temp.fasta temp1.fasta  # merge multiple contigs to a single fasta file
sed "s/>.*/>${c}/" temp2.fasta > ${c}.fasta && rm -f temp2.fasta
```  
  
Scaffolds rearrangement may take several rounds. 
After each round, each scaffold was rerun the hic-pro and plot new hic-map until the hic-map does not show any mis-break, mis-join or reverse signal.
  
The centromere is detected by viewing hic-map by showing cross-link. Each chromosome only have one centromere. Every chromosome should have two telomeres.
```
perl telomere_identification.pl hapA.scaffolds.fasta > hapA.telomere.out
perl telomere_identification.pl hapB.scaffolds.fasta > hapB.telomere.out
```
Finally, set all scaffolds as the reference genome and map all scaffolds to the genome via HiC-pro, then its hic map will supply connection signal between scaffold.
According to centromere and telomere, large scaffolds were connected to chromosomes.
