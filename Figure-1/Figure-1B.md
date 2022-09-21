# A workflow for plotting Hi-C map by HiC-Pro and HiCPlotter

## Step 1 HiC-Pro installation
```
git clone https://github.com/nservant/HiC-Pro.git
cd HiC-Pro
source ~/.bashrc.conda #enter in Conda environment
conda create -n hicpro
source activate hicpro
conda install python=2.7
conda install -y samtools bowtiew2 R
conda install -y pysam bx-python numpy scipy
#open R
install.packages(c('ggplot2','RColorBrewer'))
```
Update /PATH/TO/HiC-Pro/config-install.txt

#########################################################################
Paths and Settings  - Start editing here !
#########################################################################
PREFIX = /PATH/TO/miniconda3/Directory/
BOWTIE2_PATH = /PATH/TO/miniconda3_for_pb-assembly/envs/hicpro/bin
SAMTOOLS_PATH = /PATH/TO/miniconda3_for_pb-assembly/envs/hicpro/bin
R_PATH = /PATH/TO/miniconda3_for_pb-assembly/envs/hicpro/bin
PYTHON_PATH = //PATH/TO/miniconda3_for_pb-assembly/envs/hicpro/bin
CLUSTER_SYS = SGE #must be one of TORQUE,SGE,SLURM and LSF

```
make configure
make
```

## Step 2 Preparation for HiC-Pro
```
mkdir -p 01.ref
ln -s /PATH/TO/hapA_all.fasta ./01.ref
/PATH/TO/HiC-Pro_2.11.4/bin/utils/digest_genome.py -r mboi -o hapA_all.fasta.MboI.txt 01.ref/hapA_all.fasta
#-r define the digest enzyme
bowtie2-build --threads 20 /PATH/TO/hapA_all.fasta hapA_all
mkdir -p 02.reads
mkdir -p 02.reads/GD1913
ln -s /PATH/TO/HiC/hic_R1.fq.gz /PATH/TO/HiC/hic_R2.fq.gz ./02.reads/GD1913
samtools faidx /PATH/TO/hapA_all.fasta
awk '{print $1, $2}' hapA_all.fasta.fai > hapA_all.fasta.size
```
## Step 3 Run HiC-Pro
```
cp /PATH/TO/HiC-Pro/config-hicpro.txt ./
```
#update following items
- BOWTIE2_IDX_PATH=/PATH/TO/index_directory/
- REFERENCE_GENOME=hapA_all.fasta
- GENOME_SIZE=/PATH/TO/hapA_all.fasta.size
- GENOME_FRAGMENT=/PATH/TO/hapA_all.fasta.MboI.txt
- LIGATION_SITE=GATCGATC  #hindIII=AAGCTAGCTT MboI=GATCGATC
```
/PATH/TO/HiC-Pro_2.11.1/bin/HiC-Pro -i /PATH/TO/02.reads/ -o hicpro_output -c /PATH/TO/config-hicpro.txt
```
The outputs of HiC-Pro include 2 directoies: bowtie_restuls and hic_results. The directory Matrix under hic_results are the input for HiCPlotter.

## Step 4 HiC map visualization
To present hic map for each contig. a cycle was run as follows.
List each contig name per line in a file named "readme", and run a cycle program as follows.
```
#!/bin/bash
while read line
do
if [[ $line =~ ^# ]];then
        continue
fi
python /PATH/TO/HiCPlotter/HiCPlotter.py -f /PATH/TO/hicpro/hicpro_output/hic_results/matrix/hapA_all/iced/40000/hapA_all_40000_iced.matrix \
-o ${line} -r 40000 -tri 1 -bed /PATH/TO/hic_results/Matrix/hapA_all/raw/40000/hapA_all_40000_abs.bed -n ${line} -Chr ${line}
done<readme
```