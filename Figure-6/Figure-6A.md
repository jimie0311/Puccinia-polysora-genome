# Figure 6A Synteny Plot
We used MCScan (Python version) to plot the synteny for two haplotypes of P. polysora and haplotype A of three Puccinia species (pca203, pgt21-0 and Ptt76). For more information of MCscan (python version), pleas view https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)

## Step 1 Jcvi installtion
```
conda create -y -c bioconda -n jcvi jcvi
conda activate jcvi
pip install git://github.com/tanghaibao/jcvi.git
conda install -y -n jcvi -c bioconda bedtools emboss last seqkit
#install Latex mannually
wget https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
tar -zxvf install-tl-unx.tar.gz 
cd install-tl-20211207
perl install-tl
```
## Step 2 Obtain ortholog blocks via jcvi
put cds.fa.gz and gene.gff3.gz files of two species/haplotype pairs in four directory
```
mkdir -p ppzAB ppzA_pt76 pt76_pca203 pca203_pgt210 
```
#take the pairs of ppzA and ppzB as an example
```
#For gff3 files with RNA-seq annotated, each gene has several mrna, choose the first one
python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only ppzA.gff3.gz -o ppzA.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only ppzB.gff3.gz -o ppzB.bed
python -m jcvi.formats.fasta format ppzA.cds.fa.gz ppzA.cds
python -m jcvi.formats.fasta format ppzB.cds.fa.gz ppzB.cds
#compare the gene name in bed file and gff3 file to keep them consistent. Finally obtain ortholog blocks
python -m jcvi.compara.catalog ortolog ppzA ppzB --cscore=.99 --no_strip_names  #the output file is ppzA.ppzB.anchors
```
## Step 3 plot macro-synteny (Figure S3)
```
python -m jcvi.graphics.dotplot ppzA.ppzB.anchors
```

## Step 4 Macosynteny visualiztion of multiple species
- seqids file
each line introduces a species/haplotype with chromosome names separated by ','. 
```
chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18
Pca203_chr1_A,Pca203_chr2_A,Pca203_chr3_A,Pca203_chr4_A,Pca203_chr5_A,Pca203_chr6_A,Pca203_chr7_A,Pca203_chr8_A,Pca203_chr9_A,Pca203_chr10_A,Pca203_chr11_A,Pca203_chr12_A,Pca203_chr13_A,Pca203_chr14_A,Pca203_chr15_A,Pca203_chr16_A,Pca203_chr17_A,Pca203_chr18_A
```
-simple files
```
python -m jcvi.compara.synteny screen --minspan=30 --simple ppzA.ppzB.anchors ppzB.ppzA.anchors.simple
python -m jcvi.compara.synteny screen --minspan=30 --simple ppzA.ppzB.anchors ppzA.pt76.anchors.simple
python -m jcvi.compara.synteny screen --minspan=30 --simple ppzA.ppzB.anchors pt76.pca203.anchors.simple
python -m jcvi.compara.synteny screen --minspan=30 --simple ppzA.ppzB.anchors pca203.pgt210.anchors.simple
```
- layout file, save following parameters to a layout.file
```
#y, xstart, xend, rotation, color, label, va,  bed
.7,     .1,    .8,       0,      , ppzB, top, grape.bed
.5.5,   .1,    .8,       0,      , ppzA, top, peach.bed
.4,     .1,    .6,       0,      , pt76, bottom, cacao.bed
.2.5,   .1,    .5,       0,      , pca203, bottom, cacao.bed
.1,     .1,    .4,       0,      , pgt210, bottom, cacao.bed
#edges
e, 0, 1, ppzB.ppzA.anchors.simple
e, 1, 2, ppzA.pt76.anchors.simple
e, 2, 3, pt76.pca203.anchors.simple
e, 3, 4, pca203.pgt210.anchors.simple
```
