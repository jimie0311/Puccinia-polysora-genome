
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
