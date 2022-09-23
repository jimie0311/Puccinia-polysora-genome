# Figure 6A Synteny Plot
We used MCScan (Python version) to plot the synteny for two haplotypes of P. polysora and haplotype A of three Puccinia species (pca203, pgt21-0 and Ptt76)
## Step 1 Jcvi installtion

```
conda create -y -c bioconda -n jcvi jcvi
conda activate jcvi
pip install git://github.com/tanghaibao/jcvi.git
conda install -y -n jcvi -c bioconda bedtools emboss last seqkit
#install Latex mannually
wget https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz
tar -zxvf install-tl-unx.tar.gz -C ~/opt/biosoft/
cd install-tl-20211207
perl install-tl
![home page of Latex]()

## Step 2 Generate bed and cds files via jcvi
put cds.fa.gz and gene.gff3.gz files of two species/haplotype pairs in four directory

mkdir -
