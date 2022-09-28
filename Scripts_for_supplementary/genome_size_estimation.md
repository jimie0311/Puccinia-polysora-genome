# The workflow to estimate the genome size by jellyfish and genomescope

## Step 1 Install and run jellyfish
Pair-wies short sequencing of analyzed species, R1.clean.fq and R2.clean.fq

```
conda install jellyfish
jellyfish count -t 24 -C -m 21 -s 4G -o kmer21.out R1.clean.fq R2.clean.fq
#-s=G+G*c*e*K   G=genome size; c= coverage, e=sequence error rate; K = kmer size
jellyfish dump -c -t kmer21.out -L2 > kmer21.fasta
jellyfish histo kmer21.out -o kmer21.histo
```

## Step 2 Plot kmer distribution by genomescope
```
git clone https://github.com/tbenavi1/genomescope2.0.git 
Rscript genomescope.R -i kmer21.histo -k 21 -o test
```
This step can also run online. upload kmer21.histo to http://qb.cshl.edu/genomescope/genomescope2.0/
