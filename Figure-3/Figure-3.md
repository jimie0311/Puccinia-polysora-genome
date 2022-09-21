# Figure-3A Estimation of Structure Variation
We used Assemblytics to analyze the structure variation between two haplotypes.
Assemblytics is a web app for detecting and analyzing variants from a de novo genome assembly aligned to a reference genome.
Here is the website http://assemblytics.com
```
git clone https://github.com/MariaNattestad/Assemblytics.git
chmod a+x scripts/Assemblytics*
activate conda  #sourece ~/.bashrc.conda
conda create -n mumer4 && conda activate mumer4
nucmer --maxmatch -l 100 -c 500 REFERENCE.fa ASSEMBLY.fa --prefix OUT
gzip OUT.delta
```
upload the .delta.gz file to Assemblytics (http://assemblytics.com)

![Assemblytics Page](/Puccinia-polysora-genome/Figure-3/assemblytics.png)

For more information of mummer, pleas view http://mummer.sourceforge.net/manual/#usecases

# Figure-3B
mummer --maxmatch -t 100 -b 200 -c 200 -p hapAB hapA.fasta hapB.fasta
delta-filter -i 95 -l 1000 -1 hapAB.delta > hapAB.filter
show-coords -c -d -l -I 95 -L 10000 -r hapAB.filter > hapAB.coord
mumerplot --png -p hapAB.1000 hapAB.filter

# Figure-3C and 3D
Compare chromosomes between two haplotypes, a cycle as follows was run.
List all chromosome name in a readme file
```
  !/bin/bash
  while read line
  do
  if [[ $line =~ ^# ]];then
          continue
  fi
  nucmer --maxmatch -t 100 -b 200 -c 200 -p chr01a.${line} chr01a.fasta ${line}.fasta
  delta-filter -i 98 -l 2000 -1 chr01a.${line}.delta > chr01a.${line}.filter
  show-coords -c -d -l -I 95 -L 10000 -r ${line}.delta >${line}.coord
  mummerplot --png -p ${line}.1000 ${line}.filter
  done < readme
```
