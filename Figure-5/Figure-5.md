# The heat map of CAZyme and Protease in Figure 5
The CAZyme and Protease annotation were performed by searching GD1913.protein.fa in databases of CAZymes and MEROPS v.12.3. We perfomred in Funannotate v1.8.1.
For more infomration please veiw https://funannotate.readthedocs.io/en/latest/annotate.html
## Generate annotation file from Funannotate
```
funannotate annotate -i fun --cpus 20 --strain "puccinia_polysora" --isolate GD1913
```
The output file will be Puccinia_polysora_GD1913.annotation.txt The annotation results are in two column named as CAZYmes and MEROPS.

For other Puccinia species, download the species_protein.fasta of pca203, pgt21-0 and ptt76.
-download hmmer3.0rc2 from http://hmmer.org
-download all.hmm.ps.len, dbCAN-fam-HMMs.txt and hmmscan-parser.sh from http://csbl.bmb.uga.edu/dbCAN/
```
hmmpress dbCAN-fam-HMMs.txt
hmmscan dbCAN-fam-HMMs.txt species_protein.fasta > CAZyme_species.dbCAN
sh hummscan-parser.sh CAZyme_species.dbCAN > CAZyme_species.annotation
```

Reorgnize all output data to 5A_CAZymes.txt and 5C_proteases.txt
## Plot heat map in R

```
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(pheatmap)
cazymes <- read.table("5A_CAZymes.txt")
proteases <- read.table("5C_proteases.txt")

pdf("Figure5A.pdf", width = 3, height = 7)
pheatmap(cazymes,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white", 
         fontsize_row = 6,
         fontsize_col = 10)
dev.off()

pdf("Figure5C.pdf", width = 3, height = 7.5)
pheatmap(proteases,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white", 
         fontsize_row = 6,
         fontsize_col = 10)
dev.off()
