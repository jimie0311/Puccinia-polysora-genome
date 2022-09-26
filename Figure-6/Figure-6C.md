# The input and plotting code for Figure 6C
We calculated the Copy number variation of species-specified of genes. The basic data are species-specific orthologs in Results_XXX/Comparative_Genomics_Statistics/

## Step 1 Count the orthogroup numbers which show different levels of gene numbers: >40, 31-40, 21-30 and 11-20, for each species.
The input file (genecopy.txt) for plotting Fig.6C as below

```
Num	PpzA	Pt76	Pgt210	Pca203
11-20	18	7	12	13
11-20	16	2	8	10
11-20	19	3	3	8
11-20	11	2	6	3
11-20	9	3	7	1
11-20	8	2	5	2
11-20	9	3	2	2
11-20	9	2	2	0
11-20	6	0	0	0
11-20	7	0	0	0
21-30	4	0	3	1
21-30	7	0	1	1
21-30	2	0	1	2
21-30	1	0	1	2
21-30	10	2	1	0
21-30	2	1	1	0
21-30	4	2	1	0
21-30	2	1	0	0
21-30	0	0	0	0
21-30	0	0	0	0
31-40	4	0	1	1
31-40	1	0	1	0
31-40	2	0	0	0
31-40	1	0	0	0
31-40	1	0	0	0
31-40	3	0	0	0
31-40	1	0	0	0
31-40	2	0	0	0
31-40	0	0	0	0
31-40	0	0	0	0
>40	2	0	0	0
>40	2	0	0	0
>40	3	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	2	0	0	0
>40	2	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
>40	1	0	0	0
```
## Step 2 Plotting in R
```

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(reshape2)
library(ggplot2)

gc <- read.table("genecopy.txt", sep = "\t", header = T)
gc <- melt(gc)
colnames(gc) <- c("CopyNumber", "Species", "GeneNumber")
gc$CopyNumber <- factor(gc$CopyNumber, levels = c("11-20", "21-30", "31-40", ">40"))

ggplot(gc, aes(x = CopyNumber, y =GeneNumber)) +
  geom_boxplot(aes(fill = Species)) + theme_classic() +
  xlab("Gene copy numbers") + ylab("Orthogroup Numbers")+
  theme(legend.position = c(0.9, 0.8), legend.title = element_blank(), legend.text = element_text(face="italic"))
```
