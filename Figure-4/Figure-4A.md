# The workflow to perform clustering analysis of gene expression of secretome 
## Step 1 Generate count matrix of expression data
- index genome
```
hisat2-build -f -p 20 ./genome/GD1913_hapA.fasta genome/GD1913_hapA_index
hisat2-build -f -p 20 ./genome/GD1913_hapB.fasta genome/GD1913_hapB_index
```
- genome mapping

write the names of RNA-seq data for each time-point in sample_name.txt, run a cycle as follows

```
while read line
do
#hapA
hisat2 -p 20 -x ./genome/GD1913_hapA_index -1 data/${line}_R1.fq.gz -2 data/${line}_R2.fq.gz -S bam_files_hapA/${line}_GD1913_hapA.sam --rg-id ${line} --rg SM:${line} 2> bam_files_hapA/${line}_GD1913_hapA_mappingrate.txt
samtools sort -@ 15 ./bam_files_hapA/${line}_GD1913_hapA.sam -o ./bam_files_hapA/${line}_GD1913_hapA_sorted.bam --reference genome/GD1913_hapA.fasta
samtools index ./bam_files_hapA/${line}_GD1913_hapA_sorted.bam
#hapB
hisat2 -p 20 -x ./genome/GD1913_hapB_index -1 data/${line}_R1.fq.gz -2 data/${line}_R2.fq.gz -S bam_files_hapB/${line}_GD1913_hapB.sam --rg-id ${line} --rg SM:${line} 2> bam_files_hapB/${line}_GD1913_hapB_mappingrate.txt
samtools sort -@ 15 ./bam_files_hapB/${line}_GD1913_hapB.sam -o ./bam_files_hapB/${line}_GD1913_hapB_sorted.bam --reference genome/GD1913_hapB.fasta
samtools index bam_files_hapB/${line}_GD1913_hapB_sorted.bam
done < sample_names.txt
```
- Generate stringtie expression files from each time-point sample
```
inputA=(./bam_files_hapA/*_GD1913_hapA.bam)
for ((i=0;i<=${#inputA[@]};i++))
do
stringtie ${inputA[i]} -G ./genome/Puccinia_polysora_GD1913_hapA.gff3 -e -o ./expression/$(basename "${input[i]%_GD1913_hapA.bam}_hapA_expression") -p 10
done
inputB=(./bam_files_hapA/*_GD1913_hapB.bam)
for ((i=0;i<=${#inputB[@]};i++))
do
stringtie ${inputB[i]} -G ./genome/Puccinia_polysora_GD1913_hapB.gff3 -e -o ./expression/$(basename "${input[i]%_GD1913_hapB.bam}_hapB_expression") -p 10
done
```
- Generate count matrix
```
ls -l | grep "expression" | awk '{print $9"\t"$9}' | awk 'BEGIN{OFS="\t"}{gsub (/_hapA_expression/,"",$1); print $0}' > GFF_listA.txt
prepDE.py -i GFF_listA.txt -g ./expression/GD1913_hapA_count_matrix_genes -t ./expression/GD1913_hapA_count_matrix_transcripts -l 150
ls -l | grep "expression" | awk '{print $9"\t"$9}' | awk 'BEGIN{OFS="\t"}{gsub (/_hapB_expression/,"",$1); print $0}' > GFF_listB.txt
prepDE.py -i GFF_listB.txt -g ./expression/GD1913_hapB_count_matrix_genes -t ./expression/GD1913_hapB_count_matrix_transcripts -l 150
```

## Step 2 Plot clustering in R
```
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(DESeq2)
library(ggplot2)
library(cowplot)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# DESeq2 standardization

count_trans_hapA <- data.frame(read.delim2("Data/GD1913_hapA_count_matrix_transcripts", header = T, sep = ",", row.names = 1))
count_trans_hapB <- data.frame(read.delim2("Data/GD1913_hapB_count_matrix_transcripts", header = T, sep = ",", row.names = 1))
count_trans_hapA <- count_trans_hapA[,c(19:21, 1:18)] #move GS lane to first column
count_trans_hapB <- count_trans_hapB[,c(19:21, 1:18)] #move GS lane to first column
colnames(count_trans_hapA) <- c("GS1","GS2","GS3","dpi01-1","dpi01-2","dpi01-3","dpi02-1","dpi02-2","dpi02-3","dpi04-1","dpi04-2","dpi04-3","dpi07-1","dpi07-2","dpi07-3","dpi10-1","dpi10-2","dpi10-3","dpi14-1","dpi14-2","dpi14-3")
colnames(count_trans_hapB) <- c("GS1","GS2","GS3","dpi01-1","dpi01-2","dpi01-3","dpi02-1","dpi02-2","dpi02-3","dpi04-1","dpi04-2","dpi04-3","dpi07-1","dpi07-2","dpi07-3","dpi10-1","dpi10-2","dpi10-3","dpi14-1","dpi14-2","dpi14-3")
#load names of secreted proteins
secreted_hapA <- read.table("Data/secreted_hapA.txt")
secreted_hapB <- read.table("Data/secreted_hapB.txt")

design <- data.frame(read.delim2("Data/design.txt", header = T))

#create dds(DeseqDataSet)
ddsA <- DESeqDataSetFromMatrix(countData = count_trans_hapA, colData = design, design = ~Condition)
ddsB <- DESeqDataSetFromMatrix(countData = count_trans_hapB, colData = design, design = ~Condition)

#Pre-filtering
keepA <- rowSums(counts(ddsA)) >=10
ddsA <- ddsA[keepA,]
keepB <- rowSums(counts(ddsB)) >=10
ddsB <- ddsB[keepB,]

#de analysis
ddsA <- DESeq(ddsA)
ddsB <- DESeq(ddsB)

#rld
##choose secreted proteins for following analyze after rlog transform
rldA <- rlog(ddsA, blind=FALSE)
selectA <- vector()
for (i in 1:nrow(secreted_hapA)){ selectA <- c(selectA, which(rldA@rowRanges@partitioning@NAMES == secreted_hapA[i,1]))}

rldB <- rlog(ddsB, blind=FALSE)
selectB <- vector()
for (i in 1:nrow(secreted_hapB)){ selectB <- c(selectB, which(rldB@rowRanges@partitioning@NAMES == secreted_hapB[i,1]))}

datA <- assay(rldA)[selectA,]
datA <- datA - rowMeans(datA)
datA <- as.matrix(datA)

datB <- assay(rldB)[selectB,]
datB <- datB - rowMeans(datB)
datB <- as.matrix(datB)

#Within the Sum of Squares
wssA <- (nrow(datA)-1)*sum(apply(datA,2,var))
wssB <- (nrow(datB)-1)*sum(apply(datB,2,var))

caA <- c(rep(0,30)) #ca means change amount of Variance
caB <- c(rep(0,30))

for (i in 2:30){
  set.seed(1234)
  wssA[i] <- sum(kmeans(datA, centers=i)$withinss)}

for (i in 2:29){
  caA[i] <- abs(wssA[i+1] - wssA[i])
}
wccPlotA <- ggplot(data.frame(wssA), aes(x = c(1:30), y = wssA)) + geom_line() + geom_point() + theme_bw() + xlab("Number of Clusters") + ylab("Within groups sum of squares of A")
caPlotA <- ggplot(data.frame(caA), aes(x = c(1:30), y = caA)) + geom_line() + theme_bw() + xlab("Number of Clusters") + ylab("Change of wss A")


for (i in 2:30){
  set.seed(12345)
  wssB[i] <- sum(kmeans(datB, centers=i)$withinss)}

for (i in 2:29){
  caB[i] <- abs(wssB[i+1] - wssB[i])
}
wccPlotB <- ggplot(data.frame(wssB), aes(x = c(1:30), y = wssB)) + geom_line() + geom_point() + theme_bw() + xlab("Number of Clusters") + ylab("Within groups sum of squares of B")
caPlotB <- ggplot(data.frame(caB), aes(x = c(1:30), y = caB)) + geom_line() + theme_bw() + xlab("Number of Clusters") + ylab("Change of wss B")


plot_grid(wccPlotA, caPlotA, wccPlotB, caPlotB)
ggsave("Plot/elbow_bend_secreted_protein.pdf", last_plot(), width = 10, height = 6)

library(cluster)
library(fpc)
#k-mean clustering
clusA <- kmeans(datA, centers=7)
clusB <- kmeans(datB, centers=7)

clusplot(datA, clusA$cluster, color=TRUE, shade=TRUE, labels=4, lines=0)
clusplot(datB, clusB$cluster, color=TRUE, shade=TRUE, labels=4, lines=0)

plotcluster(datA, clusA$cluster)
plotcluster(datB, clusB$cluster)

#obtain expression file with ordered cluster
kA <- kmeans(datA,7)
m.kmeansA<- cbind(datA, kA$cluster) # combine the cluster with the matrix
oA<- order(m.kmeansA[,22]) # order the last column
m.kmeansA<- m.kmeansA[oA,] # order the matrix according to the order of the last column
kB <- kmeans(datB, 7)
m.kmeansB<- cbind(datB, kB$cluster) # combine the cluster with the matrix
oB<- order(m.kmeansB[,22]) # order the last column
m.kmeansB<- m.kmeansB[oB,] # order the matrix according to the order of the last column

#calculate mean
mm.kmeansA <- t(m.kmeansA)[0]
mm.kmeansB <- t(m.kmeansB)[0]
for (i in seq(1,21,3)){
  mm.kmeansA  <- cbind(mm.kmeansA, apply(m.kmeansA[,c(i,i+1,i+2)], MARGIN = 1, mean))
  mm.kmeansB  <- cbind(mm.kmeansB, apply(m.kmeansB[,c(i,i+1,i+2)], MARGIN = 1, mean))
}
colnames(mm.kmeansA) <- unique(design$Condition)
colnames(mm.kmeansB) <- unique(design$Condition) 

#Plot circos map
#hap A
split = factor(as.vector(m.kmeansA[,22]))
circos.clear()
circos.par(start.degree = 80, gap.after = c(3, 3, 3, 3, 3, 3, 20))
colMain = colorRamp2(as.vector(seq(from = -10, to = 6.0, length.out = 100)), 
                     rev(colorRampPalette(brewer.pal(11, "Spectral"))(100)))
colK <- c("1" = "black", 
          "2" = "black", 
          "3" = "black",
          "4" = "black",
          "5" = "black",
          "6" = "black",
          "7" = "black"
)
circos.heatmap(as.vector(split), col = colK, track.height = 0.01, split = split)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(2, "mm"), 
              CELL_META$sector.index,
              facing = "bending.inside", cex = 0.8,
              adj = c(0.5, 0), niceFacing = TRUE)
}, bg.border = NA)

circos.heatmap(mm.kmeansA[,1], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,2], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,3], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,4], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,5], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,6], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansA[,7], col = colMain, split = split, track.height = 0.07, cluster = F)
range(mm.kmeansA)
circos.clear()
lgd = Legend(col_fun = colMain)
grid.draw(lgd)

#hap B
split = factor(as.vector(m.kmeansB[,22]))
circos.clear()
circos.par(start.degree = 80, gap.after = c(3, 3, 3, 3, 3, 3, 20))
colMain = colorRamp2(as.vector(seq(from = -10, to = 6.0, length.out = 100)), 
                     rev(colorRampPalette(brewer.pal(11, "Spectral"))(100)))
colK <- c("1" = "black", 
          "2" = "black", 
          "3" = "black",
          "4" = "black",
          "5" = "black",
          "6" = "black",
          "7" = "black"
)
circos.heatmap(as.vector(split), col = colK, track.height = 0.01, split = split)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(2, "mm"), 
              CELL_META$sector.index,
              facing = "bending.inside", cex = 0.8,
              adj = c(0.5, 0), niceFacing = TRUE)
}, bg.border = NA)

circos.heatmap(mm.kmeansB[,1], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,2], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,3], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,4], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,5], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,6], col = colMain, split = split, track.height = 0.07, cluster = F)
circos.heatmap(mm.kmeansB[,7], col = colMain, split = split, track.height = 0.07, cluster = F)

circos.clear()
lgd = Legend(col_fun = colMain)
grid.draw(lgd)

#save txt files
secreted_K_A <- data.frame(m.kmeansA[,22])
secreted_K_B <- data.frame(m.kmeansB[,22])
colnames(secreted_K_A)[1] <- "K"
colnames(secreted_K_B)[1] <- "K"
write.table(secreted_K_A, file = "Data/secreted_K_A.txt", quote = F)
write.table(secreted_K_B, file = "Data/secreted_K_B.txt", quote = F)
```
