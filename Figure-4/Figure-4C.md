# RScripts used to plot avrRppC position
All related input files are uploaded in Figure-4
```
setwd("E:/RStudio/circlize_update/")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(karyoploteR)

custom.genome <- toGRanges(data.frame(chr=c("chr14A", "chr14B"), start=c(1, 1), end=c(34278376, 34814183)))
plot.params$data1height <- 80
pp$leftmargin <- 0.5
custom.cytobands <- toGRanges("chr14.txt")
kp <- plotKaryotype(plot.type = 1, genome = custom.genome, cytobands = custom.cytobands)
#divide data.lane1 to three layers
kpDataBackground(kp, r0 = 0,r1 = 0.30)
kpDataBackground(kp, r0 = 0.35,r1 = 0.65,col="#FFFACD")
kpDataBackground(kp, r0=0.70,r1=1,col="#B0E2FF")
# label scales
kpAddBaseNumbers(kp, tick.dist=1e7, tick.len=10, add.units=TRUE, digits=2, minor.ticks=TRUE, 
                 minor.tick.dist=10000, minor.tick.len=2,  cex=0.8, tick.col=NULL, minor.tick.col=NULL, clipping=TRUE)
#kpDataBackground(kp, data.panel = 2, col="#FFFACD")
# lable line logTPM data
kpAddLabels(kp, labels="log2TPM", r0=0, r1=0.3, cex=0.5)
kpAxis(kp,r0=0,r1=0.3, ymin=0, ymax=15,cex=0.7)
LineData_hapA <- read.table("chr14A_100kb.TPM") 
LineData_hapB <- read.table("chr14B_100kb.TPM")
LineData_hapA[,5:11] <- log2(LineData_hapA[,5:11]+1)
LineData_hapB[,5:11] <- log2(LineData_hapB[,5:11]+1)
LineData_hapA <-toGRanges(data.frame(chr=LineData_hapA$V1, start=LineData_hapA$V2, end=LineData_hapA$V3, GS=LineData_hapA$V5, dpi_1=LineData_hapA$V6, dpi_4=LineData_hapA$V8, dpi_7=LineData_hapA$V9, dpi_13=LineData_hapA$V11))
LineData_hapB <-toGRanges(data.frame(chr=LineData_hapB$V1, start=LineData_hapB$V2, end=LineData_hapB$V3, GS=LineData_hapB$V5, dpi_1=LineData_hapB$V6, dpi_4=LineData_hapB$V8, dpi_7=LineData_hapB$V9, dpi_13=LineData_hapB$V11))
kpLines(kp, data=LineData_hapA, ymin=0, ymax=15, y=LineData_hapA$dpi_1, data.panel=1, col="#0000FF", r0=0, r1=0.3, lwd=0.5)
kpLines(kp, data=LineData_hapA, ymin=0, ymax=15, y=LineData_hapA$dpi_7, data.panel=1, col="#FF0000", r0=0, r1=0.3, lwd=0.5)
kpLines(kp, data=LineData_hapB, ymin=0, ymax=15, y=LineData_hapB$dpi_1, data.panel=1, col="#0000FF", r0=0, r1=0.3, lwd=0.5)
kpLines(kp, data=LineData_hapB, ymin=0, ymax=15, y=LineData_hapB$dpi_7, data.panel=1, col="#FF0000", r0=0, r1=0.3, lwd=0.5)

#label lpoint of logTPM data
exprData_hapA <- read.table("eff_chr14A_TPM") 
exprData_hapB <- read.table("eff_chr14B_TPM")
exprData_hapA[,5:11] <- log2(exprData_hapA[,5:11]+1)
exprData_hapB[,5:11] <- log2(exprData_hapB[,5:11]+1)
exprData_hapA <-toGRanges(data.frame(chr=exprData_hapA$V1, start=exprData_hapA$V2, end=exprData_hapA$V3, GS=exprData_hapA$V5, dpi_1=exprData_hapA$V6, dpi_4=exprData_hapA$V8, dpi_7=exprData_hapA$V9, dpi_13=exprData_hapA$V11))
exprData_hapB <-toGRanges(data.frame(chr=exprData_hapB$V1, start=exprData_hapB$V2, end=exprData_hapB$V3, GS=exprData_hapB$V5, dpi_1=exprData_hapB$V6, dpi_4=exprData_hapB$V8, dpi_7=exprData_hapB$V9, dpi_13=exprData_hapB$V11))

kpPoints(kp, data=exprData_hapA, ymin=0, ymax=15, y=exprData_hapA$dpi_1, data.panel=1, col="#0000FF", r0=0, r1=0.3, lwd=1,cex=0.5)
kpPoints(kp, data=exprData_hapA, ymin=0, ymax=15, y=exprData_hapA$dpi_7, data.panel=1, col="#FF0000", r0=0, r1=0.3, lwd=1,cex=0.5)

kpPoints(kp, data=exprData_hapB, ymin=0, ymax=15, y=exprData_hapB$dpi_1, data.panel=1, col="#0000FF", r0=0, r1=0.3, lwd=1,cex=0.5)
kpPoints(kp, data=exprData_hapB, ymin=0, ymax=15, y=exprData_hapB$dpi_7, data.panel=1, col="#FF0000", r0=0, r1=0.3, lwd=1,cex=0.5)

# label gene density
#kpAddLabels(kp, labels="gene density", r0=0.35, r1=0.65, cex=0.7)
#kpAxis(kp,r0=0.35,r1=0.65, ymin=0, ymax=1,cex=0.7)
#chr14_gene <- read.table("chr14.gene.bed",sep="\t", header = F)
#gene <- toGRanges(chr14_gene)
#kpPlotDensity(kp, data = gene, window.size=100000, r0 = 0.35, r1 = 0.65 )

#label density secreted protein 
kpAddLabels(kp, labels="Secreted protein", r0=0.35, r1=0.65, cex=0.7)
kpAxis(kp,r0=0.35,r1=0.65, ymin=0, ymax=1,cex=0.7)
secret_protein <- read.table("chr14.secret.bed",sep="\t", header = F)
gene <- toGRanges(secret_protein)
kpPlotDensity(kp, data = gene, window.size=100000, r0 = 0.35, r1 = 0.65 )

#label repeat density
kpAxis(kp,r0=0.70,r1=1.0, ymin=0, ymax=1,cex=0.7)
chr14_repeat <- read.table("chr14.repeat.bed",sep="\t", header = F)
rep <- toGRanges(chr14_repeat)
kpPlotDensity(kp, data = rep, window.size=100000, r0=0.7,r1=1.0)
kpAddLabels(kp, labels = "repeat density", r0=0.70, r1=1.0, cex=0.5)
```
