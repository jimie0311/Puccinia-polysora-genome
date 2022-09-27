# Plot rd in Figure 7C
## Install R in Linux 
Because the vcf file is super large, we suggest to run R in linux
```
conda install R
R()
```
## calculate rD using vcfR
```
##########################
####     LD      #########
##########################

library("vcfR") #1.8.0
library("poppr") #2.7.1
library("ggplot2")
library("agricolae")
#save vcf file to Rdata
vcfFile <- read.vcfR("/PATH/TO/filtered.vcf")
vcfFile <- vcfFile[is.biallelic(vcfFile)] # Our sample is diploid and only biallelic loci are used
save(vcfFile, file = "Data/vcfFile.Rdata")

#transfer vcfFile to glFile
glFile <- vcfR2genlight(vcfFile)
save(glFile,file = "glFile.Rdata") #if the vcf file is large, please save it for next use
load("glFile.Rdata")

#set ploidy and remove NAs
ploidy(glFile) <- 2
toRemove <- is.na(glMean(glFile, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove)
glFile_rmNA <- glFile[, !toRemove] 
save(glFile_rmNA, file = "glFile_rmNA.Rdata")
load("Data/glFile_rmNA.Rdata")

# calculate IA of our population
set.seed(49)
subsnp <- samp.ia(glFile_rmNA, n.snp = 10000L, reps = 100L)

#perform simulation
###no structure divergence, fully admixed pops, in total of 394886 biallelic loci were left in our case and 77 isolates were used.
sex <- glSim(77, 394886, ploid=2, LD=T, parallel = TRUE)
### Semi-clonal 
clone_50 <- glSim(77, 394886, n.snp.struc=197443, ploid=2, LD=T, parallel = TRUE) #394886*0.5
### Most-clonal 
clone_75 <- glSim(77, 394886, n.snp.struc=296164, ploid=2, LD=T, parallel = TRUE) #394886*0.75
### Structure (clonal pops)
clone_100 <- glSim(77, 394886, n.snp.struc=394886, ploid=2, LD = T, parallel = TRUE)

## calculate IA
ia.sex <- samp.ia(sex, n.snp = 10000L, reps = 100L, threads = 4L)
## IA clone
ia.clone.50 <- samp.ia(clone_50, n.snp = 10000L, reps = 100L, threads = 4L)
## IA.semiclone
ia.clone.75 <- samp.ia(clone_75, n.snp = 10000L, reps = 100L, threads = 4L)
## IA.mostclone
ia.clone.100 <- samp.ia(clone_100, n.snp = 10000L, reps = 100L, threads = 4L)

# Summarizing data frames
d1 <- data.frame(subsnp, rep("all", length(subsnp)))
d2 <- data.frame(ia.sex, rep("sexual", length(ia.sex)))
d3 <- data.frame(ia.clone.50, rep("clone_50", length(ia.clone.50)))
d4 <- data.frame(ia.clone.75, rep("clone_75", length(ia.clone.75)))
d5 <- data.frame(ia.clone.100, rep("clone_100", length(ia.clone.100)))
colnames(d1) <- c("ia","dset")
colnames(d2) <- c("ia","dset")
colnames(d3) <- c("ia","dset")
colnames(d4) <- c("ia","dset")
colnames(d5) <- c("ia","dset")


# Normality test
# The Shapiro-Wilk’s normality test will show that whether the following observed and simulated distributions follow normality, if not follow, do simulation again
frames <- list(as.data.frame(d1), as.data.frame(d2), as.data.frame(d3), as.data.frame(d4), as.data.frame(d5))
normality <- list()
for (i in 1:length(frames)){
  normality[[i]] <- shapiro.test(frames[[i]][,'ia'])
}

# Kruskal wallis test
kruskal.test(ia ~ dset, ia.total) #Kruskal-Wallis chi-squared test
k.test <- with(ia.total, kruskal(ia, dset, group = T, p.adj = "bon"))
```
## plot distributions in Rstudio
#move ia.total.Rdata to Rstudio
```
setwd("E:/RStudio/rd/")
load("ia.total.Rdata")
library("ggplot2")
#Set theme to be classic for rest of plots
theme_set(theme_classic() + 
            theme(
              legend.text=element_text(size=11),
              axis.text=element_text(size=10),
              axis.title = element_text(size=15),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 13)
            )
)
ia.total$dset <- factor(ia.total$dset, levels = c("all", "sexual", "clone_50", "clone_75", "clone_100"))#将分组信息转化为因子，否则画图时候会横坐标顺序会乱
byTotal_plot <- ggplot(ia.total, aes(x = dset, y = log10(ia), color = as.factor(dset))) +
  geom_boxplot() +
  scale_x_discrete(labels = c("observed", "0%", "50%", "75%", "100%")) +
  #  scale_y_continuous(limits = c(0, 0.075)) +
  scale_color_manual(values = c("red",  rep("black", 4))) +
  theme(legend.position="none") +
  labs(x = "Linkage level", y = expression(italic(bar("r")["d"])))
byTotal_plot
```
