setwd("C:/Rstudio/Projects/Puccinia/snp/")
library("vcfR")
library("ggplot2")

#read VCF file
vcfFile <- read.vcfR("all_sample.vcf")
save(vcfFile, file="all_sample.Rdata")

#Subset to het positions
#This subsetting will place an NA in non-heterozygous positions in the vcfR object
vcfFile_gt <- extract.gt(vcfFile, element = 'GT')
vcfFile_hets <- is_het(vcfFile_gt)
is.na(vcfFile@gt[,-1][ !vcfFile_hets ]) <- TRUE

#Filter het positions on allele depth
vcfFile_ad <- extract.gt(vcfFile, element = 'AD')

vcfFile_allele1 <- masplit(vcfFile_ad, record = 1)
vcfFile_allele2 <- masplit(vcfFile_ad, record = 2)

#Produce quantiles at 0.15 and 0.95 probabilites
sumsFile <- apply(vcfFile_allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
?apply
#Allele 1
vcfFile_dp2 <- sweep(vcfFile_allele1, MARGIN=2, FUN = "-", sumsFile[1,])
?sweep
#allele1[dp2 < 0] <- NA
vcfFile@gt[,-1][ vcfFile_dp2 < 0 & !is.na(vcfFile@gt[,-1]) ] <- NA
vcfFile_dp2 <- sweep(vcfFile_allele1, MARGIN=2, FUN = "-", sumsFile[2,])

#allele1[dp2 > 0] <- NA
vcfFile@gt[,-1][vcfFile_dp2 > 0] <- NA

#Allele 2
vcfFile_dp2 <- sweep(vcfFile_allele2, MARGIN=2, FUN = "-", sumsFile[1,])
vcfFile@gt[,-1][ vcfFile_dp2 < 0 & !is.na(vcfFile@gt[,-1]) ] <- NA
vcfFile_dp2 <- sweep(vcfFile_allele2, MARGIN=2, FUN = "-", sumsFile[2,])
vcfFile@gt[,-1][vcfFile_dp2 > 0] <- NA

#Calculate allele balance with the new filtered vcf
vcfFile_ad <- extract.gt(vcfFile, element = 'AD')

vcfFile_allele1 <- masplit(vcfFile_ad, record = 1)
vcfFile_allele2 <- masplit(vcfFile_ad, record = 2)


vcfFile_ad1 <- vcfFile_allele1 / (vcfFile_allele1 + vcfFile_allele2)
vcfFile_ad2 <- vcfFile_allele2 / (vcfFile_allele1 + vcfFile_allele2)

#Make a tidy dataframe to make plot
vcfFile_ad1t <- tidyr::gather(tibble::as.tibble(vcfFile_ad1), "Sample", "Ab")
vcfFile_ad1t$Allele <- "ab1"
?tidyr
vcfFile_ad2t <- tidyr::gather(tibble::as.tibble(vcfFile_ad2), "Sample", "Ab")
vcfFile_ad2t$Allele <- "ab2"

vcfFile_t <- rbind(vcfFile_ad1t, vcfFile_ad2t)
vcfFile_t <- dplyr::filter(vcfFile_t, !is.na(Ab))

#Make plot
pFile <- ggplot(vcfFile_t, aes(x=Ab))
pFile <- pFile + scale_x_continuous(breaks=c(0, 1/2, 1), labels = c('0', '1/2', '1'))
pFile <- pFile + geom_histogram(data=subset(vcfFile_t, Allele == 'ab1'),fill =  "#1f78b4", binwidth = 0.02)
pFile <- pFile + geom_histogram(data=subset(vcfFile_t, Allele == 'ab2'),fill =  "#a6cee3", binwidth = 0.02)
pFile <- pFile + facet_wrap( ~ Sample, ncol=10, nrow = 9)
pFile <- pFile + theme_bw()

pFile <- pFile + theme(panel.grid.major.x=element_line(color = "#A9A9A9", size=0.4) )
pFile <- pFile + theme(panel.grid.major.y=element_line(color = "#A9A9A9", size=0.1, linetype = 'dashed') )
pFile <- pFile + theme(panel.grid.minor.x=element_blank())
pFile <- pFile + theme(panel.grid.minor.y=element_blank())
pFile <- pFile + theme(strip.text.x = element_text(size = 8))
pFile <- pFile + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=8))
pFile <- pFile + theme(axis.text.y = element_text(angle = 30, hjust = 1, size=10))

pFile <- pFile + theme(axis.title.x = element_text(size = 14, angle = 0))
pFile <- pFile + theme(axis.title.y = element_text(size = 14, angle = 90))
pFile <- pFile + xlab("Allele balance")
pFile <- pFile + ylab("Heterozygous positions")

#Save the plot
pdf("allele_balanceFile.pdf", height = 9, width = 9)
pFile
dev.off()
rlang::last_error()

