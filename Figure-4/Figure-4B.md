# A Workflow to plot Figure-4B in MER paper

## Step 1 generate bed files for intergenic spaces of all genes, secreted proteins (sp) and busco genes
You can use bedtools to calculate the 5' and 3' intergenic distance, as long as you are aware of the direction of the genes:
  
bioawk -c fastx '{print $name,length($seq)}' hapA.fa > chr.sizes
bedtools complement -i sorted.hapA.cds.bed -g chr.sizes  > hapA_intergenic_spaces.bed
bedtools closest -D a -a sorted.hapA.cds.bed -b hapA_intergenic_spaces.bed -g chr.sizes > hapA_intergenic_spaces.bed2R
grep -f sp.list hapA_intergenic_spaces.bed2R > hapA.sp.bed2R
grep -f busco.list hapA_intergenic_spaces.bed2R > hapA.busco.bed2R

busco.list and sp.list need to be coverted by dos2unix if generated from windows or meet error in grep -f


## Step 2 plot Figure 4B in R
```
setwd("E:/RStudio/effector_space/")
install.packages("hexbin")
  
#load all genes
table <- read.table('hapA_intergenic_spaces.bed2R', header = F)
colnames(table) <- c('Scaffold','Gene_start','Gene_end','Gene_name','Scaffold_inter','Intergen_start','Intergen_stop','UpDown')
table['Inter_size'] <- table$Intergen_stop - table$Intergen_start
table_up <- subset.data.frame(table, UpDown > 0)
table_down <- subset.data.frame(table, UpDown < 0)
table2 <- merge(table_up, table_down, by ='Gene_name')

#load buscos overlay in the middle lane
tableb <- read.table('hapB.busco.bed2R',header = F)
colnames(tableb) <- c('Scaffold','Gene_start','Gene_end','Gene_name','Scaffold_inter','Intergen_start','Intergen_stop','UpDown')
tableb['Inter_size'] <- tableb$Intergen_stop - tableb$Intergen_start
tableb_up <- subset.data.frame(tableb, UpDown > 0)
tableb_down <- subset.data.frame(tableb, UpDown < 0)
table2b <- merge(tableb_up, tableb_down, by ='Gene_name')

#load SPs so you can overlay on top
tablec <- read.table('hapB.sp.bed2R', header = F)
colnames(tablec) <- c('Scaffold','Gene_start','Gene_end','Gene_name','Scaffold_inter','Intergen_start','Intergen_stop','UpDown')
tablec['Inter_size'] <- tablec$Intergen_stop - tablec$Intergen_start
tablec_up <- subset.data.frame(tablec, UpDown > 0)
tablec_down <- subset.data.frame(tablec, UpDown < 0)
table2c <- merge(tablec_up, tablec_down, by ='Gene_name')

#Annotate secreted and non secreted
table2$group <- "Non SP"
table2$group <- "Non BUSCO"
table2$group[table2$Gene_name %in% table2c$Gene_name] <- "SP"
table2$group[table2$Gene_name %in% table2b$Gene_name] <- "BUSCO"
ggplot(table2, aes(x= Inter_size.x,y=Inter_size.y)) +
  geom_hex()+
  scale_fill_distiller(palette ="RdBu", na.value = "grey50", direction = -1,guide = "colourbar", name ="Gene number")+
  geom_point(data=subset(table2, group == "BUSCO"), size = 0.5, color = 'yellow', alpha=0.5)+
  geom_point(data=subset(table2, group == "SP"), size = 0.5, color = 'red', alpha=0.5)+
  geom_point(data=subset(table2, Gene_name == "FUN_023436-T1"),size = 1, color = 'black')+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' intergenic region (bp)") + xlab("3' intergenic region (bp)")+theme(axis.title = element_text(size = 10),axis.text = element_text(size=10, color="black"),legend.title=element_text(size = 8, color = 'black'), legend.text=element_text(size=8,color="black"))
```
