library(ggplot2)
library(reshape2)
library(amap)
setwd("~/Desktop/RA/data/FIltered_RA")

## deseq2 tables
table1=read.csv("sig_LRT.csv",header=T,row.names = 1)

#normalised count
normalisedCount=read.csv("normalized_counts.csv",header=T,row.names = 1)

#calculating mean
meanTable=apply(normalisedCount,1,mean)
meanTable=as.matrix(meanTable)
colnames(meanTable)=c("mean")

#merging mean with normalised count
exp_table=merge(meanTable,normalisedCount,by=0,row.names=1)
row.names(exp_table)=exp_table$Row.names
exp_table=exp_table[,c(2:170)]

# -log10p calculation
table1$log10p=-log10(table1$pvalue)
names(table1)[names(table1) == "log10p"] <- "-log10p"

##SIGNIFICANCE
table1$significance= as.factor(table1$padj < 0.05 & abs(table1$log2FoldChange) > 1.0)

#merging table1 and expression and renaming columns
merged_1=merge(table1,exp_table,by=0)
row.names(merged_1)=merged_1$Row.names
merged_1=merged_1[,c(2:177)]

#MA plots
ggplot(data=merged_1) + geom_point(aes(x=log10(mean), y= `log2FoldChange`,color=significance))

##density plot
ggplot(merged_1,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")

##volcano plot
ggplot(data=merged_1) + geom_point(aes(x=`log2FoldChange`, y= `-log10p`,color=significance))

#significant ones with log fold change>1 or log fold change < -1
significant_genes=subset(table1,table1$significance==TRUE)
write.csv(significant_genes,"significant_genes.csv")

##
library(tidyverse)
gene_list=significant_genes[,c(7,8)]
gene_list<- gene_list %>% rownames_to_column("Row.names")
gene_list=gene_list[,c(3,1,2)]
names(gene_list)[names(gene_list) == "entrezid"] <- "ENTREZID"
names(gene_list)[names(gene_list) == "Row.names"] <- "ENSEMBL"
names(gene_list)[names(gene_list) == "symbol"] <- "SYMBOL"

#GO Classification
library(clusterProfiler)
gene_1=c("8128","3048","50944","23500","80164","91646","3512","10964","1958","2354","221687","1748","140825","120114","6633","23504","55502","83985","3576","170689","2731","78995")
ggo <- groupGO(gene= gene_1,OrgDb    = 'org.Hs.eg.db',ont= "CC",level = 3,readable = TRUE)

#GO enrichment analysis
library(enrichplot)
edo <- enrichDGN(gene_1)
barplot(edo, showCategory=20)

#Pathway analysis
library(DOSE)
library(enrichplot)
edo <- enrichDGN(gene_1)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, node_label="category") 
p2 <- cnetplot(edox, node_label="gene") 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


#Heatmap
library(pheatmap)
ex=merged_1[,c(7,12:176)]
row.names(ex)=ex$symbol
ex=ex[,c(2:166)]
pheatmap(ex)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(ex, 1, cal_z_score))
pheatmap(data_subset_norm)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(ex, 1, cal_z_score))
pheatmap(data_subset_norm)
data_subset_norm=as.data.frame(data_subset_norm)
data_subset_norm1=data_subset_norm[,c(39:58,1:20,118:137)]
my_sample_col <- data.frame(sample = rep(c("Double", "Single","Triple"), c(20,20,20)))
row.names(my_sample_col) <- colnames(data_subset_norm1)
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(data_subset_norm1, annotation_col = my_sample_col,color = heat_colors)

