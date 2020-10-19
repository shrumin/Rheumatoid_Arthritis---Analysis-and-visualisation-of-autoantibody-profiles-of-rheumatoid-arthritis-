library(ggplot2)
library(reshape2)
library(amap)
setwd("~/Desktop/RA/data/FIltered_RA")

## deseq2 tables

table1=read.csv("sig_CCP_Acet_vs_CCP.csv",header=T,row.names = 1)
table1=table1[,c(2,5:8)]
table2=read.csv("sig_CCP_Carb_vs_CCP.csv",header=T,row.names = 1)
table2=table2[,c(2,5:8)]
table3=read.csv("sig_CCP_Carb_Acet_vs_CCP.csv",header=T,row.names = 1)
table3=table3[,c(2,5:8)]
table4=read.csv("sig_Negative_vs_CCP.csv",header=T,row.names = 1)
table4=table4[,c(2,5:8)]

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
table2$log10p=-log10(table2$pvalue)
names(table2)[names(table2) == "log10p"] <- "-log10p"
table3$log10p=-log10(table3$pvalue)
names(table3)[names(table3) == "log10p"] <- "-log10p"
table4$log10p=-log10(table4$pvalue)
names(table4)[names(table4) == "log10p"] <- "-log10p"

##SIGNIFICANCE
table1$significance= as.factor(table1$padj < 0.05 & abs(table1$log2FoldChange) > 1.0)
table2$significance= as.factor(table2$padj < 0.05 & abs(table2$log2FoldChange) > 1.0)
table3$significance= as.factor(table3$padj < 0.05 & abs(table3$log2FoldChange) > 1.0)
table4$significance= as.factor(table4$padj < 0.05 & abs(table4$log2FoldChange) > 1.0)

#merging all tables for master table

#merging table1 and expression and renaming columns

merged_1=merge(table1,exp_table,by=0)
row.names(merged_1)=merged_1$Row.names
merged_1=merged_1[,c(2:177)]

#merging table2 and expression and renaming columns

merged_2=merge(table2,exp_table,by=0)
row.names(merged_2)=merged_2$Row.names
merged_2=merged_2[,c(2:177)]

#merging table3 and expression and renaming columns

merged_3=merge(table3,exp_table,by=0)
row.names(merged_3)=merged_3$Row.names
merged_3=merged_3[,c(2:177)]

#merging table4 and expression and renaming columns

merged_4=merge(table4,exp_table,by=0)
row.names(merged_4)=merged_4$Row.names
merged_4=merged_4[,c(2:177)]



#MA plots
ggplot(data=merged_1) + geom_point(aes(x=log10(mean), y= `log2FoldChange`,color=significance))
ggplot(data=merged_2) + geom_point(aes(x=log10(mean), y= `log2FoldChange`,color=significance))
ggplot(data=merged_3) + geom_point(aes(x=log10(mean), y= `log2FoldChange`,color=significance))
ggplot(data=merged_4) + geom_point(aes(x=log10(mean), y= `log2FoldChange`,color=significance))

##density plot
ggplot(merged_1,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")
ggplot(merged_2,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")
ggplot(merged_3,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")
ggplot(merged_4,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")

##volcano plot
ggplot(data=merged_1) + geom_point(aes(x=`log2FoldChange`, y= `-log10p`,color=significance))
ggplot(data=merged_2) + geom_point(aes(x=`log2FoldChange`, y= `-log10p`,color=significance))
ggplot(data=merged_3) + geom_point(aes(x=`log2FoldChange`, y= `-log10p`,color=significance))
ggplot(data=merged_4) + geom_point(aes(x=`log2FoldChange`, y= `-log10p`,color=significance))
#ggplot(data=master_table) + geom_point(aes(x=`log2FoldChange_Negative_vs_CCP`, y= `-log10p_Negative_vs_CCP`,color=significance_Negative_vs_CCP))

#significant ones with log fold change>1 or log fold change < -1
significant_genes_CCP_Acet_vs_CCP=subset(table1,table1$significance==TRUE)
significant_genes_CCP_Carb_vs_CCP=subset(table2,table2$significance==TRUE)
significant_genes_CCP_Carb_Acet_vs_CCP=subset(table3,table3$significance==TRUE)
significant_genes_Negative_vs_CCP=subset(table4,table4$significance==TRUE)
#no significant genes in CCP_Carb_vs _CCP , 2 significant genes in CCP_Acet_vs_CCP, 10 
#significant genes each in CCP_Carb_Acet_vs_CCP and Negative_vs_CCP
write.csv(significant_genes_CCP_Acet_vs_CCP,"significant_genes_CCP_Acet_vs_CCP.csv")
write.csv(significant_genes_CCP_Carb_Acet_vs_CCP,"significant_genes_CCP_Carb_Acet_vs_CCP.csv")
write.csv(significant_genes_Negative_vs_CCP,"significant_genes_Negative_vs_CCP.csv")


#merging only the significant genes
merge_demo=rbind(merged_1,merged_3)
merge_demo=rbind(merge_demo,merged_4)
merge_demo=unique(merge_demo)
#rownames(merge_demo)=merge_demo$Row.names
#merge_demo=merge_demo[,c(9:176)]
#rm(merged_1,merged_2,merged_3)
merged_1=subset(merged_1,merged_1$significance==TRUE)
merged_2=subset(merged_2,merged_2$significance==TRUE)
merged_3=subset(merged_3,merged_3$significance==TRUE)
merged_4=subset(merged_4,merged_4$significance==TRUE)
#10 significant genes in CCP_Carb_Acet_vs_CCP and Negative_vs_CCP are same

##
library(tidyverse)
gene_list=merge_demo[,c(4,5)]
gene_list<- gene_list %>% rownames_to_column("Row.names")
gene_list=gene_list[,c(3,1,2)]
names(gene_list)[names(gene_list) == "entrezid"] <- "ENTREZID"
names(gene_list)[names(gene_list) == "Row.names"] <- "ENSEMBL"
names(gene_list)[names(gene_list) == "symbol"] <- "SYMBOL"

#GO Classification
library(clusterProfiler)
#GO Classification
library(clusterProfiler)
gene_1=c("3757","23504","1748","84417","140825","2354","23500","357","50944","91646","3048")
ggo <- groupGO(gene= gene_1,OrgDb    = 'org.Hs.eg.db',ont= "CC",level = 3,readable = TRUE)

#GO enrichment analysis
library(enrichplot)
edo <- enrichDGN(gene_1)
barplot(edo, showCategory=20)

#Heatmap

#CCP_Acet vs CCP
library(pheatmap)
#ex=merged_1[,c(9:176)]
ex=merged_1[,c(4,9:176)]
row.names(ex)=ex$symbol
ex=ex[,c(2:169)]
pheatmap(ex)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(ex, 1, cal_z_score))
pheatmap(data_subset_norm)
data_subset_norm=as.data.frame(data_subset_norm)
data_subset_norm1=data_subset_norm[,c(39:58,1:20)]
my_sample_col <- data.frame(sample = rep(c("CCP_Acet", "CCP"), c(20,20)))
row.names(my_sample_col) <- colnames(data_subset_norm1)
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(data_subset_norm1, annotation_col = my_sample_col,color = heat_colors)


#Triple vs CCP
library(pheatmap)
ex=merged_3[,c(4,9:176)]
row.names(ex)=ex$symbol
ex=ex[,c(2:169)]
pheatmap(ex)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(ex, 1, cal_z_score))
pheatmap(data_subset_norm)
data_subset_norm=as.data.frame(data_subset_norm)
data_subset_norm1=data_subset_norm[,c(118:137,1:20)]
my_sample_col <- data.frame(sample = rep(c("Triple", "CCP"), c(20,20)))
row.names(my_sample_col) <- colnames(data_subset_norm1)
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(data_subset_norm1, annotation_col = my_sample_col,color = heat_colors)

#negative vs CCP
library(pheatmap)
ex=merged_4[,c(4,9:176)]
row.names(ex)=ex$symbol
ex=ex[,c(2:169)]
pheatmap(ex)
cal_z_score <- function(x){(x - mean(x)) / sd(x)}
data_subset_norm <- t(apply(ex, 1, cal_z_score))
pheatmap(data_subset_norm)
data_subset_norm=as.data.frame(data_subset_norm)
data_subset_norm1=data_subset_norm[,c(139:168,1:38)]
my_sample_col <- data.frame(sample = rep(c("Negative", "CCP"), c(38,30)))
row.names(my_sample_col) <- colnames(data_subset_norm1)
library(RColorBrewer)
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(data_subset_norm1, annotation_col = my_sample_col,color = heat_colors)


