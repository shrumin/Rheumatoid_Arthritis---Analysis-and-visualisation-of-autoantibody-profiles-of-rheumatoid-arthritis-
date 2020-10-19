setwd("~/Desktop/RA/data/FIltered_RA")
dat <- read.csv(file="dat8.csv", row.names=1)

#filtering and naming genes
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
nonprotein=biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "entrezgene_id", "description", "gene_biotype"), filters='biotype', values=c("rRNA"), mart = ensembl)
#All protein coding genes
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name", "entrezgene_id", "description", "gene_biotype"), filters='biotype', values="protein_coding", mart = ensembl)

#Ribosomal proteins
ribo.proteins <- unique(t2g$ensembl_gene_id[grep("ribosomal protein", t2g$description)])


rRNA.proteins <- unique(nonprotein$ensembl_gene_id[grep("rRNA", nonprotein$biotype)])

#remove version from identifier
rownames(dat) <- sub("\\.\\d+", "", rownames(dat))

#remove any ribosomal proteins from the RNA-Seq from the counts
dat <- dat[!rownames(dat) %in% ribo.proteins,]
dat <- dat[!rownames(dat) %in% rRNA.proteins,]


#rm(dat1,dat2,dat3,dat4,dat5,dat6,dat7,myvars1,myvars2,myvars3,myvars4,myvars5,myvars6,myvars7,myvars8)
#rmdesign exp table
#colTable <- read.csv("Design-exp.csv",row.names=1)

#design exp table without acetylated

colTable <- read.csv("final_experiment_info.csv",row.names=1)


#Deseq analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=dat,colData=colTable,design= ~ group)
keep <- rowSums(counts(dds) >=10) >= 10
dds <- dds[keep,]
dds <- DESeq(dds, test="LRT",reduced = ~1)
resultsNames(dds)


#saving results

#CCP_Acet vs CCP

res_CCP_Acet_CCP <- results(dds,contrast=c("group","CCP_Acet","CCP"))
res_CCP_Acet_vs_CCP_sort=res_CCP_Acet_CCP[order(res_CCP_Acet_CCP$padj),]
res_CCP_Acet_vs_CCP_sig<- subset(res_CCP_Acet_vs_CCP_sort,padj<0.05)
res_CCP_Acet_vs_CCP_sig$symbol <- t2g$external_gene_name[match(rownames(res_CCP_Acet_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Acet_vs_CCP_sig$entrezid <- t2g$entrezgene[match(rownames(res_CCP_Acet_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Acet_vs_CCP_sig=na.omit(res_CCP_Acet_vs_CCP_sig)
write.csv(res_CCP_Acet_vs_CCP_sig, "sig_CCP_Acet_vs_CCP.csv", row.names=TRUE)


#CCP_Carb vs CCP

res_CCP_Carb_CCP <- results(dds,contrast=c("group","CCP_Carb","CCP"))
res_CCP_Carb_vs_CCP_sort=res_CCP_Carb_CCP[order(res_CCP_Carb_CCP$padj),]
res_CCP_Carb_vs_CCP_sig<- subset(res_CCP_Carb_vs_CCP_sort,padj<0.05)
res_CCP_Carb_vs_CCP_sig$symbol <- t2g$external_gene_name[match(rownames(res_CCP_Carb_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Carb_vs_CCP_sig$entrezid <- t2g$entrezgene[match(rownames(res_CCP_Carb_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Carb_vs_CCP_sig=na.omit(res_CCP_Carb_vs_CCP_sig)
write.csv(res_CCP_Carb_vs_CCP_sig, "sig_CCP_Carb_vs_CCP.csv", row.names=TRUE)

# CCP_Carb_Acet vs CCP

res_CCP_Carb_Acet_CCP <- results(dds,contrast=c("group","CCP_Carb_Acet","CCP"))
res_CCP_Carb_Acet_vs_CCP_sort=res_CCP_Carb_Acet_CCP[order(res_CCP_Carb_Acet_CCP$padj),]
res_CCP_Carb_Acet_vs_CCP_sig<- subset(res_CCP_Carb_Acet_vs_CCP_sort,padj<0.05)
res_CCP_Carb_Acet_vs_CCP_sig$symbol <- t2g$external_gene_name[match(rownames(res_CCP_Carb_Acet_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Carb_Acet_vs_CCP_sig$entrezid <- t2g$entrezgene[match(rownames(res_CCP_Carb_Acet_vs_CCP_sig), t2g$ensembl_gene_id)]
res_CCP_Carb_Acet_vs_CCP_sig=na.omit(res_CCP_Carb_Acet_vs_CCP_sig)
write.csv(res_CCP_Carb_Acet_vs_CCP_sig, "sig_CCP_Carb_Acet_vs_CCP.csv", row.names=TRUE)

#Group negative vs CCP

res_Negative_CCP <- results(dds,contrast=c("group","Negative","CCP"))
res_Negative_vs_CCP_sort=res_Negative_CCP[order(res_Negative_CCP$padj),]
res_Negative_vs_CCP_sig<- subset(res_Negative_vs_CCP_sort, padj<0.05)
res_Negative_vs_CCP_sig$symbol <- t2g$external_gene_name[match(rownames(res_Negative_vs_CCP_sig), t2g$ensembl_gene_id)]
res_Negative_vs_CCP_sig$entrezid <- t2g$entrezgene[match(rownames(res_Negative_vs_CCP_sig), t2g$ensembl_gene_id)]
res_Negative_vs_CCP_sig=na.omit(res_Negative_vs_CCP_sig)
write.csv(res_Negative_vs_CCP_sig, "sig_Negative_vs_CCP.csv", row.names=TRUE)

#immune system genes
immunegenes=read.csv("InnateDB_genes.csv",row.names = 4)
immunegenes=immunegenes[c(5)]

#immune genes in significant results of each combination
Immunereults_CCP_Acet_vs_CCP=merge(as.data.frame(res_CCP_Acet_vs_CCP_sig),immunegenes,by=0)
row.names(Immunereults_CCP_Acet_vs_CCP)=Immunereults_CCP_Acet_vs_CCP$Row.names
Immunereults_CCP_Acet_vs_CCP=Immunereults_CCP_Acet_vs_CCP[c(2:9)]
write.csv(Immunereults_CCP_Acet_vs_CCP,"Immunereults_CCP_Acet_vs_CCP.csv",row.names = TRUE)

Immunereults_CCP_Carb_vs_CCP=merge(as.data.frame(res_CCP_Carb_vs_CCP_sig),immunegenes,by=0)
row.names(Immunereults_CCP_Carb_vs_CCP)=Immunereults_CCP_Carb_vs_CCP$Row.names
Immunereults_CCP_Carb_vs_CCP=Immunereults_CCP_Carb_vs_CCP[c(2:9)]
write.csv(Immunereults_CCP_Carb_vs_CCP,"Immunereults_CCP_carb_vs_CCP.csv",row.names = TRUE)

Immunereults_CCP_Carb_Acet_vs_CCP=merge(as.data.frame(res_CCP_Carb_Acet_vs_CCP_sig),immunegenes,by=0)
row.names(Immunereults_CCP_Carb_Acet_vs_CCP)=Immunereults_CCP_Carb_Acet_vs_CCP$Row.names
Immunereults_CCP_Carb_Acet_vs_CCP=Immunereults_CCP_Carb_Acet_vs_CCP[c(2:9)]
write.csv(Immunereults_CCP_Carb_Acet_vs_CCP,"Immunereults_CCP_carb_Acet_vs_CCP.csv",row.names = TRUE)

Immunereults_Negative_vs_CCP=merge(as.data.frame(res_Negative_vs_CCP_sig),immunegenes,by=0)
row.names(Immunereults_Negative_vs_CCP)=Immunereults_Negative_vs_CCP$Row.names
Immunereults_Negative_vs_CCP=Immunereults_Negative_vs_CCP[c(2:9)]
write.csv(Immunereults_Negative_vs_CCP,"Immunereults_Negative_vs_CCP.csv",row.names = TRUE)

#normalised count
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts.csv", quote=F, col.names=NA)


