setwd("~/Desktop/RA/data/FIltered_RA")
dat <- read.csv(file="new_dat.csv", row.names=1)

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

col=read.csv("exp_info_triple_vs_negative.csv",row.names = 1)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=dat,colData=col,design= ~ group)
keep <- rowSums(counts(dds) >=10) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

#Group negative vs CCP_Carb_Acet
res_Negative_CCP_Carb_Acet <- results(dds,contrast=c("group","Negative","CCP_Carb_Acet"))
res_Negative_vs_CCP_Carb_Acet_sort=res_Negative_CCP_Carb_Acet[order(res_Negative_CCP_Carb_Acet$padj),]
res_Negative_vs_CCP_Carb_Acet_sig<- subset(res_Negative_vs_CCP_Carb_Acet_sort,padj <0.05)
res_Negative_vs_CCP_Carb_Acet_sig$symbol <- t2g$external_gene_name[match(rownames(res_Negative_vs_CCP_Carb_Acet_sig), t2g$ensembl_gene_id)]
res_Negative_vs_CCP_Carb_Acet_sig$entrezid <- t2g$entrezgene[match(rownames(res_Negative_vs_CCP_Carb_Acet_sig), t2g$ensembl_gene_id)]
res_Negative_vs_CCP_Carb_Acet_sig=na.omit(res_Negative_vs_CCP_Carb_Acet_sig)
write.csv(res_Negative_vs_CCP_Carb_Acet_sig, "sig_Negative_vs_CCP_Carb_Acet.csv", row.names=TRUE)

#normalised count
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts-2.csv", quote=F, col.names=NA)

#immune system genes
immunegenes=read.csv("InnateDB_genes.csv",row.names = 4)
immunegenes=immunegenes[c(5)]

Immunereults_Negative_vs_CCP_Carb_Acet=merge(as.data.frame(res_Negative_vs_CCP_Carb_Acet_sig),immunegenes,by=0)
row.names(Immunereults_Negative_vs_CCP_Carb_Acet)=Immunereults_Negative_vs_CCP_Carb_Acet$Row.names
Immunereults_Negative_vs_CCP_Carb_Acet=Immunereults_Negative_vs_CCP_Carb_Acet[c(2:9)]
write.csv(Immunereults_Negative_vs_CCP_Carb_Acet,"Immunereults_Negative_vs_CCP_Carb_Acet.csv",row.names = TRUE)
