setwd("~/Desktop/RA")

dat <- read.csv(file="Joint_Narrow_dat.csv", row.names=1)

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

#design exp table

colTable <- read.csv("Joint_Narrowing_exp_info.csv",row.names=1)

#Deseq analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=dat,colData=colTable,design= ~ group)
keep <- rowSums(counts(dds) >=10) >= 10
dds <- dds[keep,]
dds <- DESeq(dds,test="LRT",reduced = ~1)
resultsNames(dds)

#results
res_Joint_Narrowing<- results(dds,contrast=c("group","No_Joint_Narrowing","Joint_Narrowing"))
res_Joint_Narrowing_sort=res_Joint_Narrowing[order(res_Joint_Narrowing$padj),]
res_Joint_Narrowing_sig<- subset(res_Joint_Narrowing_sort,padj<0.05)
res_Joint_Narrowing_sig$symbol <- t2g$external_gene_name[match(rownames(res_Joint_Narrowing_sig), t2g$ensembl_gene_id)]
res_Joint_Narrowing_sig$entrezid <- t2g$entrezgene[match(rownames(res_Joint_Narrowing_sig), t2g$ensembl_gene_id)]
res_Joint_Narrowing_sig=na.omit(res_Joint_Narrowing_sig)
write.csv(res_Joint_Narrowing_sig, "sig_No_Joint_Narrowing_vs_Joint_Narrowing.csv", row.names=TRUE)

#normalised count
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts_Joint_Narrowing.csv", quote=F, col.names=NA)

#plot count of HLA-DQB1 gene
plotCounts(dds, gene="ENSG00000179344", intgroup="group")
plotCounts(dds, gene="ENSG00000179344", intgroup="group", transform=FALSE)

