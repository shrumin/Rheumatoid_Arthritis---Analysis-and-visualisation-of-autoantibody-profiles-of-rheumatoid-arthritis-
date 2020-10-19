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

#Result Extraction
res_LRT <- results(dds)
res_sort=res_LRT[order(res_LRT$padj),]
res_LRT_sig<- subset(res_sort,padj<0.05)
res_LRT_sig$symbol <- t2g$external_gene_name[match(rownames(res_LRT_sig), t2g$ensembl_gene_id)]
res_LRT_sig$entrezid <- t2g$entrezgene[match(rownames(res_LRT_sig), t2g$ensembl_gene_id)]
res_LRT_sig=na.omit(res_LRT_sig)
write.csv(res_LRT_sig, "sig_LRT.csv", row.names=TRUE)

#normalised count
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts.csv", quote=F, col.names=NA)
