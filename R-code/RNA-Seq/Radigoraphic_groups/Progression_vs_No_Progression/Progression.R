setwd("~/Desktop/RA")

dat <- read.csv(file="Progression_data.csv", row.names=1)

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

#design exp table without acetylated

colTable <- read.csv("Progression_exp_info.csv",row.names=1)


#Deseq analysis

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=dat,colData=colTable,design= ~ group)
keep <- rowSums(counts(dds) >=10) >= 10
dds <- dds[keep,]
dds <- DESeq(dds,test="LRT",reduced = ~1)
resultsNames(dds)

#results
res_Progression<- results(dds,contrast=c("group","Progession","No_Progession"))
res_Progression_sort=res_Progression[order(res_Progression$padj),]
res_Progression_sig<- subset(res_Progression_sort,padj<0.05)
res_Progression_sig$symbol <- t2g$external_gene_name[match(rownames(res_Progression_sig), t2g$ensembl_gene_id)]
res_Progression_sig$entrezid <- t2g$entrezgene[match(rownames(res_Progression_sig), t2g$ensembl_gene_id)]
res_Progression_sig=na.omit(res_Progression_sig)
write.csv(res_Progression_sig, "sig_Progression_vs_No_Progression.csv", row.names=TRUE)

#normalised count
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, file="normalized_counts_Progression.csv", quote=F, col.names=NA)


