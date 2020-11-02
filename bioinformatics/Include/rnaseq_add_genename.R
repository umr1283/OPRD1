#!/usr/bin/Rscript

library(biomaRt)

argv <- commandArgs(TRUE)
samplename <- argv[1]
host <- argv[2]
dataset <- argv[3]

ensembl <- useMart(host=host,biomart="ENSEMBL_MART_ENSEMBL",dataset=dataset)
#listDatasets(ensembl)
genefile <- paste(samplename,".genes.results", sep = "")
transcriptfile <- paste(samplename,".isoforms.results", sep = "")
geneout <- paste(samplename,".genes_genename.results", sep = "")
transcriptout <- paste(samplename,".isoforms_genename.results", sep = "")


counting <- read.delim(genefile, header = TRUE)
genemap<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values=counting$gene_id,mart=ensembl)
idx<-match(counting$gene_id,genemap$ensembl_gene_id)
counting$external_gene_name<-genemap$external_gene_name[idx]
write.csv(as.data.frame(counting),file=geneout)

counting <- read.delim(transcriptfile, header = TRUE)
genemap<-getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","external_gene_name"),filters="ensembl_transcript_id",values=counting$transcript_id,mart=ensembl)
idx<-match(counting$transcript_id,genemap$ensembl_transcript_id)
counting$external_gene_name<-genemap$external_gene_name[idx]
write.csv(as.data.frame(counting),file=transcriptout)
