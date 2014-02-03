library("biomaRt")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="ensembl.org")

gene.biotypes=getBM(attributes=c('ensembl_gene_id', 'external_gene_id', 'gene_biotype'), mart=ensembl)

write.table(gene.biotypes, "annotation.txt", sep="\t", row.names=FALSE, quote=TRUE)
