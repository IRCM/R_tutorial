library("DESeq2")
library("WriteXLS")
library("gdata")

##################
# DESeq analysis #
##################

# Create output directories
dir.create("../../deseq")
dir.create("../../deseq/all_genes")

sampleTable <- read.table("deseqDesign.txt", header=TRUE)

# Create a DESeqDataSet object.
dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="../../htseqcount", design= ~ condition)

#############################################
# All samples                               #
# Calculate and save the normalized counts. #
#############################################

dds <- DESeq(dds)

normalized.counts <- counts( dds, normalized=TRUE ) 
output <- as.data.frame(normalized.counts)
WriteXLS("output", ExcelFileName = "../../deseq/all_genes/normalized_counts.xls", row.names=TRUE, BoldHeaderRow=TRUE)

# Calculate and save the normalized counts.
normalized.counts <- counts( dds, normalized=TRUE ) 
output <- as.data.frame(normalized.counts)
WriteXLS("output", ExcelFileName = "../../deseq/all_genes/normalized_counts.xls", row.names=TRUE, BoldHeaderRow=TRUE)

# Heatmap of the normalized count table
library("RColorBrewer")
library("gplots")

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30] 
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

dir.create("../../deseq/plots")
pdf("../../deseq/plots/heatmap_30_most_expressed_genes.pdf")
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none",dendrogram="none", trace="none", margin=c(12,12))
dev.off()

colData(dds)$condition <- relevel(colData(dds)$condition, "total")
dds <- DESeq(dds)

##############
# cyto,total #
##############

res <- as.data.frame(results(dds, "condition_cyto_vs_total"))

MeanA <- rowMeans(counts(dds,normalized=TRUE)[,colData(dds)$condition == "cyto"])
MeanB <- rowMeans(counts(dds,normalized=TRUE)[,colData(dds)$condition == "total"])

res <- data.frame(Mean_cyto=MeanA, Mean_total=MeanB, log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, pvalue=res$pvalue, padj=res$padj)
res <- res[order(res$padj),]

WriteXLS("res", ExcelFileName = "../../deseq/all_genes/cyto_vs_total.xls", row.names=TRUE, BoldHeaderRow=TRUE)

