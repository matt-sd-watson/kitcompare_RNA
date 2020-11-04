library(Rsubread)
library(DESeq2)
library(data.table)
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("genefilter")
library(EDASeq)
library(ggplot2)
library(rmarkdown)
library(geneplotter)
library(plyr)
library(LSD)
library(gplots)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(fdrtool)
library(org.Mm.eg.db)
library(BiocStyle)
library(org.Mm.eg.db)

## note that this code continues from featureCounts_rna_type_coverage
## and uses various objects that were created in that environment

kits <- c(rep("Kit A", 6), rep(c("Kit B", "Kit C"),6))
col_data <- data.frame(sample = sample_names, kit = kits)
all_norm_data <- rbind(tpms_kit_a, tpms_kit_b, tpms_kit_c)
norm_data_for_deseq <- all_norm_data[,1:3]
norm_data_for_deseq <- spread(norm_data_for_deseq, key = sample, value = counts)
norm_data_for_deseq <- norm_data_for_deseq[,-1]
row.names(norm_data_for_deseq) <- unique(all_norm_data$gene)
norm_data_for_deseq <- norm_data_for_deseq[complete.cases(norm_data_for_deseq),]
deseqdata <- DESeqDataSetFromMatrix(countData = norm_data_for_deseq, colData = col_data, design = ~kit)
nrow(deseqdata)


differential_expression <- DESeq(deseqdata)
summary(differential_expression)
# view the comparisons between the kits that are determined by deseq2
resultsNames(differential_expression)

vsd <- vst(differential_expression)
head(assay(vsd), 10)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

poisd <- PoissonDistance(t(counts(deseqdata)))
pois_matrix <- as.matrix(poisd$dd)
rownames(pois_matrix) <- col_data$sample
colnames(pois_matrix) <- col_data$sample
pheatmap(pois_matrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


pca_plot <- plotPCA(vsd, intgroup = "kit")
pca_plot + geom_text(aes(label = name), position = position_nudge(y = -3), size = 2.5)

res <- results(differential_expression)
res

ordered_p_value <- res[order(res$pvalue),]
summary(ordered_p_value)


plotMA(res, ylim=c(min(res$log2FoldChange),max(res$log2FoldChange)), alpha = 0.05)

# plot ma with padj less than the threshold of 0.05
plotMA(res[res$padj > 0.05 & !is.na(res$padj),])

topgenes <- row.names(ordered_p_value[1:10,])
with(res[topgenes, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2)
  text(baseMean, log2FoldChange, topgenes, pos=2, col="red")
})

topGene <- rownames(res)[which.min(res$padj)]
topGene
# get the top 10 genes by ordered p value, not using padj
top10genes <- row.names(ordered_p_value[1:10,])

column_data <- data.frame(kit = col_data$kit)
row.names(column_data) <- col_data$sample
column_data

lwid=c(1,5)
lhei=c(1,5)

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

colnames(mat) <- col_data$sample
pheatmap(mat, annotation_col = column_data,
         margins=c(10,10))


res$gene_significance <- ifelse(res$padj < 0.05, "significant", "not significant")
data_frame <- as.data.frame(res)
qplot(log2FoldChange, -log(res$padj), data = data_frame, col = gene_significance)

table(res$gene_significance)

# further analysis on normalization using deseq2
# https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html#gene-ontology-enrichment-analysis

deseq2table <- estimateSizeFactors(deseqdata)
sizeFactors(deseq2table)



multidensity(counts(deseq2table, normalized = T)[c(1,10,26,44,56,72),],
              xlab="mean counts")

# dispersion estimates

dispersion <- estimateDispersions(deseq2table)
plotDispEsts(dispersion)

# results of differential expression using Benjamin–Hochberg correction for false discovery rate
# find how many genes are differentially expressed with the correction
wald_test_significance <-  nbinomWaldTest(dispersion)
gene_significance_with_bh <- results(wald_test_significance, pAdjustMethod = "BH")
table(gene_significance_with_bh$padj < 0.05)

plotMA(gene_significance_with_bh)
ordered_p_value <- gene_significance_with_bh[order(gene_significance_with_bh$pvalue),]
topgenes <- row.names(ordered_p_value[1:5,])
with(gene_significance_with_bh[topgenes, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=1, lwd=2)
  text(baseMean, log2FoldChange, topgenes, pos=2, col="red")
})

# annotation information for differentially expressed genes as measured with BH correction
sigGenes <- rownames(subset(gene_significance_with_bh, padj < 0.05))


anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=rownames(gene_significance_with_bh), 
                              columns=c("SYMBOL","SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")


anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))
nrow(anSig)
# see the annotation information for 25 of the differentially expressed genes
sample_n(anSig, 25)

