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
library(scater)
library(Rqc)

dir <- "/Users/mattsdwatson/VIB_proj_3425/fastq/"
setwd(dir)
fastq_files <- list.files(pattern = ".fastq.gz$", recursive = TRUE,
                        include.dirs = FALSE)
fastq_files

dir <- "/Users/mattsdwatson/VIB_proj_3425/merged_bam/"
setwd(dir)
bam_files <- list.files(pattern = ".bam$", recursive = TRUE,
                        include.dirs = FALSE)
bam_files
bfs <- BamFileList(bam_files)

all_kit <- c(rep("kita", 6), rep(c("kitc", "kitb"),6))
colour <- ifelse(all_kit == "kita", "red", ifelse(all_kit == "kitb", "green", "blue"))
plotQuality(bfs, col = colour, lty = 1)
legend(locator(1), unique(all_kit), fill = unique(colour))

col_data <- data.frame(sample = sample_names, kit = all_kit)
col_data
elementMetadata(fastq) <- col_data
fastq

ShortRead::plotReadQuality(bfs)

# random sample from kit a, plot nucleotide distribution
plotNtFrequency(bfs[[3]], main = "Nucleotide Frequency By Cycle, 31-AR68a_S3_merged,
                Kit A")
plotNtFrequency(bfs[[8]], main = "Nucleotide Frequency By Cycle, AR56a_S17_merged,
                Kit B")
plotNtFrequency(bfs[[9]], main = "Nucleotide Frequency By Cycle, AR61a_HT_S37_merged,
                Kit C")


meanVarPlot(count_expression_set, log=TRUE)
biasPlot(count_expression_set, "gc")

MDPlot(count_expression_set)

fData(count_expression_set) <- fData(count_expression_set)[complete.cases(fData(count_expression_set)), ]

pattern = c(rep("red", 4), rep("blue", 4), rep("green", 4))
boxplot(count_expression_set, col = pattern, las = 2, cex = 0.7)

dataWithin <- withinLaneNormalization(count_expression_set,"gc")
dataNorm <- betweenLaneNormalization(dataWithin, which="full")
boxplot(dataNorm)

biasPlot(dataNorm, "gc", log=TRUE)

counts(dataNorm)

deseq_object <- DESeqDataSetFromMatrix(countData = counts(dataNorm),
                                       colData = pData(dataNorm),
                                       design = ~ kit)
dds <- DESeq(deseq_object)
res <- results(dds)
res
resultsNames(dds)
vsd <- vst(dds)

head(assay(vsd), 10)

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

poisd <- PoissonDistance(t(counts(deseq_object)))
pois_matrix <- as.matrix(poisd$dd)
rownames(pois_matrix) <- col_data$sample
pheatmap(pois_matrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

pca_plot <- plotPCA(vsd, intgroup = "kit")
pca_plot


ordered_p_value <- res[order(res$pvalue),]
summary(ordered_p_value)


plotMA(res, ylim=c(min(res$log2FoldChange),max(res$log2FoldChange)), alpha = 0.05)


topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 25)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

column_data <- data.frame(kit = col_data$kit)
row.names(column_data) <- col_data$sample
column_data

colnames(mat) <- col_data$sample
pheatmap(mat, annotation_col = column_data,
         margins=c(10,10))

res$gene_significance <- ifelse(res$padj < 0.05, "significant", "not significant")
table(res$gene_significance)


# plot the percentage of mapped reads for each kit
alignment_stats <- read.table(file.choose(), header = FALSE, fill = TRUE,
                              sep = "\t")

sample_indices_alignment <- seq(1, 214, 3)
percentage_alignment_indices <- seq(2, 215, 3)
samples_alignment <- data.frame(samples = alignment_stats$V1[sample_indices_alignment],
                                percentage = alignment_stats$V1[percentage_alignment_indices])

kitb_kitc <- c(rep("Kit C",4),rep("Kit B",4))
kitb_kitc <- rep(kitb_kitc, 6)
samples_alignment$kit <- c(rep("Kit A",24),kitb_kitc)
samples_alignment <- samples_alignment[,-1]
samples_alignment$percentage <- as.numeric(samples_alignment$percentage)

ggplot(samples_alignment, aes(x = kit, y = percentage)) + geom_boxplot(fill = c("red",
                                                                                "green",
                                                                                "blue")) +
  ggtitle("Alignment Percentage with STAR, All Kits") +
  xlab("Kit") + ylab("Percentage of Mapped reads (%)")
