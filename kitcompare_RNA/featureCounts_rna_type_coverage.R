library(Rsubread)
library(DESeq2)
library(biomaRt)
library(data.table)
library(EDASeq)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(stringi)
library(xlsx)
library(RColorBrewer)
library(scater)
library(Rsome)
library(gridExtra)
library(grid)


dir <- "/Users/mattsdwatson/VIB_proj_3425/merged_bam/"
setwd(dir)
bam_files <- list.files(pattern = ".bam$", recursive = TRUE,
                        include.dirs = FALSE)
bam_files
  
short_file_names <- c(list.files(dir))
complete_names <- short_file_names[-c(73,74)]
complete_names

# split the file names to retain just the unique identifier and not the file extension
sample_names <- as.factor(str_split_fixed(bam_files, "_m", 2)[,1])
kits <- c(rep("kit_a", 6), rep(c("kit_c", "kit_b"), 6))
experimental_set_up <- data.frame(sample = sample_names, kit = kits)
experimental_set_up <- experimental_set_up[order(experimental_set_up$kit),]
rownames(experimental_set_up) <- 1:nrow(experimental_set_up)

annotation_file <- "/Users/mattsdwatson/star/index/vib_annotations/exp2323-genes.gtf"
gene_path <- "/Users/mattsdwatson/VIB_proj_3425/mm10_gene_length_gc_content.txt"
gene_table <- read.table(gene_path, header = TRUE, fill = TRUE, sep = ",")
colnames(gene_table) <- c("gc", "start", "end", "gene")
gene_table$length <- gene_table$end - gene_table$start + 1
# kit a 
features_kit_a <- featureCounts(files = bam_files[1:6],
                               annot.ext = annotation_file,
                               isGTFAnnotationFile = TRUE,
                               useMetaFeatures = TRUE,
                               strandSpecific = 2)

feature_counts_kit_a <- as.data.frame(features_kit_a$counts)
colnames(feature_counts_kit_a) <- sample_names[1:6]

#re format the columns and rows so that the samples can be plotted with ggplot stacked barplots
features_kit_a_format <- cbind(rep(row.names(feature_counts_kit_a), 6), gather(feature_counts_kit_a))

#re-format the sample name in the re-formatted column
features_kit_a_format$key <- as.factor(str_split_fixed(features_kit_a_format$key, "_m", 2)[,1])
colnames(features_kit_a_format) <- c("gene", "sample", "counts")

annotations_file <- "/Users/mattsdwatson/star/index/vib_annotations/exp2323-genes.tab.txt"
annotations <- read.table(annotations_file)
annotations <- annotations[,c(1,10)]
annotations <- annotations[-1,]
colnames(annotations) <- c("gene", "rna_type")
annotations$rna_type <- as.factor(annotations$rna_type)
total_possible_annotations <- data.frame(annotations %>% group_by(rna_type) 
                                         %>% summarise(counts = n()))

feature_frame_kit_a <- merge(features_kit_a_format, annotations, by.x = "gene",
                       by.y = "gene")

grouped_by_feature_a <- as.data.frame(subset(feature_frame_kit_a, select = c(sample, counts, rna_type)) %>% group_by(sample, rna_type) %>%
                                      summarise(cum_counts = sum(counts)))

# look at the mirna genes to evaluate the length
# since we set the minimum read length after filtering to 35 bp, we
# expect to see precursor mirna transcripts that are 35 np or longer
# note that mature mirnas are processed to a maximum length of 25 bp,
# which would not be detetcable using our filtering length
only_mirna <- subset(feature_frame_kit_a, rna_type == "miRNA")
mirna_genes <- subset(gene_table, gene_table$gene %in% only_mirna$gene)
min(mirna_genes$length, na.rm = TRUE)
max(mirna_genes$length, na.rm = TRUE)

only_mirna <- merge(only_mirna, mirna_genes, by.x = "gene", by.y = "gene")
only_mirna <- only_mirna[,c(1,2,3,8)]
only_mirna <- only_mirna[order(-only_mirna$counts),]
grid.table(head(only_mirna, 25))

# merge all of the brewer pallettes to create a large, distinct ones for all of
# the genomic features for plotting
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot by both total counts as well as feature percentage
ggplot(grouped_by_feature_a, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("Cumulative RNA Type Coverage, Kit A, No Filtering/Normalization") +
  xlab("Sample") + ylab("Total Counts") + theme(legend.title = element_text(size = 7), 
                                         legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

ggplot(grouped_by_feature_a, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "fill", stat = "identity") + coord_flip() +
  ggtitle("RNA Type Percent Coverage, Kit A, No Filtering/Normalization") +
  xlab("Sample") + ylab("Percentage") + theme(legend.title = element_text(size = 7), 
                                                legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

# look at the relative expression for kit a

genes_for_a <- subset(gene_table, gene_table$gene %in% row.names(feature_counts_kit_a))
genes_for_kit_a <- genes_for_a[,c(1,5)]
row.names(genes_for_kit_a) <- genes_for_a$gene

feature_counts_kit_a <- feature_counts_kit_a[ order(row.names(feature_counts_kit_a)), ]
genes_for_kit_a <- genes_for_kit_a[ order(row.names(genes_for_kit_a)), ]

feature_counts_kit_a <- feature_counts_kit_a[row.names(feature_counts_kit_a) %in% genes_for_a$gene,]
row.names(genes_for_kit_a) <- row.names(feature_counts_kit_a)

expression_set_kita <- newSeqExpressionSet(counts = as.matrix(feature_counts_kit_a),
                                            phenoData = data.frame(kit = rep("kita", 6), 
                                                                   row.names = colnames(feature_counts_kit_a)),
                                            featureData = data.frame(genes_for_kit_a),
                                            row.names = colnames(feature_counts_kit_a))

dataWithin_kita <- withinLaneNormalization(expression_set_kita, "gc", which = "full")
dataNorm_kita <- betweenLaneNormalization(dataWithin_kita, which="full",
                                                round = TRUE)

# calculate the tpms for kit a to perform filtering based on different expression values

feature_counts_kit_a <- as.data.frame(assayData(dataNorm_kita)$normalizedCounts)
features_kit_a_format <- cbind(rep(row.names(feature_counts_kit_a), 6), gather(feature_counts_kit_a))
colnames(features_kit_a_format) <- c("gene", "sample", "counts")
gene_information_kit_a <- as.data.frame(genes_for_kit_a)
gene_information_kit_a <- data.frame(gene = row.names(genes_for_kit_a),
                                     gc = genes_for_kit_a$gc,
                                     length = genes_for_kit_a$length)
features_kit_a_format <- merge(features_kit_a_format, gene_information_kit_a, by.x = "gene",
                               by.y = "gene")
colnames(features_kit_a_format) <- c("gene", "sample", "counts", "gc", "length")

tpms_kit_a <- as.data.frame(features_kit_a_format %>% group_by(sample))
tpms_kit_a$counts <- as.numeric(tpms_kit_a$counts)
total_counts_kit_a <- as.data.frame(tpms_kit_a %>% group_by(sample) %>% summarise(total = sum(counts,
                                                                                na.rm = TRUE)))

tpms_kit_a <- merge(tpms_kit_a, total_counts_kit_a, by.x = "sample", by.y = "sample",
                    na.rm = TRUE)
nrow(tpms_kit_a)
tpms_kit_a$tpm <- tpms_kit_a$counts / (tpms_kit_a$total / 1000000)
tpms_kit_a <- tpms_kit_a[complete.cases(tpms_kit_a),]
nrow(tpms_kit_a)
tpms_kit_a <- tpms_kit_a[tpms_kit_a$tpm >= 1,]
nrow(tpms_kit_a)

tpms_kit_a <- tpms_kit_a %>% group_by(gene) %>% filter(n()==6)
nrow(tpms_kit_a)
tpms_kit_a <- as.data.frame(tpms_kit_a)
length(unique((tpms_kit_a$gene)))

to_spread <- subset(tpms_kit_a, select = c(gene, sample, counts))
spread <- spread(to_spread, key = sample, value = counts)
par(mar=c(8, 6, 3, 3))
EDASeq::plotRLE(as.matrix(spread[,-1]), las = 2, col = "blue", ylim = c(-4, 4),
                main = "RLE, Kit A, Normalized Counts",
                ylab = "Log2 fold change wrt median")
EDASeq::plotRLE(as.matrix(feature_counts_kit_a), las = 2, col = "blue", ylim = c(-5.5, 5.5),
                main = "RLE, Kit A, Raw/Non-Normalized Counts",
                ylab = "Log2 fold change wrt median")

# find the tpm for each gene to perform filtering

# kit b
# do not set strandedness for either kit b or c
kit_b_indices <- seq(8,18,2)
features_kit_b <- featureCounts(files = bam_files[kit_b_indices],
                                annot.ext = annotation_file,
                                isGTFAnnotationFile = TRUE,
                                useMetaFeatures = TRUE,
                                strandSpecific = 0)

feature_counts_kit_b <- as.data.frame(features_kit_b$counts)
colnames(feature_counts_kit_b) <- sample_names[kit_b_indices]

#re format the columns and rows so that the samples can be plotted with ggplot stacked barplots
features_kit_b_format <- cbind(rep(row.names(feature_counts_kit_b), 6), gather(feature_counts_kit_b))

#re-format the sample name in the re-formatted column
features_kit_b_format$key <- as.factor(str_split_fixed(features_kit_b_format$key, "_m", 2)[,1])
colnames(features_kit_b_format) <- c("gene", "sample", "counts")

feature_frame_kit_b <- merge(features_kit_b_format, annotations, by.x = "gene",
                             by.y = "gene")

grouped_by_feature_b <- as.data.frame(subset(feature_frame_kit_b, select = c(sample, counts, rna_type)) %>% group_by(sample, rna_type) %>%
                                      summarise(cum_counts = sum(counts)))

# merge all of the brewer palettes to create a large, distinct ones for all of
# the genomic features
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot by both total counts as well as feature percentage
ggplot(grouped_by_feature_b, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("Cumulative RNA Type Coverage, Kit B, No Filtering/Normalization") +
  xlab("Sample") + ylab("Total Counts") + theme(legend.title = element_text(size = 7), 
                                                legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

ggplot(grouped_by_feature_b, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "fill", stat = "identity") + coord_flip() +
  ggtitle("RNA Type Percent Coverage, Kit B, No Filtering/Normalization") +
  xlab("Sample") + ylab("Percentage") + theme(legend.title = element_text(size = 7), 
                                              legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

genes_for_b <- subset(gene_table, gene_table$gene %in% row.names(feature_counts_kit_b))
genes_for_kit_b <- genes_for_b[,c(1,5)]
row.names(genes_for_kit_b) <- genes_for_b$gene

feature_counts_kit_b <- feature_counts_kit_b[ order(row.names(feature_counts_kit_b)), ]
genes_for_kit_b <- genes_for_kit_b[ order(row.names(genes_for_kit_b)), ]

feature_counts_kit_b <- feature_counts_kit_b[row.names(feature_counts_kit_b) %in% genes_for_b$gene,]
row.names(genes_for_kit_b) <- row.names(feature_counts_kit_b)

expression_set_kitb <- newSeqExpressionSet(counts = as.matrix(feature_counts_kit_b),
                                           phenoData = data.frame(kit = rep("kitb", 6), 
                                                                  row.names = colnames(feature_counts_kit_b)),
                                           featureData = data.frame(genes_for_kit_b),
                                           row.names = colnames(feature_counts_kit_b))

# perform normalization procedure for kit a separate from the other kits
# between lane and within lane normalization for library size and gc content

dataWithin_kitb <- withinLaneNormalization(expression_set_kitb,"gc", which = "full")
dataNorm_kitb <- betweenLaneNormalization(dataWithin_kitb, which="full")

# calculate the tpms for kit a to perform filtering based on different expression values

feature_counts_kit_b <- as.data.frame(assayData(dataNorm_kitb)$normalizedCounts)
features_kit_b_format <- cbind(rep(row.names(feature_counts_kit_b), 6), gather(feature_counts_kit_b))
colnames(features_kit_b_format) <- c("gene", "sample", "counts")
gene_information_kit_b <- as.data.frame(genes_for_kit_b)
gene_information_kit_b <- data.frame(gene = row.names(gene_information_kit_b),
                                     gc = gene_information_kit_b$gc,
                                     length = gene_information_kit_b$length)
features_kit_b_format <- merge(features_kit_b_format, gene_information_kit_b, by.x = "gene",
                               by.y = "gene")
colnames(features_kit_b_format) <- c("gene", "sample", "counts", "gc", "length")

tpms_kit_b <- as.data.frame(features_kit_b_format %>% group_by(sample))
tpms_kit_b$counts <- as.numeric(tpms_kit_b$counts)
total_counts_kit_b <- as.data.frame(tpms_kit_b %>% group_by(sample) %>% summarise(total = sum(counts,
                                                                                              na.rm = TRUE)))

tpms_kit_b <- merge(tpms_kit_b, total_counts_kit_b, by.x = "sample", by.y = "sample",
                    na.rm = TRUE)
nrow(tpms_kit_b)
tpms_kit_b$tpm <- tpms_kit_b$counts / (tpms_kit_b$total / 1000000)
tpms_kit_b <- tpms_kit_b[complete.cases(tpms_kit_b),]
tpms_kit_b <- tpms_kit_b[tpms_kit_b$tpm >= 1,]
nrow(tpms_kit_b)
tpms_kit_b <- tpms_kit_b %>% group_by(gene) %>% filter(n()==6)
nrow(tpms_kit_b)
tpms_kit_b <- as.data.frame(tpms_kit_b)
length(unique((tpms_kit_b$gene)))

to_spread <- subset(tpms_kit_b, select = c(gene, sample, counts))
spread <- spread(to_spread, key = sample, value = counts)
par(mar=c(8, 6, 3, 3))
EDASeq::plotRLE(as.matrix(spread[,-1]), las = 2, col = "red", ylim = c(-3, 3),
                main = "RLE, Kit B, Normalized Counts",
                ylab = "Log2 fold change wrt median")
EDASeq::plotRLE(as.matrix(feature_counts_kit_b), las = 2, col = "red", ylim = c(-6, 6),
                main = "RLE, Kit B, Raw/Non-Normalized Counts",
                ylab = "Log2 fold change wrt median")

# kit c
# do not set strandedness for either kit b or c
kit_c_indices <- seq(7,17,2)
features_kit_c <- featureCounts(files = bam_files[kit_c_indices],
                                annot.ext = annotation_file,
                                isGTFAnnotationFile = TRUE,
                                useMetaFeatures = TRUE,
                                strandSpecific = 0)

feature_counts_kit_c <- as.data.frame(features_kit_c$counts)
colnames(feature_counts_kit_c) <- sample_names[kit_c_indices]

#re format the coumns and rows so that the samples can be plotted with ggplot stacked barplots
features_kit_c_format <- cbind(rep(row.names(feature_counts_kit_c), 6), gather(feature_counts_kit_c))

#re-format the sample name in the re-formatted column
features_kit_c_format$key <- as.factor(str_split_fixed(features_kit_c_format$key, "_m", 2)[,1])
colnames(features_kit_c_format) <- c("gene", "sample", "counts")

feature_frame_kit_c <- merge(features_kit_c_format, annotations, by.x = "gene",
                             by.y = "gene")

grouped_by_feature_c <- as.data.frame(subset(feature_frame_kit_c, select = c(sample, counts, rna_type)) %>% group_by(sample, rna_type) %>%
                                        summarise(cum_counts = sum(counts)))

# merge all of the brewer pallettes to create a large, distinct ones for all of
# the genomic features
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot by both total counts as well as feature percentage
ggplot(grouped_by_feature_c, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA Type Percent Coverage, Kit C, No Filtering/Normalization") +
  xlab("Sample") + ylab("Total Counts") + theme(legend.title = element_text(size = 7), 
                                                legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

ggplot(grouped_by_feature_c, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "fill", stat = "identity") + coord_flip() +
  ggtitle("RNA Type Percent Coverage, Kit C, No Filtering/Normalization") +
  xlab("Sample") + ylab("Percentage") + theme(legend.title = element_text(size = 7), 
                                              legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

genes_for_c <- subset(gene_table, gene_table$gene %in% row.names(feature_counts_kit_c))
genes_for_kit_c <- genes_for_c[,c(1,5)]
row.names(genes_for_kit_c) <- genes_for_c$gene

feature_counts_kit_c <- feature_counts_kit_c[ order(row.names(feature_counts_kit_c)), ]
genes_for_kit_c <- genes_for_kit_c[ order(row.names(genes_for_kit_c)), ]

feature_counts_kit_c <- feature_counts_kit_c[row.names(feature_counts_kit_c) %in% genes_for_c$gene,]
row.names(genes_for_kit_c) <- row.names(feature_counts_kit_c)

expression_set_kitc <- newSeqExpressionSet(counts = as.matrix(feature_counts_kit_c),
                                           phenoData = data.frame(kit = rep("kitc", 6), 
                                                                  row.names = colnames(feature_counts_kit_c)),
                                           featureData = data.frame(genes_for_kit_c),
                                           row.names = colnames(feature_counts_kit_c))

# perform normalization procedure for kit a separate from the other kits
# between lane and within lane normalization for library size and gc content

dataWithin_kitc <- withinLaneNormalization(expression_set_kitc,"gc", which = "full")
dataNorm_kitc <- betweenLaneNormalization(dataWithin_kitc, which="full")

# calculate the tpms for kit a to perform filtering based on different expression values

feature_counts_kit_c <- as.data.frame(assayData(dataNorm_kitc)$normalizedCounts)
features_kit_c_format <- cbind(rep(row.names(feature_counts_kit_c), 6), gather(feature_counts_kit_c))
colnames(features_kit_c_format) <- c("gene", "sample", "counts")
gene_information_kit_c <- as.data.frame(genes_for_kit_c)
gene_information_kit_c <- data.frame(gene = row.names(gene_information_kit_c),
                                     gc = gene_information_kit_c$gc,
                                     length = gene_information_kit_c$length)
features_kit_c_format <- merge(features_kit_c_format, gene_information_kit_c, by.x = "gene",
                               by.y = "gene")
colnames(features_kit_c_format) <- c("gene", "sample", "counts", "gc", "length")

tpms_kit_c <- as.data.frame(features_kit_c_format %>% group_by(sample))
tpms_kit_c$counts <- as.numeric(tpms_kit_c$counts)
total_counts_kit_c <- as.data.frame(tpms_kit_c %>% group_by(sample) %>% summarise(total = sum(counts,
                                                                                              na.rm = TRUE)))

tpms_kit_c <- merge(tpms_kit_c, total_counts_kit_c, by.x = "sample", by.y = "sample",
                    na.rm = TRUE)
nrow(tpms_kit_c)
tpms_kit_c$tpm <- tpms_kit_c$counts / (tpms_kit_c$total / 1000000)
tpms_kit_c <- tpms_kit_c[complete.cases(tpms_kit_c),]
tpms_kit_c <- tpms_kit_c[tpms_kit_c$tpm >= 1,]
nrow(tpms_kit_c)
tpms_kit_c <- tpms_kit_c %>% group_by(gene) %>% filter(n()==6)
nrow(tpms_kit_c)
tpms_kit_c <- as.data.frame(tpms_kit_c)
length(unique((tpms_kit_c$gene)))

to_spread <- subset(tpms_kit_c, select = c(gene, sample, counts))
spread <- spread(to_spread, key = sample, value = counts)
par(mar=c(8, 6, 3, 3))
EDASeq::plotRLE(as.matrix(spread[,-1]), las = 2, col = "green", ylim = c(-3, 3),
                main = "RLE, Kit C, Normalized Counts",
                ylab = "Log2 fold change wrt median")
EDASeq::plotRLE(as.matrix(feature_counts_kit_c), las = 2, col = "green", ylim = c(-6, 6),
                main = "RLE, Kit C, Raw/Non-Normalized Counts",
                ylab = "Log2 fold change wrt median")

# produce the total annotation coverage for each kit using the normalized counts

norm_counts_kit_c <- subset(tpms_kit_c, select = c(sample, gene, counts))

features_format_normalized <- merge(norm_counts_kit_c, annotations, by.x = "gene",
                             by.y = "gene")

grouped_by_feature_normalized <- as.data.frame(subset(features_format_normalized, select = c(sample, counts, rna_type)) %>% group_by(sample, rna_type) %>%
                                        summarise(cum_counts = sum(counts)))
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot by both total counts as well as feature percentage
ggplot(grouped_by_feature_normalized, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("Cumulative RNA Type Coverage, Kit C, Filtered/Normalized") +
  xlab("Sample") + ylab("Total Counts") + theme(legend.title = element_text(size = 7), 
                                                legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

ggplot(grouped_by_feature_normalized, aes(x = reorder(sample, cum_counts), y = cum_counts, fill = rna_type)) +
  geom_bar(position = "fill", stat = "identity") + coord_flip() +
  ggtitle("RNA Type Percent Coverage, Kit C, Filtered/Normalized") +
  xlab("Sample") + ylab("Percentage") + theme(legend.title = element_text(size = 7), 
                                              legend.text = element_text(size = 5.5)) +
  scale_fill_manual(values = col_vector)

# produce a bar chart of the different assignment categories generated by featureCounts
assignment_categories <- features_kit_a$stat$Status

kit_a_categories <- cbind(assignment_categories, features_kit_a$stat)
kita_average <- rowMeans(features_kit_a$stat[,2:7])
kita_average
kita_average_percentage <- kita_average/sum(kita_average)
kita_average_percentage

kit_b_categories <- cbind(assignment_categories, features_kit_b$stat)
kitb_average <- rowMeans(features_kit_b$stat[,2:7])
kitb_average
kitb_average_percentage <- kitb_average/sum(kitb_average)
kitb_average_percentage

kit_c_categories <- cbind(assignment_categories, features_kit_c$stat)
kitc_average <- rowMeans(features_kit_c$stat[,2:7])
kitc_average
kitc_average_percentage <- kitc_average/sum(kitc_average)
kitc_average_percentage

feature_summary <- cbind(assignment_categories, kita_average_percentage*100, kitb_average_percentage*100, kitc_average_percentage*100)
colnames(feature_summary) <- c("Genomic_Feature", "Kit_A", "Kit_B",
                               "Kit_C")
feature_summary

kita_to_plot <- cbind(assignment_categories, kita_average_percentage*100)
kitb_to_plot <- cbind(assignment_categories, kitb_average_percentage*100)
kitc_to_plot <- cbind(assignment_categories, kitc_average_percentage*100)

# use a ggplot to visualize the average of the genomic features
# do not include features for which there are no counts
features_to_plot <- as.data.frame(rbind(kita_to_plot, kitb_to_plot, kitc_to_plot))
features_to_plot
colnames(features_to_plot) <- c("feature", "percentage")
features_to_plot$kit <- c(rep("Kit A", 14), rep("Kit B", 14), rep("Kit C", 14))
features_to_plot <- features_to_plot[features_to_plot$percentage > 0,]
features_to_plot$percentage <- as.numeric(features_to_plot$percentage)
ggplot(data = features_to_plot, aes(fill=kit, y=percentage, x=feature)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=round(percentage, 2)), position=position_dodge(width=0.9), vjust=-0.25)

# plot the individual coverage profiles for each kit- if the feature is covered by 
# at least one of the samples in the kit, retain it

ind_cov_kit_a <- data.frame(gene = row.names(feature_counts_kit_a))
ind_cov_kit_a <- cbind(ind_cov_kit_a, feature_counts_kit_a[,1:6])
row.names(ind_cov_kit_a) <- c()
ind_cov_kit_a$present <- ifelse(rowSums(ind_cov_kit_a[,2:7]) >= 1, 1, 0)
kit_a_coverage <- merge(ind_cov_kit_a, annotations, by.x = "gene",
                        by.y = "gene")

# get a tally of how many of each feature is captured by kit a 
# from all of the possible features that are found in the gene tab txt file
kit_a_features_covered <- as.data.frame(subset(kit_a_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present)))
kit_a_features_covered <- merge(data.frame(kit_a_features_covered), total_possible_annotations,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_a_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_a_features_covered[,2:3] <- as.numeric(unlist(kit_a_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_possible_annotations$rna_type,
                              number_captured = kit_a_features_covered[,3] - kit_a_features_covered[,2])
coverage_plot_kit_a <- rbind(kit_a_features_covered[,1:2], number_captured)
coverage_plot_kit_a$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_a

ggplot(coverage_plot_kit_a, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit A") +
  xlab("RNA Type") + ylab("Count")

# individual feature coverage for kit b

ind_cov_kit_b <- data.frame(gene = row.names(feature_counts_kit_b))
ind_cov_kit_b <- cbind(ind_cov_kit_b, feature_counts_kit_b[,1:6])
row.names(ind_cov_kit_b) <- c()
ind_cov_kit_b$present <- ifelse(rowSums(ind_cov_kit_b[,2:7]) >= 1, 1, 0)
kit_b_coverage <- merge(ind_cov_kit_b, annotations, by.x = "gene",
                        by.y = "gene")

# get a tally of how many of each feature is captured by kit a 
# from all of the possible features that are found in the gene tab txt file
kit_b_features_covered <- as.data.frame(subset(kit_b_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present)))
kit_b_features_covered <- merge(data.frame(kit_b_features_covered), total_possible_annotations,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_b_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_b_features_covered[,2:3] <- as.numeric(unlist(kit_b_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_possible_annotations$rna_type,
                              number_captured = kit_b_features_covered[,3] - kit_b_features_covered[,2])
coverage_plot_kit_b <- rbind(kit_b_features_covered[,1:2], number_captured)
coverage_plot_kit_b$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_b

ggplot(coverage_plot_kit_b, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit B") +
  xlab("RNA Type") + ylab("Count")

# individual feature coverage for kit c

ind_cov_kit_c <- data.frame(gene = row.names(feature_counts_kit_c))
ind_cov_kit_c <- cbind(ind_cov_kit_c, feature_counts_kit_c[,1:6])
row.names(ind_cov_kit_c) <- c()
ind_cov_kit_c$present <- ifelse(rowSums(ind_cov_kit_c[,2:7]) >= 1, 1, 0)
kit_c_coverage <- merge(ind_cov_kit_c, annotations, by.x = "gene",
                        by.y = "gene")

# get a tally of how many of each feature is captured by kit a 
# from all of the possible features that are found in the gene tab txt file
kit_c_features_covered <- as.data.frame(subset(kit_c_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present)))
kit_c_features_covered <- merge(data.frame(kit_c_features_covered), total_possible_annotations,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_c_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_c_features_covered[,2:3] <- as.numeric(unlist(kit_c_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_possible_annotations$rna_type,
                              number_captured = kit_c_features_covered[,3] - kit_c_features_covered[,2])
coverage_plot_kit_c <- rbind(kit_c_features_covered[,1:2], number_captured)
coverage_plot_kit_c$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_c

ggplot(coverage_plot_kit_c, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit C") +
  xlab("RNA Type") + ylab("Count")

