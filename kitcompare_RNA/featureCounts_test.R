library(Rsubread)
library(DESeq2)
library(data.table)
library(EDASeq)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(stringi)
library(xlsx)

dir <- "/Users/mattsdwatson/VIB_proj_3425/star_align/"
setwd(dir)
bam_files <- list.files(pattern = ".bam$", recursive = TRUE,
                        include.dirs = FALSE)
bam_files
annotation_file <- "/Users/mattsdwatson/star/index/vib_annotations/exp2323-genes.gtf"

# kit a 
features_kit_a <- featureCounts(files = bam_files[1:24],
                                annot.ext = annotation_file,
                                isGTFAnnotationFile = TRUE,
                                useMetaFeatures = TRUE,
                                strandSpecific = 2,
                                ignoreDup = TRUE)

features_kit_a <- featureCounts(files = bam_files[25:72],
                                annot.ext = annotation_file,
                                isGTFAnnotationFile = TRUE,
                                useMetaFeatures = TRUE,
                                strandSpecific = 0)

features_test <- featureCounts(files = file.choose(),
                                annot.ext = annotation_file,
                                isGTFAnnotationFile = TRUE,
                                useMetaFeatures = TRUE,
                                strandSpecific = 0,
                                ignoreDup = TRUE)


library(yeastRNASeq)
data(geneLevelData)

mat <- as.matrix(geneLevelData)

data <- newSeqExpressionSet(mat,
                            phenoData=AnnotatedDataFrame(
                              data.frame(conditions=factor(c("mut", "mut", "wt", "wt")),
                                         row.names=colnames(geneLevelData))))


plotRLE(mat, col = c(2, 2, 3, 3))

example_sce <- mockSCE()
example_sce <- logNormCounts(example_sce)

plotRLE(example_sce, colour_by = "Mutation_Status", style = "minimal")

plotRLE(example_sce, colour_by = "Mutation_Status", style = "full",
        outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)


test_frame <- data.frame(person = c("john", "ann", "john", "ann", "ann"),
                         salary = c(3, 4, 5, 6, 7))

test_frame %>% group_by(person) %>% summarise(total = sum(salary))

