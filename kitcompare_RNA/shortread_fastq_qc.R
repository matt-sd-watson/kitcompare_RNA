library(ShortRead)
library(Rqc)
library(seqplots)
fastq_dir = "/Users/mattsdwatson/VIB_proj_3425/fastq/"
setwd(fastq_dir)
fastq_files <- list.files(pattern = ".fastq.gz$", recursive = TRUE,
                          include.dirs = TRUE)
fastq_files

qa_summary <- qa(fastq_files)
quality_frame <- qa_summary[["readQualityScore"]]
plot(density(quality_frame$quality))

read_fastq <- readFastq(fastq_files)
