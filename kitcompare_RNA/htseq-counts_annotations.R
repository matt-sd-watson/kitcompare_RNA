library(dplyr)
library(stringi)
library(tidyverse)
library(ggplot2)

count_file_dir <- "/Users/mattsdwatson/VIB_proj_3425/htseq-count"
count_files <- list.files(count_file_dir, full.names = TRUE)
count_files                 

count_data <- list()
raw_counts <- list()
for (i in 1:18) {
  sample_name_stripped <- strsplit(count_files[i], "(/)|(__)")[[1]][6]
  htseq_counts <- read.table(count_files[i], fill = T)
  string_location <- as.data.frame(stri_locate_first(htseq_counts$V1, regex = "ENSMUSG"))
  start <- which(string_location$start == 1)[1]
  htseq_counts <- htseq_counts[start[1]:nrow(htseq_counts),1:2]
  colnames(htseq_counts) <- c("gene", "counts")
  rownames(htseq_counts) <- 1:nrow(htseq_counts)
  annotations_file <- "/Users/mattsdwatson/star/index/vib_annotations/exp2323-genes.tab.txt"
  annotations <- read.table(annotations_file)
  annotations <- annotations[,c(1,10)]
  annotations <- annotations[-1,]
  colnames(annotations) <- c("gene", "rna_type")
  annotations$rna_type <- as.factor(annotations$rna_type)
  
  # get the total number of possible annotations per features for gene counting
  total_features <- as.data.frame(annotations %>% group_by(rna_type) %>% summarise(total_features = n()))
  
  merged_counts <- merge(htseq_counts, annotations, by.x = "gene", by.y = "gene")
  merged_counts$counts <- as.numeric(merged_counts$counts)
  
  raw_counts[[i]] <- merged_counts[,1:2]
  
  grouped_by_feature <- as.data.frame(subset(merged_counts, select = c(counts, rna_type)) %>% group_by(rna_type) %>%
                                        summarise(cum_counts = sum(counts)))
  
  sample_counts <- data.frame(sample = grouped_by_feature$cum_counts)
  colnames(sample_counts) <- sample_name_stripped
  count_data[[i]] <- sample_counts
}

master_list <- cbind(grouped_by_feature$rna_type, do.call("cbind", count_data))
master_list
par(mfrow = c(1, 3),     
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(2, 9.5, 2, 2), # space for one row of text at ticks and to separate plots
    mgp = c(0, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA)    # allow content to protrude into outer margin (and beyond)

for (i in c(2:4)) {
  data <- master_list[,i]
  title_name <- colnames(master_list)[i]
  barplot(data, horiz = TRUE, names.arg = master_list[,1], las = 2, cex.axis = 1,
          cex.names = 0.8, col = "blue", main = title_name)
}

# get the coverage of individual features- if the features is present in any of the samples
# from one kit, it is counted

kit_a_coverage <- raw_counts[c(1:6)] %>% reduce(right_join, by = "gene")
kit_a_coverage$present <- ifelse(rowSums(kit_a_coverage[,2:7]) >= 1, 1, 0)
kit_a_coverage <- merge(kit_a_coverage, annotations, by.x = "gene",
                        by.y = "gene")

# get a tally of how many of each feature is captured by kit a 
# from all of the possible features that are found in the gene tab txt file
kit_a_features_covered <- subset(kit_a_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present))
kit_a_features_covered <- merge(data.frame(kit_a_features_covered), total_features,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_a_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_a_features_covered[,2:3] <- as.numeric(unlist(kit_a_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_features$rna_type,
                              number_captured = kit_a_features_covered[,3] - kit_a_features_covered[,2])
coverage_plot_kit_a <- rbind(kit_a_features_covered[,1:2], number_captured)
coverage_plot_kit_a$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_a

ggplot(coverage_plot_kit_a, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit A") +
  xlab("RNA Type") + ylab("Count")

#repeat for kit b and c
kit_b_indices <- seq(8, 18, 2)
kit_b_coverage <- raw_counts[c(kit_b_indices)] %>% reduce(right_join, by = "gene")
kit_b_coverage$present <- ifelse(rowSums(kit_b_coverage[,2:7]) >= 1, 1, 0)
kit_b_coverage <- merge(kit_b_coverage, annotations, by.x = "gene",
                        by.y = "gene")


kit_b_features_covered <- subset(kit_b_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present))
kit_b_features_covered <- merge(data.frame(kit_b_features_covered), total_features,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_b_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_b_features_covered[,2:3] <- as.numeric(unlist(kit_b_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_features$rna_type,
                              number_captured = kit_b_features_covered[,3] - kit_b_features_covered[,2])
coverage_plot_kit_b <- rbind(kit_b_features_covered[,1:2], number_captured)
coverage_plot_kit_b$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_b

ggplot(coverage_plot_kit_b, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit B") +
  xlab("RNA Type") + ylab("Count")


kit_c_indices <- seq(7, 17, 2)
kit_c_coverage <- raw_counts[c(kit_c_indices)] %>% reduce(right_join, by = "gene")
kit_c_coverage$present <- ifelse(rowSums(kit_c_coverage[,2:7]) >= 1, 1, 0)
kit_c_coverage <- merge(kit_c_coverage, annotations, by.x = "gene",
                        by.y = "gene")

kit_c_features_covered <- subset(kit_c_coverage, select = c(present, rna_type)) %>%
  group_by(rna_type) %>% summarise(tally = sum(present))
kit_c_features_covered <- merge(data.frame(kit_c_features_covered), total_features,
                                by.x = "rna_type", by.y = "rna_type")

colnames(kit_c_features_covered) <- c("rna_type", "number_captured", "total_features")
kit_c_features_covered[,2:3] <- as.numeric(unlist(kit_c_features_covered[,2:3]))
number_captured <- data.frame(rna_type = total_features$rna_type,
                              number_captured = kit_c_features_covered[,3] - kit_c_features_covered[,2])
coverage_plot_kit_c <- rbind(kit_c_features_covered[,1:2], number_captured)
coverage_plot_kit_c$coverage <- c(rep("covered", 43), rep("not covered", 43))
coverage_plot_kit_c

ggplot(coverage_plot_kit_c, aes(x = reorder(rna_type, number_captured), y = number_captured, fill = coverage)) +
  geom_bar(position = "stack", stat = "identity") + coord_flip() +
  ggtitle("RNA type feature coverage, Kit C") +
  xlab("RNA Type") + ylab("Count")



