working_dir=/Users/mattsdwatson/VIB_proj_3425/star_align
final_dir=/Users/mattsdwatson/VIB_proj_3425/merged_bam

ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE41_S32_L001_R1_001/ARE41_S32_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_S32_L002_R1_001/ARE41_S32_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_S32_L003_R1_001/ARE41_S32_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_S32_L004_R1_001/ARE41_S32_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE41_S32_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR56a_HT_S36_L001_R1_001/AR56a_HT_S36_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_HT_S36_L002_R1_001/AR56a_HT_S36_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_HT_S36_L003_R1_001/AR56a_HT_S36_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_HT_S36_L004_R1_001/AR56a_HT_S36_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR56a_HT_S36_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE44_HT_S41_L001_R1_001/ARE44_HT_S41_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_HT_S41_L002_R1_001/ARE44_HT_S41_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_HT_S41_L003_R1_001/ARE44_HT_S41_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_HT_S41_L004_R1_001/ARE44_HT_S41_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE44_HT_S41_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/31-AR68a_S3_L001_R1_001/31-AR68a_S3_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/31-AR68a_S3_L002_R1_001/31-AR68a_S3_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/31-AR68a_S3_L003_R1_001/31-AR68a_S3_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/31-AR68a_S3_L004_R1_001/31-AR68a_S3_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/31-AR68a_S3_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR68a_S22_L001_R1_001/AR68a_S22_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_S22_L002_R1_001/AR68a_S22_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_S22_L003_R1_001/AR68a_S22_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_S22_L004_R1_001/AR68a_S22_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR68a_S22_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE41_HT_S40_L001_R1_001/ARE41_HT_S40_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_HT_S40_L002_R1_001/ARE41_HT_S40_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_HT_S40_L003_R1_001/ARE41_HT_S40_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE41_HT_S40_L004_R1_001/ARE41_HT_S40_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE41_HT_S40_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/41-ARE41_S5_L001_R1_001/41-ARE41_S5_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/41-ARE41_S5_L002_R1_001/41-ARE41_S5_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/41-ARE41_S5_L003_R1_001/41-ARE41_S5_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/41-ARE41_S5_L004_R1_001/41-ARE41_S5_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/41-ARE41_S5_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE36_HT_S39_L001_R1_001/ARE36_HT_S39_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_HT_S39_L002_R1_001/ARE36_HT_S39_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_HT_S39_L003_R1_001/ARE36_HT_S39_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_HT_S39_L004_R1_001/ARE36_HT_S39_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE36_HT_S39_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR61a_HT_S37_L001_R1_001/AR61a_HT_S37_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_HT_S37_L002_R1_001/AR61a_HT_S37_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_HT_S37_L003_R1_001/AR61a_HT_S37_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_HT_S37_L004_R1_001/AR61a_HT_S37_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR61a_HT_S37_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/44-ARE44_S6_L001_R1_001/44-ARE44_S6_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/44-ARE44_S6_L002_R1_001/44-ARE44_S6_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/44-ARE44_S6_L003_R1_001/44-ARE44_S6_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/44-ARE44_S6_L004_R1_001/44-ARE44_S6_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/44-ARE44_S6_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR61a_S19_L001_R1_001/AR61a_S19_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_S19_L002_R1_001/AR61a_S19_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_S19_L003_R1_001/AR61a_S19_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR61a_S19_L004_R1_001/AR61a_S19_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR61a_S19_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/36-ARE36_S4_L001_R1_001/36-ARE36_S4_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/36-ARE36_S4_L002_R1_001/36-ARE36_S4_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/36-ARE36_S4_L003_R1_001/36-ARE36_S4_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/36-ARE36_S4_L004_R1_001/36-ARE36_S4_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/36-ARE36_S4_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR56a_S17_L001_R1_001/AR56a_S17_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_S17_L002_R1_001/AR56a_S17_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_S17_L003_R1_001/AR56a_S17_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR56a_S17_L004_R1_001/AR56a_S17_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR56a_S17_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/28-AR61a_S2_L001_R1_001/28-AR61a_S2_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/28-AR61a_S2_L002_R1_001/28-AR61a_S2_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/28-AR61a_S2_L003_R1_001/28-AR61a_S2_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/28-AR61a_S2_L004_R1_001/28-AR61a_S2_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/28-AR61a_S2_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE44_S35_L001_R1_001/ARE44_S35_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_S35_L002_R1_001/ARE44_S35_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_S35_L003_R1_001/ARE44_S35_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE44_S35_L004_R1_001/ARE44_S35_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE44_S35_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/ARE36_S27_L001_R1_001/ARE36_S27_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_S27_L002_R1_001/ARE36_S27_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_S27_L003_R1_001/ARE36_S27_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/ARE36_S27_L004_R1_001/ARE36_S27_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/ARE36_S27_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/26-AR56a_S1_L001_R1_001/26-AR56a_S1_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/26-AR56a_S1_L002_R1_001/26-AR56a_S1_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/26-AR56a_S1_L003_R1_001/26-AR56a_S1_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/26-AR56a_S1_L004_R1_001/26-AR56a_S1_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/26-AR56a_S1_merged.bam
ulimit -H -c unlimited
java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MergeSamFiles \
I=$working_dir/AR68a_HT_S38_L001_R1_001/AR68a_HT_S38_L001_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_HT_S38_L002_R1_001/AR68a_HT_S38_L002_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_HT_S38_L003_R1_001/AR68a_HT_S38_L003_R1_001Aligned.sortedByCoord.out.bam \
I=$working_dir/AR68a_HT_S38_L004_R1_001/AR68a_HT_S38_L004_R1_001Aligned.sortedByCoord.out.bam \
O=$final_dir/AR68a_HT_S38_merged.bam
