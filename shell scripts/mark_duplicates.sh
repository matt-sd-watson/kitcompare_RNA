bam_dir=/Users/mattsdwatson/VIB_proj_3425/merged_bam
for file in $(ls $bam_dir)
do
    SAMPLE=`basename $file`
    java -jar /Users/mattsdwatson/picard/picard-2.23.4/picard.jar MarkDuplicates \
      I=$bam_dir/${SAMPLE} \
      O=/Users/mattsdwatson/VIB_proj_3425/duplicates_marked/${SAMPLE} \
      M=/Users/mattsdwatson/VIB_proj_3425/duplicates_marked/marked_dup_metrics.txt \
      USE_JDK_DEFLATER=true 

done