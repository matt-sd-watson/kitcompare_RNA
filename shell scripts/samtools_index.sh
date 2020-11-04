RUN_PATH=/Users/mattsdwatson/VIB_proj_3425/merged_bam/
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    samtools index ${SAMPLE}
done