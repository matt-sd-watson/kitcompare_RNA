RUN_PATH=$1
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    FILE="$(cut -d'.' -f1 <<<"$SAMPLE")"	
    STAR --genomeDir /Users/mattsdwatson/star/index/vib_reference/ --readFilesIn /Users/mattsdwatson/VIB_proj_3425/fastq_trim/${SAMPLE} \
    --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /Users/mattsdwatson/VIB_proj_3425/star_align/${FILE}/${FILE} \
    --alignIntronMin 50 --alignIntronMax 500000 \
    --sjdbGTFfile /Users/mattsdwatson/star/index/vib_annotations/exp2323-genes.gtf \
    --outSAMprimaryFlag OneBestScore --twopassMode Basic

done
