starting_dir=/Users/mattsdwatson/VIB_proj_3425/fastq_pp
final_dir=/Users/mattsdwatson/VIB_proj_3425/fastq_trim
for file in $(ls $starting_dir)
do
    SAMPLE=`basename $file`
    cutadapt -a 'TTTTTTTTTTTTTTT' -a 'AAAAAAAAAAAAAAA' -n 20 -m 35 -O 10 \
    $starting_dir/${SAMPLE} \
    -o $final_dir/${SAMPLE}
done