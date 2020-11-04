RUN_PATH=$1
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    mkdir /Users/mattsdwatson/VIB_proj_3425/fastqc_reports/${SAMPLE}
    fastqc -t 5 ${SAMPLE} -o /Users/mattsdwatson/VIB_proj_3425/fastqc_reports/${SAMPLE}
done

