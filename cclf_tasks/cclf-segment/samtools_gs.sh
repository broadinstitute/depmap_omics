
export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

for bam in $(cut -f 2 -d, CCLF.bam_list | sed 1d); do
echo $(basename bam) $(gsutil -u broad-qqin ls $bam)
#echo $(basename bam) $(samtools view -c -q 20 $bam)
done
