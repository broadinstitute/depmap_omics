#!/bin/bash -ex

for bam in $(cut -f 1 ../all.called.seg | sed 1d | sort -k 1 -u); do
	export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
	echo $bam $(samtools view -q 30 -c gs://cclebams/hg38_wes/${bam}.hg38.bam) >> all.called.seg.tsv
done
