#!/bin/bash -ex

gsutil cp gs://fc-secure-bd7b8bc9-f665-4269-997e-5a402088a369/submissions/8b0e39a7-6875-4241-8913-78369b4b6173/aggregate_CN_segments_wrkflw/27a8b507-c60d-41c6-92f4-674e34b20d57/call-aggregate_CN_segments/all.called.seg all.called.seg.wgs

awk '$27>0.98' hg38_superdups | sed 1d | cut -f 2-4 > hg38_superdups98 

#for seg in latest_segments/23Q2.called.seg 23Q2.called.seg all.called.seg all.called.seg.wgs; do
#echo $seg
#awk '{print $2"\t"$3"\t"$4"\t"$1}' $seg  | sed 1d | grep -v "4.1e+07" | grep -v "1.2e+07" | grep -v "e+07" | bedtools intersect -wa -u -a - -b CNV_and_centromere_blacklist.hg38liftover.bed | cut -f 4 | sort -k 1 | uniq -c
#done

#for seg in latest_segments/23Q2.called.seg 23Q2.called.seg all.called.seg all.called.seg.wgs; do
#echo $seg
#awk '{print $2"\t"$3"\t"$4"\t"$1}' $seg  | sed 1d | grep -v "4.1e+07" | grep -v "1.2e+07" | grep -v "e+07" | bedtools intersect -wa -u -a - -b hg38_superdups98 | cut -f 4 | sort -k 1 | uniq -c
#done
