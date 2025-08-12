#!/bin/bash

inputvcf=$1
sampleid=$2
hg38_mapp=$3
cds_bed=$4
bam_file=$5
hg38_genome_ref=$6
read_depth_threshold=$7
base_qual_threshold=$8

### Denominator = number of bases genomewide covered by at least read_depth_threshold and median base quality > base_qual_threshold

total_covered_bases_for_denominator=`samtools mpileup -q 60 -f ${hg38_genome_ref}  ${bam_file} | \
gawk -vread_depth_threshold=${read_depth_threshold} -vqual_threshold=${base_qual_threshold}  'BEGIN { for(i=0;i<256;i++) {ord[sprintf("%c",i)] = i} total_bases_covered=0;}
{  if ($4 >= read_depth_threshold) {
    count = 0
    for (i=1; i<=length($6); i++) {
        if (ord[substr($6,i,1)] - 33 >= qual_threshold)
            count++
    }
    if (count >= 0.5*$4) {printf("%s\t%d\t%d\n",$1,$2,$2) >> "mapped_bases.bed"; total_bases_covered++}
}} END{print total_bases_covered}'`

echo ",Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2" > header.csv
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%AD]\t%AF\t%RS\t%gnomADe_AF\t%gnomADg_AF\n' ${inputvcf} | \
awk -vsampleid=${sampleid} '{ 
	ad_split=split($6,arr,","); 
	total_reads = arr[1]+arr[2]; 
	alt_reads = arr[2]; 
	gsub("chr","",$1)
	ref=$3
	alt=$4
	start=$2
	if (length($3) > length($4)) { 
		ind=index($3,$4)
		ref=substr($3,ind+1)
		start=start+ind
		alt="-"
	} else if (length($3) < length($4)) { 
		ind=index($4,$3)
		ref="-"
		start=start+ind
		alt=substr($4,ind+1)
	}
	if (total_reads >= 10 && alt_reads >= 5 && ($9 == "." || $9 < 0.000001) && ($10 == "." || $10 < 0.000001)) {
		printf("%s\t%d\t%d\t%s\t%s\t%s\n",$1,start,start,ref,alt,sampleid)
	}
}' > ${sampleid}.somatic_mutations_called.bed 

awk '{printf("%s,Unknown,%s,%s,%s,%s,%s\n",NR,$6,$1,$2,$4,$5)}' ${sampleid}.somatic_mutations_called.bed > ${sampleid}.somatic_maf_noheader.csv

awk -F"," '{print $3"\t"$4"\t"$4}' ${sampleid}.somatic_maf_noheader.csv > ${sampleid}.somatic_mutations_called.bed
bedtools intersect ${sampleid}.somatic_mutations_called.bed -b mappable_bases.bed -u > ${sampleid}.somatic_mutations_called_in_mappable_regions.bed
n_mutations_in_mappable_regions=`cat ${sampleid}.somatic_mutations_called_in_mappable_regions.bed | wc -l`
tmb=$(awk -v n_mutations="$n_mutations_in_mappable_regions" -v mappable_genome_size="$total_covered_bases_for_denominator" \
#    'BEGIN {tmb = 1e6 * n_mutations / mappable_genome_size; print tmb}')

# Prep format to call signatureanalyzer - comma separated maf format
cat header.csv ${sampleid}.somatic_maf_noheader.csv > ${sampleid}.somatic_maf.csv

# TMB CDS only
bedtools intersect -a ${sampleid}.somatic_mutations_called_in_mappable_regions.bed -b ${cds_bed} -u > ${sampleid}.somatic_mutations_in_mappable_cds.bed
n_mutations_cds_mappable_regions=`cat ${sampleid}.somatic_mutations_in_mappable_cds.bed | wc -l`
bedtools intersect -a mappable_bases -b ${cds_bed} -u > mappable_bases_in_cds.bed
cds_size=`cat mappable_bases_in_cds.bed | wc -l`
tmb_cds=$(awk -v n_mutations_cds="$n_mutations_cds" -v cds_size="$cds_size" 'BEGIN {tmb = 1e6 * n_mutations_cds / cds_size; print tmb}')

echo ${tmb} > tmb.out
echo ${tmb_cds} > tmb_cds.out
echo ${n_mutations_in_mappable_regions} > n_mutations.out
echo ${n_mutations_cds_mappable_regions} > n_mutations_cds.out
