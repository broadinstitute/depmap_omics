inputvcf=$1
sampleid=$2
hg38_mapp=$3
cds_bed=$4
awk '{ if ($4 > 0.99) print }' ${hg38_mapp} > ${hg38_mapp}.high_map.bed
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
	if (alt_reads >= 5 && ($9 == "." || $9 < 0.000001) && ($10 == "." || $10 < 0.000001)) {
		printf("%s,Unknown,%s,%s,%s,%s,%s\n",NR,sampleid,$1,start,ref,alt)
	}
}' > ${sampleid}.somatic_maf_noheader.csv

# Prep format to call signatureanalyzer - comma separated maf format
cat header.csv ${sampleid}.somatic_maf_noheader.csv > ${sampleid}.somatic_maf.csv

# TMB whole genome
mappable_genome_size=$(awk 'BEGIN{sum_bases=0} {sum_bases += ($3 - $2)} END{print sum_bases}' "${hg38_mapp}.high_map.bed")
n_mutations=$(wc -l < "${sampleid}.somatic_maf_noheader.csv")
tmb=$(awk -v n_mutations="$n_mutations" -v mappable_genome_size="$mappable_genome_size" \
    'BEGIN {tmb = 1e6 * n_mutations / mappable_genome_size; print tmb}')

# TMB CDS only
awk -F"," '{ print "chr"$4"\t"$5"\t"$5}' ${sampleid}.somatic_maf_noheader.csv > ${sampleid}.somatic_maf.bed
n_mutations_cds=`bedtools intersect -a ${sampleid}.somatic_maf.bed -b ${cds_bed} -u | wc -l`
cds_size=$(awk 'BEGIN{sum=0} {sum += ($3 - $2)} END {print sum}' "${cds_bed}")
tmb_cds=$(awk -v n_mutations="$n_mutations" -v cds_size="$cds_size" 'BEGIN {tmb = 1e6 * n_mutations / cds_size; print tmb}')

echo ${tmb} > tmb.out
echo ${tmb_cds} > tmb_cds.out
echo ${n_mutations} > n_mutations.out
echo ${n_mutations_cds} > n_mutations_cds.out
