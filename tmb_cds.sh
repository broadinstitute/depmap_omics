inputvcf=$1
sampleid=$2
cds_bed=$3
bedtools intersect -a ${inputvcf} -b ${cds_bed} -header > ${sampleid}.cds_only.vcf
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t[%AD]\t%AF\t%RS\t%gnomADe_AF\t%gnomADg_AF\n' ${sampleid}.cds_only.vcf | \
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
}' > ${sampleid}.cds_only.maf
cds_size=$(awk 'BEGIN{sum=0} {sum += ($3 - $2)} END {print sum}' "${cds_bed}")
n_mutations=$(wc -l < "${sampleid}.cds_only.maf")
tmb=$(awk -v n_mutations="$n_mutations" -v cds_size="$cds_size" 'BEGIN {tmb = 1e6 * n_mutations / cds_size; print tmb}')
echo "$tmb"

