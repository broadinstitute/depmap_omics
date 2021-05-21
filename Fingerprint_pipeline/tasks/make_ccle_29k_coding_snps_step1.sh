#!/bin/bash -l
# Download and install https://github.com/naumanjaved/fingerprint_maps
# Install GATK and bcftools
# Download example RNA-seq bam to select expressed SNPs and save as example_rna.bam.
# I used gs://cclebams/rnasq_hg38/CDS-QcwFpK.Aligned.sortedByCoord.out.bam as an example.
# Then with python 2.7 enviroment run this bash script for all chromosomes except for Y
CHR="22"
echo Running chromosome ${CHR} SNPs
echo Downloading chromosome ${CHR} SNPs
cd vcfs_1k_genomes
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
cd ..
echo Indexing chromosome ${CHR} SNPs
bcftools index vcfs_1k_genomes/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
echo Filtering chromosome ${CHR} SNPs by MAF
bcftools view vcfs_1k_genomes/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -i "AF>.1&AF<.9&VT='SNP'" --max-alleles 2 -o filtered_chr${CHR}.vcf -O v
echo Getting RNA-seq coverage for ${CHR} SNPs
java -Xmx6g -jar ~/Library/gatk/gatk-package-4.2.0.0-local.jar AnnotateVcfWithBamDepth -I example_rna.bam -V filtered_chr${CHR}.vcf -O rna_coverage_filtered_chr${CHR}.vcf
echo Filtering chromosome ${CHR} SNPs by RNA-seq coverage
bcftools view rna_coverage_filtered_chr${CHR}.vcf -i "BAM_DEPTH>10" -o coding_filtered_chr${CHR}.vcf -O v
echo Running build_fingerprint_maps for chromosome ${CHR} SNPs
python ~/Data/Fingerprinting/fingerprint_maps/build_fingerprint_maps.py \
 --recomb_file genetic_map_b37/genetic_map_chr${CHR}_combined_b37.txt --chrom ${CHR} \
 --VCF_file coding_filtered_chr${CHR}.vcf  --LD_script ldsc/ldsc.py --int_dir intermediates \
 --out_dir output_75_20 --prune_cutoff 0.20 --clump_cutoff 0.75
