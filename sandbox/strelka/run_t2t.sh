#!/bin/bash -ex

#../minimap2/minimap2 -x sr -d chm13v2.0.mmi chm13v2.0.fa.gz

#samtools sort -n -@8 CDS-9ImTq5.hg38.bam | \
#samtools fastq - -1 CDS-9ImTq5_1.fq -2 CDS-9ImTq5_2.fq > /dev/null

#samtools sort -n -@8 AGILENT_10237.bam | \
#samtools fastq - -1 AGILENT_10237_1.fq -2 AGILENT_10237_2.fq > /dev/null

#../minimap2/minimap2 -t 8 -a -x sr chm13v2.0.mmi AGILENT_10237_1.fq  AGILENT_10237_2.fq | \
#	samtools fixmate -u -m - - | \
#	samtools sort -u -@2 -T /tmp/example_prefix - | \
#	samtools markdup -@8 --reference chm13v2.fa - /output/AGILENT_10237_chm13.bam

#../minimap2/minimap2 -t 8 -a -x sr chm13v2.0.mmi CDS-9ImTq5_1.fq -2 CDS-9ImTq5_2.fq | \
#	samtools fixmate -u -m - - | \
#	samtools sort -u -@2 -T /tmp/example_prefix2 - | \
#	samtools markdup -@8 --reference chm13v2.fa - /output/CDS-9ImTq5_chm13.bam


#samtools index /output/CDS-9ImTq5_chm13.bam
#samtools index /output/AGILENT_10237_chm13.bam

#strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
#	--normalBam /output/AGILENT_10237_chm13.bam \
#	--tumorBam /output/CDS-9ImTq5_chm13.bam \
#	--referenceFasta chm13v2.0.fa \
#	--runDir CDS-9ImTq5_vs_AGILENT_10237_chm

#wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain

#docker run -v $(pwd):$(pwd) quay.io/biocontainers/picard:2.26.10--hdfd78af_0 picard CreateSequenceDictionary R=$(pwd)/chm13v2.0.fa O=$(pwd)/chm13v2.0.fa.dict


#Failed
#docker run -v $(pwd):$(pwd) quay.io/biocontainers/picard:2.26.10--hdfd78af_0 picard -Xmx8072M LiftoverVcf I=$(pwd)/CDS-9ImTq5_vs_AGILENT_10237_chm/results/variants/somatic.snvs.vcf.gz O=$(pwd)/lifted_over.vcf CHAIN=$(pwd)/chm13v2-grch38.chain REJECT=$(pwd)/rejected_variants.vcf R=$(pwd)/chm13v2.0.fa
#all failed
#docker run -v $(pwd):$(pwd) quay.io/biocontainers/picard:2.26.10--hdfd78af_0 picard -Xmx8072M LiftoverVcf I=$(pwd)/CDS-9ImTq5_vs_AGILENT_10237_chm/results/variants/somatic.snvs.vcf.gz O=$(pwd)/lifted_over.vcf CHAIN=$(pwd)/chm13v2-hg38.over.chain REJECT=$(pwd)/rejected_variants.vcf R=$(pwd)/chm13v2.0.fa
#liftOver oldFile map.chain newFile unMapped
#./liftOver somatic.snvs.vcf $(pwd)/chm13v2-hg38.over.chain test.vcf test_un.vcf

#Direct search
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/liftover/chm13v2.0_ClinVar20220313.vcf.gz
