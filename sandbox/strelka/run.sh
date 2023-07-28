#!/bin/bash -ex

wget -c https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2

strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam AGILENT_10237.bam \
	--tumorBam CDS-9ImTq5.hg38.bam \
	--referenceFasta Homo_sapiens_assembly38.fasta \
	--runDir CDS-9ImTq5_vs_AGILENT_10237

docker run -v $(pwd):$(pwd) docker.io/ensemblorg/ensembl-vep vep -i $(pwd)/test.vcf -o $(pwd)/test_ann.vcf --species human --database

samtools tview -p chr21:43104346 CDS-9ImTq5.hg38.bam
samtools tview -p chr21:43104346 AGILENT_10237.bam
