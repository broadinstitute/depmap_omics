#!/bin/bash 

bcftools norm -m- ~/vep_test/CDS-9ImTq5.hg38-unfiltered.vcf -w 10000 -f ~/vep_test/Homo_sapiens_assembly38.fasta.gz -O b -o ~/vep_test/CDS-9ImTq5.hg38-unfiltered.bcf

~/bin/echtvar anno -e gnomad.v3.1.2.echtvar.v2.zip ~/vep_test/CDS-9ImTq5.hg38-unfiltered.bcf test_out.bcf

