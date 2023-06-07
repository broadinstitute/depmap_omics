#!/bin/bash -ex

python3 process_cosmic_data.py CosmicCodingMuts.hg38.normal.vcf.gz CosmicMutantExport_hg38.tsv.gz

~/bin/echtvar encode \
            cosmic_coding_ev.zip \
            cosmic.json \
            CosmicCodingMuts.hg38.normal_GRCh38_processed.vcf.gz

bcftools norm -m- ~/vep_test/CDS-9ImTq5.hg38-unfiltered.vcf -w 10000 -f ~/vep_test/Homo_sapiens_assembly38.fasta.gz -O b -o ~/vep_test/CDS-9ImTq5.hg38-unfiltered.bcf

~/bin/echtvar anno -e gnomad.v3.1.2.echtvar.v2.zip ~/vep_test/CDS-9ImTq5.hg38-unfiltered.bcf test_out.bcf

~/bin/echtvar anno -e ~/echtvar/gnomad.v3.1.2.echtvar.v2.zip CosmicCodingMuts.hg38.normal_GRCh38_processed.vcf.gz CosmicCodingMuts.hg38.normal_GRCh38_processed.bcf

