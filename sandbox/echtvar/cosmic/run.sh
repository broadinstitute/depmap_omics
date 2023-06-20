#!/bin/bash -ex

# only need to run once (now saved in google bucket):
# python3 process_cosmic_data.py CosmicCodingMuts.hg38.normal.vcf.gz CosmicMutantExport_hg38.tsv.gz

# ~/bin/echtvar encode \
#             cosmic_coding_ev.zip \
#             cosmic.json \
#             CosmicCodingMuts.hg38.normal_GRCh38_processed.vcf.gz

bcftools norm -m- /home/ubuntu/depmap_omics/CDS-9ImTq5.filtered.vcf.gz -w 10000 -f /home/ubuntu/depmap_omics/Homo_sapiens_assembly38.fasta -O b -o ~/CDS-9ImTq5.hg38-unfiltered.bcf

~/bin/echtvar anno -e gnomad.v3.1.2.echtvar.v2.zip ~/CDS-9ImTq5.hg38-unfiltered.bcf test_out.bcf

~/bin/echtvar anno -e ~/echtvar/gnomad.v3.1.2.echtvar.v2.zip CosmicCodingMuts.hg38.normal_GRCh38_processed.vcf.gz CosmicCodingMuts.hg38.normal_GRCh38_processed.bcf

# 1. make a test subset of data
# 2. make a docker file
# 3. build a docker image (docker build . -t depmapomics:test)
# 4. make a wdl file that optionally uses this docker image
# 5. write a test in the style of branch pgm-sample-test
# 6. poetry install
# 7. set python env in vscode
# 8. run pytest
