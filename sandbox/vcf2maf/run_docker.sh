#!/bin/bash -ex

chmod 777 /home/ubuntu/depmap_omics/sandbox/vcf2maf/

docker run -v $(pwd):$(pwd) vcf2maf perl vcf2maf/vcf2maf.pl --input-vcf $(pwd)/test.vcf --output-maf $(pwd)/test.maf --ref $(pwd)/Homo_sapiens_assembly38.fasta.gz --vep-path /opt/conda/envs/vep/bin/ --vep-data $(pwd)/vep_data/ --ncbi-build GRCh38 --tmp-dir $(pwd)/

#/opt/conda/envs/vep/bin/
#/opt/conda/envs/vep/bin/ 
