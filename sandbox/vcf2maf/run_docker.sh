#!/bin/bash -ex

chmod 777 /home/ubuntu/depmap_omics/sandbox/vcf2maf/

#/home/ubuntu/mambaforge/envs/vep/bin/perl '/home/ubuntu/mambaforge/envs/vep/bin//vep' --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick --pick_order canonical,tsl,biotype,rank,ccds,length --dir '/home/ubuntu/depmap_omics/sandbox/vcf2maf/vep_data/' --fasta '/home/ubuntu/depmap_omics/sandbox/vcf2maf/Homo_sapiens_assembly38.fasta.gz' --format vcf --input_file '/home/ubuntu/depmap_omics/sandbox/vcf2maf/test.vcf' --output_file '/home/ubuntu/depmap_omics/sandbox/vcf2maf//test.vep.vcf' --offline --pubmed --fork 4 --polyphen b --af --af_1kg --af_gnomad --regulatory --force_overwrite
#
#/home/ubuntu/mambaforge/envs/vep/bin/perl '/home/ubuntu/mambaforge/envs/vep/bin//vep' --species homo_sapiens --cache --assembly GRCh38 --no_progress --no_stats --everything --dir /home/ubuntu/depmap_omics/sandbox/vcf2maf/vep_data/ --input_file test.vcf --output_file test2_vep.vcf --force_overwrite --offline --fasta Homo_sapiens_assembly38.fasta.gz --fork 8 --vcf --pick 

#for i in {1..20}; do
#    #echo $i
#    #perl vcf2maf/vcf2maf.pl --input-vcf $(pwd)/test2_vep.vcf --output-maf $(pwd)/test${i}.maf --ref $(pwd)/Homo_sapiens_assembly38.fasta.gz --vep-path /home/ubuntu/mambaforge/envs/vep/bin/ --vep-data $(pwd)/vep_data/ --ncbi-build GRCh38 --tmp-dir $(pwd)/ --inhibit-vep
#    diff $(pwd)/test1.maf $(pwd)/test${i}.maf
#done


for i in {1..20}; do
	#perl vcf2maf/vcf2maf.pl --input-vcf $(pwd)/test.vcf  --output-maf $(pwd)/test${i}.maf --ref $(pwd)/Homo_sapiens_assembly38.fasta.gz --vep-path /home/ubuntu/mambaforge/envs/vep/bin/ --vep-data $(pwd)/vep_data/ --ncbi-build GRCh38 --tmp-dir $(pwd)/
	docker run -v $(pwd):$(pwd) us.gcr.io/cds-docker-containers/vcf2maf:test perl vcf2maf/vcf2maf.pl --input-vcf $(pwd)/test.vcf  --output-maf $(pwd)/test${i}.maf --ref $(pwd)/Homo_sapiens_assembly38.fasta.gz --vep-path /opt/conda/envs/vep/bin/ --vep-data $(pwd)/vep_data/ --ncbi-build GRCh38 --tmp-dir $(pwd)/
	#/home/ubuntu/mambaforge/envs/vep/bin/perl '/home/ubuntu/mambaforge/envs/vep/bin//vep' --species homo_sapiens --cache --assembly GRCh38 --no_progress --no_stats --everything --dir /home/ubuntu/depmap_omics/sandbox/vcf2maf/vep_data/ --input_file test.vcf --output_file test.vep${i}.vcf --force_overwrite --offline --fasta Homo_sapiens_assembly38.fasta.gz --fork 8 --vcf --pick --mane
	mv test.vep.vcf test.vep${i}.vcf
done

#/opt/conda/envs/vep/bin/
#/opt/conda/envs/vep/bin/ 
