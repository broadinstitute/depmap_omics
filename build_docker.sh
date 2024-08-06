docker build . -t tcga_depmapomics:latest --no-cache
docker tag tcga_depmapomics:latest us-docker.pkg.dev/depmap-omics/public/tcga_vcf_to_depmap:latest
docker push us-docker.pkg.dev/depmap-omics/public/tcga_vcf_to_depmap:latest