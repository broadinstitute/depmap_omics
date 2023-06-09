#!/bin/bash -ex

#CDS-bT5QTI DEPMAP
#CDS-rPCQiq CCLF
#CDS-Y5dmwC CCLF

#dsub \
#  --provider google-cls-v2 \
#  --project broad-qqin \
#  --regions us-central1 \
#  --logging gs://broad-snakemake-bucket/logging/ \
#  --output OUT=gs://broad-snakemake-bucket/output/out_demo.txt \
#  --input BAM=gs://cclebams/hg38_wes/CDS-00rz9N.hg38.bam \
#  --image quay.io/cancercollaboratory/dockstore-tool-samtools-index \
#  --command 'samtools flagstat ${BAM} 2>&1> ${OUT}' \
#  --user-project broad-qqin \
#  --wait

#gcloud config list && gcloud config set compute/region ""
#export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` samtools view -c $bam


#gcloud iam service-accounts create "dsub-broad-cds"
#gcloud iam service-accounts list
#gcloud projects add-iam-policy-binding broad-qqin --member serviceAccount:dsub-broad-cds@broad-qqin.iam.gserviceaccount.com --role roles/serviceusage.serviceUsageConsumer
#gsutil iam ch serviceAccount:dsub-broad-cds@broad-qqin.iam.gserviceaccount.com:roles/storage.objectAdmin gs://broad-snakemake-bucket
#gsutil iam ch serviceAccount:dsub-broad-cds@broad-qqin.iam.gserviceaccount.com:roles/storage.objectAdmin gs://artifacts.broad-qqin.appspot.com
#gsutil iam ch serviceAccount:602149243264-compute@developer.gserviceaccount.com:roles/storage.objectAdmin gs://artifacts.broad-qqin.appspot.com
#gsutil iam ch serviceAccount:sparkles-qy61ynrp5c@broad-qqin.iam.gserviceaccount.com:roles/storage.objectAdmin gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324
#gsutil iam ch serviceAccount:602149243264-compute@developer.gserviceaccount.com:roles/storage.objectAdmin gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324

#activate terra workspace access https://support.terra.bio/hc/en-us/articles/7448594459931
#cd firecloud-tools/
#./run.sh scripts/register_service_account/register_service_account.py  -j ../../../sandbox/cnv/broad-qqin.json -e qqin@broadinstitute.org
#cd ..

dsub \
  --provider google-cls-v2 \
  --project broad-qqin \
  --regions us-central1 \
  --logging gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324/logging \
  --output OUT=gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324/test/out_sa.txt \
  --input BAM=gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324/submissions/d2dfc36f-b916-44c6-9241-7c202ce8f82c/WGS_preprocessing/c37707b7-3bfa-4fc8-843e-8870abe43880/call-PreProcessingForVariantDiscovery_GATK4/PreProcessingForVariantDiscovery_GATK4/650cfe5a-b95f-4026-856a-d0b79eab2345/call-GatherBamFiles/CDS-1aJRXx.hg38.bam \
  --image quay.io/cancercollaboratory/dockstore-tool-samtools-index \
  --service-account sparkles-qy61ynrp5c@broad-qqin.iam.gserviceaccount.com \
  --command 'samtools view -c -q 20 ${BAM} 2>&1> ${OUT}' \
  --user-project broad-qqin \
  --wait

  ##--service-account 602149243264-compute@developer.gserviceaccount.com \
  
##gcloud projects add-iam-policy-binding depmap-omics --member serviceAccount:602149243264-compute@developer.gserviceaccount.com --role roles/storage.admin
##gcloud projects add-iam-policy-binding broad-qqin --member=serviceAccount:snakemake-test@broad-qqin.iam.gserviceaccount.com --role=roles/compute.networkUser
