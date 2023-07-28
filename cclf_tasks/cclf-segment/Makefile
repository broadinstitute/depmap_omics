

download:
	#gsutil cp gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324/submissions/00e5937e-8471-4f52-b43f-9665bfa40698/aggregate_CN_segments_wrkflw/db8a9c91-f4d7-4d21-a04e-10ff066ce5ed/call-aggregate_CN_segments/23Q2.called.seg .
	#gsutil cp gs://fc-secure-d2a2d895-a7af-4117-bdc7-652d7d268324/bdd767bd-f535-4716-88cb-0faeff5e5ec7/aggregate_CN_segments_wrkflw/a816ac7e-2fc6-4165-b98e-ec39975ab9fa/call-aggregate_CN_segments/attempt-2/all.called.seg .
	gsutil cp gs://fc-secure-bd7b8bc9-f665-4269-997e-5a402088a369/submissions/67a3964c-cb59-4321-bef0-8e4894ee59f7/aggregate_CN_segments_wrkflw/7ad1a42c-9e54-408e-b70a-a9e78ef6dcf3/call-aggregate_CN_segments/all.called.seg all.called.wgs.seg


mkdir: auth
	gsutil mb -l us-east1 gs://broad-snakemake-bucket
	# not working
	# gsutil acl ch -u "snakemake-test@broad-qqin.iam.gserviceaccount.com":W gs://broad-snakemake-bucket
	#602149243264-compute@developer.gserviceaccount.com
	#gsutil acl ch -u "602149243264-compute@developer.gserviceaccount.com":W gs://broad-snakemake-bucket
	#gsutil acl ch -u "602149243264-compute@developer.gserviceaccount.com":R gs://broad-snakemake-bucket
	# gcloud storage buckets add-iam-policy-binding  gs://broad-snakemake-bucket --member serviceAccount:602149243264-compute@developer.gserviceaccount.com --role=roles/storage.admin
	# This is a must-have step for snakemake
	# due to snakemake do not have --service-account option yet
	gcloud storage buckets add-iam-policy-binding  gs://broad-snakemake-bucket --member serviceAccount:602149243264-compute@developer.gserviceaccount.com --role=roles/storage.objectAdmin

list:
	gcloud iam service-accounts list

create:
	gcloud iam service-accounts create snakemake-test \
	    --description="snakemake runner" \
	    --display-name="snakemake"
	#grant permission here https://console.cloud.google.com/iam-admin/iam?project=broad-qqin
	for role in storage.objectCreator storage.objectViewer lifesciences.workflowsRunner lifesciences.admin iam.serviceAccountUser storage.admin; do \
	    gcloud projects add-iam-policy-binding broad-qqin --member serviceAccount:"snakemake-test@broad-qqin.iam.gserviceaccount.com" --role roles/$$role ; \
	done

auth:
	gcloud auth activate-service-account --key-file snakemake-credential.json

create-key: create
	gcloud iam service-accounts keys create snakemake-credential.json \
	  --iam-account=snakemake-test@broad-qqin.iam.gserviceaccount.com

run:
	# gcloud auth login
	# gcloud services enable lifesciences.googleapis.com
	# python upload_google_storage.py broad-snakemake-bucket/data/ data/
	snakemake --preemption-default 10 -j 3 --google-lifesciences --default-remote-prefix broad-snakemake-bucket --use-conda --google-lifesciences-region us-east1

local:
# not work on Mac
	snakemake -j 4

manual:
	# https://cloud.google.com/life-sciences/docs/process-genomic-data
	gsutil mb gs://broad-qqin-life-sciences/
	gsutil ls gs://broad-qqin-life-sciences/
        # need to be run in shell below
	gcloud beta lifesciences pipelines run \
	--labels test-lifesci \
	--command-line 'samtools index ${BAM} ${BAI}' \
	--docker-image "gcr.io/cloud-lifesciences/samtools" \
	--inputs BAM=gs://genomics-public-data/NA12878.chr20.sample.bam \
	--outputs BAI=gs://broad-qqin-life-sciences/NA12878.chr20.sample.bam.bai --service-account-email=snakemake-test@broad-qqin.iam.gserviceaccount.com

delete:
	gsutil rm gs://broad-qqin-life-sciences/NA12878.chr20.sample.bam.bai
	gsutil rb gs://broad-qqin-life-sciences

describe:
	gcloud beta lifesciences operations wait 4854895555446362058
	gcloud beta lifesciences operations describe projects/602149243264/locations/us-central1/operations/4854895555446362058
        gsutil ls gs://broad-qqin-life-sciences/NA12878.chr20.sample.bam.bai
