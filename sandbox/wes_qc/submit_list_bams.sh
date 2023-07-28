#!/bin/bash 

set -o errexit
set -o nounset

readonly MY_PROJECT=${1}
readonly MY_BUCKET_PATH=${2}

readonly CONTAINER_PROJECT="$(echo "${MY_PROJECT}" | sed 's_:_/_')"
readonly OUTPUT_ROOT="${MY_BUCKET_PATH}/DNA_QC"
readonly SCRIPT_DIR="$(dirname "${0}")"

# Build the docker image
gcloud builds submit "${SCRIPT_DIR}" \
  --tag="gcr.io/${CONTAINER_PROJECT}/miniwdl"

# cut -f 1 ../all.called.seg | sed 1d | sort -k 1 -u
# dsub \
#   --provider google-cls-v2 \
#   --project "${MY_PROJECT}" \
#   --zones "us-central1-*" \
#   --logging "${OUTPUT_ROOT}/logging/" \
#   --disk-size 100 \
#   --name "miniwdl" \
#   --image "gcr.io/${CONTAINER_PROJECT}/miniwdl" \
#   --tasks "${SCRIPT_DIR}/submit_list.tsv" \
#   --command 'fastqc ${INPUT_BAM} --outdir=$(dirname ${OUTPUT_FILES})' \
#   --wait
