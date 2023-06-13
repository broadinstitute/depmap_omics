
# b012f426-0347-4433-8d16-f4298a134614 is the id for alt-aware bams
# mkdir -p clean
# touch clean/cclf.23Q2.alt.json
# gsutil cp -r clean gs://broad-qqin-sparkles/BamMetrics/
# gsutil rm gs://broad-qqin-sparkles/BamMetrics/clean/cclf.23Q2.alt.json

croo gs://broad-qqin-sparkles/BamMetrics/b012f426-0347-4433-8d16-f4298a134614/metadata.json --out-def-json gs://broad-qqin-sparkles/BamMetrics/clean/cclf.23Q2.alt.json --out-dir gs://broad-qqin-sparkles/BamMetrics/clean 

