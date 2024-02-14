
init: list 
	# first generate a service account json key file
	# here we used sparkle spray to generate key 
	# https://github.com/ENCODE-DCC/caper/tree/master/scripts/gcp_caper_server
	bash create_instance.sh caper-runner broad-qqin broad-qqin.json gs://broad-qqin-sparkles -b 100

list:
	gcloud compute instances list --format 'value(name)'

stop:
	gcloud compute instances stop caper-runner --zone us-central1-a --project broad-qqin

start:
	gcloud compute instances start caper-runner --zone us-central1-a --project broad-qqin
	sleep 60
	gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="cd /opt/caper && screen -dmS caper_server bash -c 'caper server > caper_server.log 2>&1'"
	#  --backend-file gcp.conf
	sleep 10
			    
delete:
	gcloud compute instances delete caper-runner --zone us-central1-a --project broad-qqin
	gcloud compute disks delete caper-runner --zone us-central1-a --project broad-qqin

test_submit:
	gcloud beta compute scp ~/Workspace/depmap_omics/sandbox/cnv/hello.wdl qqin@caper-runner:hello.wdl --zone us-central1-a --project broad-qqin
	gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper submit hello.wdl"

status:
	gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper list"

ssh:
	gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin

compute-ssh:
	gcloud compute config-ssh	

submit:
	rsync -auviz tasks workflow.wdl input.json* input_cclf*.json input_cclf_all.json input_cclf_all_noalt.json caper-runner.us-central1-a.broad-qqin:
	#gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper submit workflow.wdl -i input.json"
	#gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper submit workflow.wdl -i input_cclf.json"
	#gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper submit workflow.wdl -i input_cclf_all.json"
	gcloud beta compute ssh --zone us-central1-a caper-runner --project broad-qqin --command="caper submit workflow.wdl -i input_cclf_all_noalt.json"

