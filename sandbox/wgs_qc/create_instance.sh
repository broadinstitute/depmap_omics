#!/bin/bash
set -eo pipefail

if [[ $# -lt 1 ]]; then
  echo "Automated shell script to create Caper server instance with PostgreSQL on Google Cloud."
  echo
  echo "Usage: ./create_instance.sh [INSTANCE_NAME] [GCP_PRJ]"
  echo "                            [GCP_SERVICE_ACCOUNT_SECRET_JSON_FILE] [GCP_OUT_DIR]"
  echo "                            <OPTIONAL_ARGUMENTS>"
  echo
  echo "Positional arguments:"
  echo "  [INSTANCE_NAME]: New instance's name."
  echo "  [GCP_PRJ]: Your project's ID on Google Cloud Platform. --gcp-prj in Caper."
  echo "  [GCP_SERVICE_ACCOUNT_KEY_JSON_FILE]: Service account's secret key JSON file. --gcp-service-account-key-json in Caper."
  echo "  [GCP_OUT_DIR]: gs:// bucket dir path for outputs. --gcp-out-dir in Caper."
  echo
  echo "Optional arguments for Caper:"
  echo "  -l, --gcp-loc-dir: gs:// bucket dir path for localization."
  echo "  --gcp-region: Region for Google Life Sciences API. us-central1 by default. CHECK SUPPORTED REGIONS. This is different from --zone, which is used for instance creation only. us-central1 by default."
  echo "  --postgresql-db-ip: localhost by default."
  echo "  --postgresql-db-port: 5432 by default."
  echo "  --postgresql-db-user: cromwell by default."
  echo "  --postgresql-db-password: cromwell by default."
  echo "  --postgresql-db-name: cromwell by default."
  echo
  echo "Optional arguments for instance creation (gcloud compute instances create):"
  echo "  -z, --zone: Zone. Check available zones: gcloud compute zones list. us-central1-a by default."
  echo "  -m, --machine-type: Machine type. Check available machine-types: gcloud compute machine-types list. n1-standard-4 by default."
  echo "  -b, --boot-disk-size: Boot disk size. Use a suffix for unit. e.g. GB and MB. 100GB by default."
  echo "  -u, --username: Username (super user) used for transferring key file to the instance. ubuntu by default."
  echo "  --boot-disk-type: Boot disk type. pd-standard (Standard persistent disk) by default."
  echo "  --image: Image. Check available images: gcloud compute images list. ubuntu-2004-focal-v20220118 by default."
  echo "  --image-project: Image project. ubuntu-os-cloud by default."
  echo "  --tags: Tags to apply to the new instance. caper-server by default."
  echo "  --startup-script: Startup script CONTENTS (NOT A FILE). These command lines should sudo-install screen, Java, PostgreSQL, Python3 and pip3. DO NOT INSTALL CAPER HERE. some apt-get command lines by default."
  echo

  if [[ $# -lt 4 ]]; then
    echo "Define all positional arguments."
  fi
  exit 1
fi

# parse opt args first.
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -l|--gcp-loc-dir)
      GCP_LOC_DIR="$2"
      shift
      shift
      ;;
    --gcp-region)
      GCP_REGION="$2"
      shift
      shift
      ;;
    --postgresql-db-ip)
      POSTGRESQL_DB_IP="$2"
      shift
      shift
      ;;
    --postgresql-db-port)
      POSTGRESQL_DB_PORT="$2"
      shift
      shift
      ;;
    --postgresql-db-user)
      POSTGRESQL_DB_USER="$2"
      shift
      shift
      ;;
    --postgresql-db-password)
      POSTGRESQL_DB_PASSWORD="$2"
      shift
      shift
      ;;
    --postgresql-db-name)
      POSTGRESQL_DB_NAME="$2"
      shift
      shift
      ;;
    -z|--zone)
      ZONE="$2"
      shift
      shift
      ;;
    -m|--machine-type)
      MACHINE_TYPE="$2"
      shift
      shift
      ;;
    -b|--boot-disk-size)
      BOOT_DISK_SIZE="$2"
      shift
      shift
      ;;
    -u|--username)
      USERNAME="$2"
      shift
      shift
      ;;
    --boot-disk-type)
      BOOT_DISK_TYPE="$2"
      shift
      shift
      ;;
    --image)
      IMAGE="$2"
      shift
      shift
      ;;
    --image-project)
      IMAGE_PROJECT="$2"
      shift
      shift
      ;;
    --tags)
      TAGS="$2"
      shift
      shift
      ;;
    --startup-script)
      STARTUP_SCRIPT="$2"
      shift
      shift
      ;;
    -*)
      echo "Wrong parameter: $1."
      shift
      exit 1
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

# restore pos args.
set -- "${POSITIONAL[@]}"

# parse pos args.
INSTANCE_NAME="$1"
GCP_PRJ="$2"
GCP_SERVICE_ACCOUNT_KEY_JSON_FILE="${3/#\~/$HOME}"
GCP_OUT_DIR="$4"

# set defaults for opt args. (caper)
if [[ -z "$GCP_LOC_DIR" ]]; then
  GCP_LOC_DIR="$GCP_OUT_DIR"/.caper_tmp
fi
if [[ -z "$GCP_REGION" ]]; then
  GCP_REGION=us-central1
fi
if [[ -z "$POSTGRESQL_DB_IP" ]]; then
  POSTGRESQL_DB_IP=localhost
fi
if [[ -z "$POSTGRESQL_DB_PORT" ]]; then
  POSTGRESQL_DB_PORT=5432
fi
if [[ -z "$POSTGRESQL_DB_USER" ]]; then
  POSTGRESQL_DB_USER=cromwell
fi
if [[ -z "$POSTGRESQL_DB_PASSWORD" ]]; then
  POSTGRESQL_DB_PASSWORD=cromwell
fi
if [[ -z "$POSTGRESQL_DB_NAME" ]]; then
  POSTGRESQL_DB_NAME=cromwell
fi

# set defaults for opt args. (gcloud)
if [[ -z "$ZONE" ]]; then
  ZONE=us-central1-a
fi
if [[ -z "$MACHINE_TYPE" ]]; then
  MACHINE_TYPE=n1-standard-4
fi
if [[ -z "$BOOT_DISK_SIZE" ]]; then
  BOOT_DISK_SIZE=100GB
fi
if [[ -z "$USERNAME" ]]; then
  USERNAME=ubuntu
fi
if [[ -z "$BOOT_DISK_TYPE" ]]; then
  BOOT_DISK_TYPE=pd-standard
fi
if [[ -z "$IMAGE" ]]; then
  IMAGE=ubuntu-2004-focal-v20220118
fi
if [[ -z "$IMAGE_PROJECT" ]]; then
  IMAGE_PROJECT=ubuntu-os-cloud
fi
if [[ -z "$TAGS" ]]; then
  TAGS=caper-server
fi
if [[ -z "$STARTUP_SCRIPT" ]]; then
  STARTUP_SCRIPT="""
sudo apt-get update
sudo apt-get -y install screen python3 python3-pip default-jre postgresql postgresql-contrib acl
"""
fi

# validate all args.
if [[ -z "$GCP_PRJ" ]]; then
  echo "[GCP_PRJ] is not valid."
  exit 1
fi
if [[ "$GCP_OUT_DIR" != gs://* ]]; then
  echo "[GCP_OUT_DIR] should be a GCS bucket path starting with gs://"
  exit 1
fi
if [[ "$GCP_LOC_DIR" != gs://* ]]; then
  echo "-l, --gcp-loc-dir should be a GCS bucket path starting with gs://"
  exit 1
fi
if [[ ! -f "$GCP_SERVICE_ACCOUNT_KEY_JSON_FILE" ]]; then
  echo "[GCP_SERVICE_ACCOUNT_KEY_JSON_FILE] does not exists."
  exit 1
fi
if [[ "$POSTGRESQL_DB_IP" == localhost && "$POSTGRESQL_DB_PORT" != 5432 ]]; then
  echo "--postgresql-db-port should be 5432 for locally installed PostgreSQL (--postgresql-db-ip localhost)."
  exit 1
fi

# constants for files/params on instance.
GCP_AUTH_SH="/etc/profile.d/gcp-auth.sh"
CAPER_CONF_DIR=/opt/caper
ROOT_CAPER_CONF_DIR=/root/.caper
GLOBAL_CAPER_CONF_FILE="$CAPER_CONF_DIR/default.conf"
REMOTE_KEY_FILE="$CAPER_CONF_DIR/service_account_key.json"

# prepend more init commands to the startup-script
STARTUP_SCRIPT="""#!/bin/bash
### make caper's directories
sudo mkdir -p $CAPER_CONF_DIR
sudo mkdir -p $CAPER_CONF_DIR/local_loc_dir $CAPER_CONF_DIR/local_out_dir

### set default permission on caper's directories
sudo chmod 777 -R $CAPER_CONF_DIR
sudo setfacl -R -d -m u::rwX $CAPER_CONF_DIR
sudo setfacl -R -d -m g::rwX $CAPER_CONF_DIR
sudo setfacl -R -d -m o::rwX $CAPER_CONF_DIR

### make caper conf file
cat <<EOF > $GLOBAL_CAPER_CONF_FILE
# caper
backend=gcp
no-server-heartbeat=True
# cromwell
max-concurrent-workflows=300
max-concurrent-tasks=1000
# local backend
local-out-dir=$CAPER_CONF_DIR/local_out_dir
local-loc-dir=$CAPER_CONF_DIR/local_loc_dir
# gcp backend
gcp-prj=$GCP_PRJ
gcp-region=$GCP_REGION
gcp-out-dir=$GCP_OUT_DIR
gcp-loc-dir=$GCP_LOC_DIR
gcp-service-account-key-json=$REMOTE_KEY_FILE
use-google-cloud-life-sciences=True
# metadata DB
db=postgresql
postgresql-db-ip=$POSTGRESQL_DB_IP
postgresql-db-port=$POSTGRESQL_DB_PORT
postgresql-db-user=$POSTGRESQL_DB_USER
postgresql-db-password=$POSTGRESQL_DB_PASSWORD
postgresql-db-name=$POSTGRESQL_DB_NAME
EOF
sudo chmod +r $GLOBAL_CAPER_CONF_FILE

### soft-link conf file for root
sudo mkdir -p $ROOT_CAPER_CONF_DIR
sudo ln -s $GLOBAL_CAPER_CONF_FILE $ROOT_CAPER_CONF_DIR

### google auth shared for all users
sudo touch $GCP_AUTH_SH
echo \"gcloud auth activate-service-account --key-file=$REMOTE_KEY_FILE\" > $GCP_AUTH_SH
echo \"mkdir -p ~/.caper\" >> $GCP_AUTH_SH
echo \"ln -s /opt/caper/default.conf ~/.caper/ 2> /dev/null | true\" >> $GCP_AUTH_SH
echo \"export GOOGLE_APPLICATION_CREDENTIALS=$REMOTE_KEY_FILE\" >> $GCP_AUTH_SH

$STARTUP_SCRIPT
"""

# append more init commands to the startup-script
STARTUP_SCRIPT="""$STARTUP_SCRIPT
### init PostgreSQL for Cromwell
sudo -u postgres createuser root -s
sudo createdb $POSTGRESQL_DB_NAME
sudo psql -d $POSTGRESQL_DB_NAME -c \"create extension lo;\"
sudo psql -d $POSTGRESQL_DB_NAME -c \"create role $POSTGRESQL_DB_USER with superuser login password '$POSTGRESQL_DB_PASSWORD'\"

### upgrade pip and install caper croo
sudo python3 -m pip install --upgrade pip
sudo pip install --ignore-installed caper croo
"""

echo "$(date): Google auth with service account key file."
gcloud auth activate-service-account --key-file="$GCP_SERVICE_ACCOUNT_KEY_JSON_FILE"
export GOOGLE_APPLICATION_CREDENTIALS="$GCP_SERVICE_ACCOUNT_KEY_JSON_FILE"

echo "$(date): Making a temporary startup script..."
echo "$STARTUP_SCRIPT" > tmp_startup_script.sh

echo "$(date): Creating an instance..."
gcloud --project "$GCP_PRJ" compute instances create \
  "$INSTANCE_NAME" \
  --boot-disk-size="$BOOT_DISK_SIZE" \
  --boot-disk-type="$BOOT_DISK_TYPE" \
  --machine-type="$MACHINE_TYPE" \
  --zone="$ZONE" \
  --image="$IMAGE" \
  --image-project="$IMAGE_PROJECT" \
  --tags="$TAGS" \
  --metadata-from-file startup-script=tmp_startup_script.sh
echo "$(date): Created an instance successfully."

echo "$(date): Deleting the temporary startup script..."
rm -f tmp_startup_script.sh

while [[ $(gcloud --project "$GCP_PRJ" compute instances describe "${INSTANCE_NAME}" --zone "${ZONE}" --format="value(status)") -ne "RUNNING" ]]; do
    echo "$(date): Waiting for 20 seconds for the instance to spin up..."
    sleep 20
done

echo "$(date): If key file transfer fails for several times then manually transfer it to $REMOTE_KEY_FILE on the instance."
echo "$(date): Transferring service account key file to the instance..."
until gcloud --project "$GCP_PRJ" compute scp "$GCP_SERVICE_ACCOUNT_KEY_JSON_FILE" "$USERNAME"@"$INSTANCE_NAME":"$REMOTE_KEY_FILE" --zone="$ZONE"; do
  echo "$(date): Key file transfer failed. Retrying in 20 seconds..."
  sleep 20
done
echo "$(date): Transferred a key file to instance successfully."

echo "$(date): Waiting for the instance finishing up installing Caper..."
until gcloud --project "$GCP_PRJ" compute ssh --zone="$ZONE" "$USERNAME"@"$INSTANCE_NAME" --command="caper -v"; do
  echo "$(date): Caper has not been installed yet. Retrying in 40 seconds..."
  sleep 40
done
echo "$(date): Finished installing Caper on the instance. Ready to run Caper server."

echo "$(date): Spinning up Caper server..."
gcloud --project "$GCP_PRJ" compute ssh --zone="$ZONE" "$USERNAME"@"$INSTANCE_NAME" --command="cd $CAPER_CONF_DIR && sudo screen -dmS caper_server bash -c \"sudo caper server > caper_server.log 2>&1\""
sleep 60
until gcloud --project "$GCP_PRJ" compute ssh --zone="$ZONE" "$USERNAME"@"$INSTANCE_NAME" --command="caper list"; do
  echo "$(date): Caper server has not been started yet. Retrying in 60 seconds..."
  sleep 60
done
echo
echo "$(date): Caper server is up and ready to take submissions."
echo "$(date): You can find Caper server log file at $CAPER_CONF_DIR/caper_server.log."
echo "$(date): Cromwell's STDERR will be written to $CAPER_CONF_DIR/cromwell.out*."
echo
echo "$(date): Use the following command line to SSH to the instance."
echo
echo "gcloud beta compute ssh --zone $ZONE $INSTANCE_NAME --project $GCP_PRJ"
echo
