
for job in $(gcloud beta lifesciences operations list --filter "DONE!=True" | sed 1d | cut -f 1 -d" "); do
    gcloud beta lifesciences operations cancel $job
done
