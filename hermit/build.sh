DOCKER_IMG=us-central1-docker.pkg.dev/cds-docker-containers/docker/depmap_omics-hermit-env:v2
docker build . --platform linux/amd64 -t ${DOCKER_IMG} && \
 docker push ${DOCKER_IMG}
 

