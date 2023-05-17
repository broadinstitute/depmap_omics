DOCKER_IMG=us.gcr.io/depmap-omics/hermit-dev-env:v1
docker build . --platform linux/amd64 -t ${DOCKER_IMG} && \
 docker push ${DOCKER_IMG}
 

