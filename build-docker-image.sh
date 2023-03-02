#!/bin/bash

if [ "$1" = "" ]; then
    echo "Missing name to tag image with"
    exit 1
fi
IMAGE_TAG="$1"

set -ex

# Build Docker image
docker build . \
  -t ${IMAGE_TAG} \

