## Build container
#FROM quay.io/condaforge/mambaforge:latest as conda
#
#RUN mamba info --envs && \
#    mamba install -y miniwdl
#
#RUN apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
#
#RUN miniwdl --version
#
#RUN docker run hello-world
#
#LABEL maintainer="depmap-omics team"
#
FROM ubuntu:20.04 as deps
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
RUN apt-get -qq update && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y \
    python3-pip python3-setuptools tzdata wget zip git-core default-jre jq shellcheck docker.io
RUN pip3 install -U pip  # due to infamous pyca/cryptography#5771

# add and become 'wdler' user -- it's useful to run the test suite as some arbitrary uid, because
# the runner has numerous file permissions-related constraints
RUN useradd -ms /bin/bash -u 1337 wdler
USER wdler
WORKDIR /home/wdler

ENV PYTHONPATH $PYTHONPATH:/home/wdler/.local/lib/python3.6
ENV PATH $PATH:/home/wdler/.local/bin

RUN pip install miniwdl
RUN miniwdl run_self_test
