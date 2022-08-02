FROM python:3.9-slim-buster

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"

# Installing apt packages
RUN apt-get update && apt-get -y upgrade && apt-get install -y build-essential curl git libbz2-dev libcurl3-dev liblzma-dev libgsl-dev libncurses5-dev wget zip && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Installing msisensor2
RUN git clone https://github.com/niu-lab/msisensor2.git
RUN cd msisensor2
RUN chmod +x msisensor2

# Fixing path issues
ENV PATH="/msisensor2:${PATH}"
