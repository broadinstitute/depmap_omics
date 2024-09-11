FROM ubuntu:24.04

# ----- start of hermit-specific requirements ------
#
# The requirements for a docker image to work with hermitcrab:
# 1. There must be a user named "ubuntu" with uid 2000
# 2. They must be able to sudo without entering a password
# 3. sshd must be installed 

RUN apt-get update && \
  apt-get install -y openssh-server sudo
RUN mkdir /run/sshd && \
  deluser ubuntu && \
  adduser ubuntu --disabled-password --uid 2000 --gecos "" && \
  echo 'ubuntu ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers 

RUN mkdir /home/ubuntu/.ssh && chmod go-rx /home/ubuntu/.ssh && chown ubuntu:ubuntu /home/ubuntu/.ssh

# ----- end of hermit-specific requirements ------

# install docker (just the command line tool, not the service) from Docker's package repo
RUN sudo apt-get install -y apt-transport-https ca-certificates curl software-properties-common && \
   curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
   echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null && \
   apt-get update && \
   apt-get install -y docker-ce-cli

# set this to avoid building image hanging due to packages prompting for configuration info
ENV DEBIAN_FRONTEND=noninteractive

# Install python 3.9 from PPA repo. Probably unnecessary, but thought it'd be good to be explict about version of python
# being installed. Can update the below to switch to whichever version of python one wants
RUN add-apt-repository -y ppa:deadsnakes/ppa && \
  apt-get update -y && apt-get install -y vim nano joe python3.9 python3.9-dev python3.9-venv

# install miniwdl in a python 3.9 virtual env
RUN mkdir -p /opt/miniwdl && \
    python3.9 -m venv /opt/miniwdl && \
    /opt/miniwdl/bin/python -m pip install miniwdl
    
RUN ln -s /opt/miniwdl/bin/miniwdl /usr/local/bin/miniwdl

# also install misc other useful utilities
# git because we'll need to checkout the code
# npm because pyright uses it to install some of its dependencies
# screen because it's generally useful
RUN apt-get install -y git npm screen

# install poetry. This is perhaps a strangely round about way, but it seems that the python3.9 install
# above doesn't include pip, but if we create a virtualenv, we get pip in that environment
RUN python3.9 -m venv /opt/poetry && /opt/poetry/bin/pip install poetry==1.8.2
RUN ln -s /opt/poetry/bin/poetry /usr/local/bin/poetry

# install google sdk per the debian instructions
RUN apt-get install apt-transport-https ca-certificates gnupg
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-cli -y
