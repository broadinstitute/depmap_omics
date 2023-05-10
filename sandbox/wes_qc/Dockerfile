# copy from https://raw.githubusercontent.com/broadinstitute/hermitcrab/main/dockerimage/miniwdl/Dockerfile
FROM ubuntu

# ----- start of hermit-specific requirements ------
#
# The requirements for a docker image to work with hermitcrab:
# 1. There must be a user named "ubuntu" with uid 2000
# 2. They must be able to sudo without entering a password
# 3. sshd must be installed 

RUN apt-get update && \
  apt-get install -y openssh-server sudo && \
  mkdir /run/sshd && \
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
